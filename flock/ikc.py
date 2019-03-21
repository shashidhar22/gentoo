import os
import sys
import glob
import mmap
import math
import struct
import pandas as pd
from collections import OrderedDict
from collections import namedtuple
from bunting.kmerize import Kmer
from bunting.seq import Seq

class IKC:
    # TODO: Search function for k-mer counts from minimizer bins
    # TODO: Managing k-mer distribution within minimizer bins to ensure uniform distribution of k-mers in each bin
    # TODO: IKC writer module


    def __init__(self, ikc_path):
        ikc_file = open(ikc_path, 'r+b')
        self.ikc_map = mmap.mmap(ikc_file.fileno(), 0)
        self.ikc_name = os.path.basename(ikc_path)
        self.magic = self.ikc_map[0:15]
        self.version = self.ikc_map[15]
        self.reserved = self.ikc_map[16:23]
        self.min_size = self.ikc_map[23]
        self.kmer_size = struct.unpack('>i', self.ikc_map[24:28])[0]
        self.mask = struct.unpack('>i', self.ikc_map[28:32])[0]
        self.index_offset = struct.unpack('>q', self.ikc_map[32:40])[0]
        self.meta_offset = struct.unpack('>q', self.ikc_map[40:48])[0]
        self.minimizer_bins = self.getMinimizerBins()
        return

    def toString(self, kdec):
        kmer = ''
        base_dict = {0x00: 'A', 0x01: 'C', 0x02: 'G', 0x03: 'T'}
        while len(kmer) != self.kmer_size:
            kdec, base = divmod(kdec, 4)
            base = base_dict[base]
            kmer += base
        return(kmer[::-1])

    def findKmer(self, bin, kmer):
        start = bin.start
        stop = bin.stop

    def getMinimizerBins(self):
        bins = OrderedDict()
        for segments in range(self.index_offset, self.meta_offset, 12):
            minimizer = struct.unpack('>i',
                                    self.ikc_map[segments:segments+4])[0]
            offset_start = struct.unpack('>q',
                                    self.ikc_map[segments+4:segments+12])[0]
            if segments + 12 == self.meta_offset:
                offset_end = self.index_offset
            else:
                offset_end = struct.unpack('>q',
                                    self.ikc_map[segments+16:segments+24])[0]
            kmer_size = math.ceil(self.kmer_size/4) + 4
            bins[minimizer] = range(offset_start, offset_end, kmer_size )
        return(bins)

    def streamKmer(self):
        for minimizer in self.minimizer_bins:
            kmer_bin = self.minimizer_bins[minimizer]
            #print(minimizer, kmer_bin)
            for segments in kmer_bin:
                ksize = kmer_bin.step-4
                #print(minimizer, kmer_bin, segments, ksize)
                kmer = int.from_bytes(self.ikc_map[segments: segments+ksize],
                        byteorder='big', signed=True)
                count = struct.unpack('>i',
                        self.ikc_map[segments+ksize: segments+ksize+4])[0]
                kseq = self.toString(kmer)
                yield(kmer, kseq, count)

    def getMeta(self, kdec):
        kmer_obj = Seq(kdec, self.kmer_size)
        kmer_rec = Kmer(kmer_obj, self.kmer_size, self.min_size)
        kmer = kmer_rec.kmer
        minimizer = kmer_rec.minimizer
        try:
            kmer_bin = self.minimizer_bins[minimizer]
            ## TODO: Improve kmer search in k-mer bin, a binary search can be implemented instead of linear search
            for segments in kmer_bin:
                ksize = kmer_bin.step -4
                skmer = int.from_bytes(self.ikc_map[segments: segments+ksize],
                        byteorder='big', signed=True)
                if skmer == kmer:
                    meta_loc = struct.unpack('>i',
                            self.ikc_map[segments+ksize: segments+ksize+4])[0]
                    meta_len = struct.unpack('>i',
                            self.ikc_map[meta_loc: meta_loc + 4])[0]
                    #print(meta_loc, meta_len, self.ikc_map[meta_loc+4: meta_loc+meta_len+4])
                    meta_info = self.ikc_map[meta_loc+4: meta_loc+meta_len+4].decode('utf-8')
                    return(meta_info)
        except KeyError:
            pass
        return(None)

    def getCounts(self, kdec):
        kmer_obj = Seq(kdec, self.kmer_size)
        kmer_rec = Kmer(kmer_obj, self.kmer_size, self.min_size)
        kmer = kmer_rec.kmer
        minimizer = kmer_rec.minimizer
        try:
            kmer_bin = self.minimizer_bins[minimizer]
            count_dict = OrderedDict()
            for segments in kmer_bin:
                ksize = kmer_bin.step - 4
                skmer = int.from_bytes(self.ikc_map[segments: segments+ksize],
                           byteorder='big', signed=True)
                if skmer == kmer:
                    scount = struct.unpack('>i',
                               self.ikc_map[segments+ksize: segments+ksize+4])[0]
                    return(scount)
        except KeyError:
            pass
        return(0)

    def kmersPerBins(self):
        sample = os.path.splitext(self.ikc_name)[0]
        bin_file = open('{0}_{1}.bs'.format(sample, self.min_size), 'w')
        for minimizers in self.minimizer_bins:
            kmer_bins = self.minimizer_bins[minimizers]
            kpb = 0
            for segments in kmer_bins:
                kpb += 1
            bin_file.write('{0}\t{1}\n'.format(minimizers, kpb))
        bin_file.close()
        #count_file.close()
        return


class IKCWriter:

    def __init__(self, kmer_list, count_list, meta_list, out_path, sample_name=None,
        magic='Idx_Kmer_Class', version=2, reserved=None, kmer_size=31, min_size=15,
        min_mask=0):

        self.sample_name = sample_name
        self.magic = magic
        self.version = version
        self.reserved = reserved
        self.kmer_size = kmer_size
        self.out_path = out_path
        self.min_size = min_size
        self.min_mask = min_mask
        self.kmer_list = kmer_list
        self.count_list = count_list
        self.meta_list = meta_list
        self.meta_bins, self.meta_info = self.getMetaBins()

    def getMetaBins(self):
        eclass_file = open('{0}/test.ec'.format(self.out_path), 'w')
        meta_bins = OrderedDict()
        meta_info = OrderedDict()
        meta_rev = OrderedDict()
        index = 1
        for kmer, ec_list in self.meta_list.items():
            eclass_file.write('{2}\t{0}\t{1}\n'.format(kmer, ec_list, index))
            try:
                meta_bins[kmer] = meta_rev[ec_list]
            except KeyError:
                meta_rev[ec_list] = index
                meta_info[index] = ec_list
                meta_bins[kmer] = index
            index += 1
            #if ec_list in meta_rev:
            #    if
            #    continue
            #    meta_info[index] = ec_list
            #    meta_bins[kmer] = index
            #    index += 1
        eclass_file.close()
        return(meta_bins, meta_info)

    def getKmersPerBin(self):
        ## TODO: Fix the kmer issue
        kmerizer = Kmer(self.min_size)
        min_list = list()
        converter = Kmer(self.kmer_size)
        for kmers in kmer_list:
            fwd_seq = converter.toString(kmers)
            rev_seq = converter.revComp(fwd_seq)
            minimizer = kmerizer.getMinimizer(fwd_seq, rev_seq)
            if minimizer in kmers:
                continue
            else:
                min_list.append(minimizer)
        return(sorted(min_list))

    def sortKmers(self):
        sorted_kmers = dict()
        kmerizer = Kmer(self.min_size)
        for i in range(len(self.kmer_list)):
            kmers = self.kmer_list[i]
            count = self.count_list[i]
            converter = Kmer(self.kmer_size)
            fwd_kmer = converter.toString(kmers)
            rev_kmer = converter.revComp(fwd_kmer)
            minimizer = kmerizer.getMinimizer(fwd_kmer, rev_kmer)
            try:
                sorted_kmers['Minimizer'].append(minimizer)
                sorted_kmers['Kmer'].append(kmers)
                sorted_kmers['Count'].append(count)
            except KeyError:
                sorted_kmers['Minimizer'] = [minimizer]
                sorted_kmers['Kmer'] = [kmers]
                sorted_kmers['Count'] = [count]

        sorted_kmers = pd.DataFrame(sorted_kmers).sort_values(by=['Minimizer',
            'Kmer'])
        min_list = sorted_kmers['Minimizer'].value_counts(sort=False)
        kmer_list = sorted_kmers['Kmer'].tolist()
        sorted_kmers.to_csv('{0}/test.kcsv'.format(self.out_path))
        return(min_list, kmer_list)


    def writeIKC(self):
        min_list, kmer_list = self.sortKmers()
        meta_bins, meta_info = self.getMetaBins()
        ikc_file = '{0}/test.ikc'.format(os.path.abspath(self.out_path))
        ikc_handle = open(ikc_file, 'wb')
        print(self.magic)
        ikc_handle.write(struct.pack('>15s', self.magic.encode('utf-8')))
        ikc_handle.write(struct.pack('>b', self.version))
        ikc_handle.write(struct.pack('>7x'))
        ikc_handle.write(struct.pack('>b', self.min_size))
        ikc_handle.write(struct.pack('>i', self.kmer_size))
        ikc_handle.write(struct.pack('>i', self.min_mask))
        data_offset = 80
        index_offset = (len(kmer_list) * (math.ceil(self.kmer_size/4) + 4)) + data_offset
        meta_offset = (len(min_list) * 12) + index_offset
        ikc_handle.write(struct.pack('>q', index_offset))
        ikc_handle.write(struct.pack('>q', meta_offset))
        ikc_handle.write(struct.pack('>32x'))
        index = 80
        kmerize = Kmer(self.min_size)
        kmer_offsets = dict()
        meta_offsets = dict()
        for kmers in kmer_list:
            converter = Kmer(self.kmer_size)
            fwd_kmers = converter.toString(kmers)
            rev_kmers = converter.revComp(fwd_kmers)
            minimizer = kmerize.getMinimizer(fwd_kmers, rev_kmers)
            if minimizer not in kmer_offsets:
                kmer_offsets[minimizer] = index
            bit_length = math.ceil(self.kmer_size/4)
            ikc_handle.write(kmers.to_bytes(bit_length, byteorder='big'))
            index += bit_length
            ikc_handle.write(struct.pack('>q', meta_offset))
            index += 8
            metadata = len(meta_info[meta_bins[kmers]])
            meta_offset += metadata + 4
            meta_offsets[meta_info[meta_bins[kmers]]] = meta_offset

        for minimizer in min_list.index:
            ikc_handle.write(struct.pack('>i', minimizer))
            #index += 4
            ikc_handle.write(struct.pack('>q', kmer_offsets[minimizer]))
            #index += 8

        meta_sorted = sorted(meta_offsets.items(), key=lambda kv: kv[1])
        for ec_list in meta_offsets:
            ec = ec_list
            ec_length = len(ec)
            print(ec_length, ec)
            ikc_handle.write(struct.pack('>i', ec_length))
            ikc_handle.write(struct.pack('>{0}s'.format(ec_length), ec.encode('utf-8')))

        ikc_handle.close()
        return


        #struct.pack_into('>q', ikc_handle, 32, self)


    def printDetails(self):
        print('File name : {0}'.format(self.ikc_name))
        print('Magic : {0}'.format(str(self.magic, 'utf-8')))
        print('File version : {0}'.format(self.version))
        print('Minimizer size : {0}'.format(self.min_size))
        print('K-mer size : {0}'.format(self.kmer_size))
        print('K-mer mask : {0}'.format(self.mask))
        return
