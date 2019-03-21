import os
import sys
import glob
import heapq
import logging
import subprocess
import numpy as np
from bunting.seq import Seq
from bunting.fasta import Fasta
from bunting.fastq import Fastq
from collections import namedtuple

class Kmer:

    def __init__(self, seq_obj, ksize=21, msize=15):
        self.seq_obj = seq_obj
        self.ksize = ksize
        self.msize = msize


    def kmerize(self, seq, ksize):
        mask = (1 << (ksize*2)) -1
        frag = 0
        bit_counter = 1
        for nuc in seq:
            frag  = frag << 2
            if nuc == 'a' or nuc == 'A':
                frag = frag | 0x00
            elif nuc == 'c' or nuc == 'C':
                frag = frag | 0x01
            elif nuc == 'g' or nuc == 'G':
                frag = frag | 0x02
            elif nuc == 't' or nuc == 'T' or nuc == 'u' or nuc == 'U':
                frag = frag | 0x03
            else:
                bit_counter = 0

            if bit_counter == ksize:
                kmer = frag & mask
                kseq = self.toString(kmer)
                krev = self.revComp(kseq)
                yield(kmer, kseq, krev)
            else:
                bit_counter += 1

    def getMinimizer(self, fwd_seq, rev_seq):
        fwd_kmers = set([kmer for kmer, kseq, krev in self.kmerize(fwd_seq, self.msize)])
        rev_kmers = set([kmer for kmer, kseq, krev in self.kmerize(rev_seq, self.msize)])
        minimizer = min(fwd_kmers | rev_kmers)
        return(minimizer)

    def revComp(self, seq):
        rev_seq = ''
        rev_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        for base in seq[::-1]:
            rev_seq += rev_dict[base]
        return(rev_seq)

    def toString(self, kdec):
        kmer = ''
        base_dict = {0x00: 'A', 0x01: 'C', 0x02: 'G', 0x03: 'T'}
        while len(kmer) != self.ksize:
            kdec, base = divmod(kdec, 4)
            base = base_dict[base]
            kmer += base
        return(kmer[::-1])

    def stream(self):
        for sequences in self.seq_obj.read:
            seq = sequences.seq
            header = sequences.header
            Kmer = namedtuple('Kmer', ['kmer', 'kseq', 'rseq', 'minimizer',
                    'origin', 'ksize', 'msize', 'header'])
            for kmer, kseq, krev in self.kmerize(seq, self.ksize):
                minimizer = self.getMinimizer(kseq, krev)
                krec = Kmer(kmer, kseq, krev, minimizer, header, self.ksize,
                        self.msize, header)
                yield(krec)

    def heapSort(self, items):
        heapq.heapify(items)
        items[:] = [heapq.heappop(items) for i in range(len(items))]
        return(items)

    def counter(self, items):
        count_list = list()
        kmer_list = list()
        prev = None
        for pos, item in enumerate(items, start=1):
            if item != prev:
               count_list.append(1)
               kmer_list.append(item)
               prev = item
            elif item == prev:
               count_list[-1] += 1
        return(kmer_list, count_list)

    def pyKanalyze(self, seq_obj=None):
        if seq_obj is None:
            seq_obj = self.seq_obj
        kmer_list = list()
        #logging.warning('Starting k-mer counting')
        #t1 = time.time()
        for records in seq_obj.read:
            for kmer, kseq, krev in self.kmerize(records.seq, self.ksize):
                kmer_list.append(kmer)
        #t2 = time.time()
        #logging.warning('Found all k-mers, sorting k-mers')
        kmer_list = self.heapSort(kmer_list)
        #t3 = time.time()
        #   logging.warning('{0} kmers sorted'.format(len(kmer_list)))
        kmer_list, count_list = self.counter(kmer_list)


        return(kmer_list, count_list)

    def writeCounts(self, kmer_list, count_list, out_path):
        if os.path.isdir(out_path):
            out_path = '{0}/count.pkc'.format(out_path)
        count_file = open(out_path, 'w', buffering=1048576000)
        for i in range(len(kmer_list)):
            count_file.write('{0}\t{1}\n'.format(kmer_list[i], count_list[i]))
        count_file.close()
        return(out_path)

    def getFastqPaths(self, fastq_path):
        prepper = Prepper(fastq_path, None)
        filenames = prepper.prepInputs()
        return(filenames)

    def kanalyze(self, fastq_path, out_file, format='fastq'):

        if format == 'fastq':
            fastq_files = self.getFastqPaths(fastq_path)
            for fastq in fastq_files:
                kcmd = ['count', '-k', '{0}'.format(self.ksize), '-m', 'ikc',
                    '-o', out_file, '-rcanonical'] + fastq.files
                krun = subprocess.Popen(kcmd, stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL, shell=False)
                krun.wait()
                if krun.returncode != 0:
                    print('KAnalyze exited with return code : {0}'.format(
                            krun.returncode))
                    print('Run command : \n {0}'.format(' '.join(kcmd)))
            return
        else:
            kcmd = ['count', '-k', '{0}'.format(self.ksize), '-m', 'ikc', '-o',
                out_file, fastq_path]
            krun = subprocess.Popen(kcmd, stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL, shell=False)
            krun.wait()

        if krun.returncode != 0:
            print('KAnalyze exited with return code : {0}'.format(
                    krun.returncode))
            print('Run command : \n {0}'.format(' '.join(kcmd)))
            return(None)
        else:
            return(out_file)

    def kanalyzeByRead(self, out_path):
        out_path = os.path.abspath(out_path)
        out_list = list()
        out_file = '{0}/read_kmers.kc'.format(out_path)
        kc_handle = open(out_file, 'w')
        for records in self.seq_obj.read:
            kmers = self.kmerize(records.seq, self.ksize)
            if records.type == 'fasta':
                kc_handle.write('>{0}\n'.format(records.header))
            else:
                kc_handle.write('{0}\n'.format(records.header))
            for kmer, kseq, krev in kmers:
                kc_handle.write('{0}\n'.format(kmer))
        kc_handle.close()

    def streamByRecord(self):
        for records in self.seq_obj.read:
            kmers = self.kmerize(records.seq, self.size)
            for kmer_rec in kmers:
                yield kmer_rec

    def countByRead(self, out_path):
        out_path = os.path.abspath(out_path)
        out_list = list()
        out_file = '{0}/readCount_kmers.kc'.format(out_path)
        kc_handle = open(out_file, 'w')
        for records in self.seq_obj.read:
            seq = Seq(records.seq)
            kmer_list, count_list = self.pyKanalyze(seq)
            kc_handle.write('{0}\n'.format(records.header))
            for kmer, count in zip(kmer_list, count_list):
                kc_handle.write('{0}\t{1}\n'.format(kmer,count))
        kc_handle.close()


    def kanalyzeByContig(self, out_path):
        out_path = os.path.abspath(out_path)
        out_list = list()
        out_file = open('{0}/fasta_kmers.kc'.format(out_path), 'w')
        for records in self.seq_obj.read:
            seq = Seq(records.seq)
            header = records.header
            kmer_list, count_list = self.pyKanalyze(seq)
            out_file.write('>{0}\n'.format(records.header))
            for kmer, count in zip(kmer_list, count_list):
                out_file.write('{0}\t{1}\n'.format(kmer, count))
        out_file.close()
        return('{0}/fasta_kmers.kc'.format(out_path))
