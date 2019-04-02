import os
import sys
import re
import glob
import logging
import subprocess
import numpy as np
from collections import namedtuple
from flock.fastq import Fastq
from itertools import groupby


class Identifier:
    """ The Identifier class contains a set of modules that recognize the
    type of fastq file from the fastq the fastq header. The class is initialzed
    by passing a record to the constructor"""

    def __init__(self, record):
        """Initialize the Identifier class from a fastq record"""
        self.rec = record


    def isIlluminaOld(self):
        """Identify fastq headers in the following format
         @HWUSI-EAS100R:6:73:941:1973#0/1"""
        header_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


    def isIlluminaNew(self):
        """Idetify fastq headers in the following format
        #@D00468:24:H8ELMADXX:1:1101:1470:2237 1:N:0:2"""
        header_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s\d:\w+:\w+:\w*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


    def isSraOld(self):
        """Identify fastq headers in the following format
        #@SRR037455.1 HWI-E4_6_30ACL:4:1:0:29 length=35
        #@SRR902931.1 HWI-ST1384:61:D1DJ4ACXX:8:1101:1240:2015 length=50"""
        header_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isSraNew(self):
        """Identify fastq headers in the following format
        #@SRR1873770.5 DH1DQQN1:437:HACT2ADXX:1:2204:8270:58140 length=150"""
        header_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isENA(self):
        """Identify fastq headers in the following format
        # ERR161234.14 14 length=100"""
        header_regex = re.compile('@[\w\.]+ \d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isPacbio(self):
        """Identify fastq headers in the following format
        #@m160113_152755_42135_c100906712550000001823199104291667_s1_p0/15/7044_26271"""
        header_regex = re.compile('@\w+/\d+/\d+_\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


class Metrics:
    """The Metrics class calculates read length and quality metrics for a given
    fastq file. The class is initialized with the following inputs:
      1. fastq (str): Path to the fastq file"""

    def __init__(self, fastq):
        self.fastq = fastq

    def avgReadLen(self):
        """Given a fastq file, calculate average read length"""
        fastq_reader = Fastq(self.fastq, 'phred33')
        total_length = 0
        total_reads = 0
        for lines in fastq_reader.read():
            total_length += len(lines.seq)
            total_reads += 1
            if total_reads >= 100:
                break

        avg_length = total_length/float(total_reads)
        return(avg_length)

class Prepper:
    """Prepper class create the config dictionary for a given study. The class
    constructor takes the following input parameters:
      1. input_path (str) : Path to input directory or sra accession list
      2. sra_path (str) : Path to fastq-dump executable """

    def __init__(self, input_path, sra_path):
        self.input_path = os.path.abspath(input_path)
        self.sra_path = sra_path
        self.logger = logging.getLogger('Kookaburra.prepInputs')

    def downloadSRA(self):
        """Give a SRA accession list, download all the associated fastq files"""
        out_dir = os.path.dirname(self.input_path)
        sra_list = open(self.input_path)
        for accessions in sra_list:
            self.logger.debug('Downloading : {0}'.format(accessions))
            accessions = accessions.strip()
            fqd_cmd = [self.sra_path, '-f', '-3', '-O', out_dir, accessions]
            fqd_run = subprocess.Popen(fqd_cmd, shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            fqd_run.wait()
            if fqd_run.returncode != 0:
                self.logger.error('Could not download {0}'.format(accessions))
                self.logger.error(' '.join(fqd_cmd))
            else:
                self.logger.info('Downladed complete: {0}'.format(accessions))
        return(out_dir)

    def getFastqPaths(self):
        """Given a directory path, extract all the fastq file names"""
        filenames = list()
        for subdir, dirname, files in os.walk(self.input_path):
            for filename in files:
                if ('.fastq' in filename or '.fastq.gz' in filename or
                    'fq' in filename or 'fq.gz' in filename):
                    filepath = subdir + os.sep + filename
                    filenames.append(filepath)
        return(filenames)

    def getReadPresence(self, file_name, minimum=1):
        """Given a fastq file, check for the presence of a minimum number of
        reads"""
        reader = Fastq(file_name, 'phred33')
        read_number = 0
        for rec in reader.read():
            read_number += 1
            if read_number >= minimum:
                return(True)
            else:
                return(False)

    def parseMaRS(self, file_name):
        """Given a fastq file, check is the file name follows the MaRS regex:
        <Year>{2}<Country>{2}<Site>{2}<Day of treatment>{2}<Treatment>{1}
        <SampleID>{4}<Genus and species>{2}<Strain type>{1}<Markers>{3}<Rep>{1}
        If the regex is found, the module decodes the sample information and
        returns a namedtuple for the sample"""

        mars_regex = ('(?P<Year>[0-9x]{2})(?P<Country>[A-Zx]{2})'
                      '(?P<Site>[A-Zx]{2})(?P<DT>[0-9]{2})'
                      '(?P<Treatment>[A-Kx]{1})(?P<SID>[0-9]{4})'
                      '(?P<GS>[a-zA-Z]{2})(?P<ST>[BFPTSx]{1})'
                      '(?P<Markers>[0-9]{3})(?P<Rep>[0-9]{1})')
        mars_groups = re.match(mars_regex, file_name)
        sample_info  = namedtuple('Sample', ['Year', 'Country', 'Site',
            'TreatmentDay', 'Treatment', 'ID', 'Genus', 'Type', 'Markers',
            'Replicate'])

        if not mars_groups:
            marker_list = np.array(['PfK13', 'PfCRT', 'PfMDR', 'MT', 'CytB',
                'PfDHPS', 'PfDHFR'])
            sample = sample_info(np.nan, 'xx', 'xx', 'xx', 'xx', 'xxxx', 'xx',
                'x', marker_list, 1)
            return(sample)
        if mars_groups.group('Year') == 'xx':
            year = np.nan
        else:
            year = int('20' + mars_groups.group('Year'))
        country = mars_groups.group('Country')
        site = mars_groups.group('Site')
        treatDate = mars_groups.group('DT')
        treatment = mars_groups.group('Treatment')
        sampleID = mars_groups.group('SID')
        genusSpecies = mars_groups.group('GS')
        type_dict = {'B': 'Blood', 'F': 'Filter blood spots', 'P': 'Plasma',
            'T': 'Tissue', 'S': 'Stool', 'x': 'Unknown'}
        type = type_dict[mars_groups.group('ST')]
        marker_list = np.array(['PfK13', 'PfCRT', 'PfMDR', 'MT', 'CytB',
            'PfDHPS', 'PfDHFR'])
        marker_found = list()
        markers = int(mars_groups.group('Markers'))
        while len(marker_found) < 8:
            markers, rem = divmod(markers, 2)
            marker_found.append(rem)

        marker_found = np.array(marker_found[::-1][1:])
        marker_list = np.extract(marker_found, marker_list)
        rep = int(mars_groups.group('Rep'))
        sample = sample_info(year, country, site, treatDate, treatment,
            sampleID, genusSpecies, type, marker_list, rep)
        return(sample)

    def prepInputs(self):
        """Given a input directory path or sra accession list, iterate through
        the list and create a sample record of reach file. Each sample record is
        added to a dictionary with the sample name as the key and sample record
        as the value"""
        if os.path.isfile(self.input_path):
            self.logger.info('Found SRA accession list,'
                            'Will download files from SRA')
            self.input_path = self.downloadSRA()
        files = self.getFastqPaths()
        experiment = dict()
        for fastq in files:
            reader = Fastq(fastq, 'phred33')
            Sample = namedtuple('Sample', ['sample', 'libname', 'library',
            'files', 'prep', 'paired', 'year', 'country', 'site',
            'treatmentDay', 'treatment', 'iD', 'genus', 'type', 'markers',
            'replicate'])
            readPresence = self.getReadPresence(fastq)
            if not readPresence:
                self.logger.warning('Sample doesn\'t contain minimum number of required reads; skipping sample : {0}'.format(fastq))
                continue
            rec = next(reader.read())
            identifier = Identifier(rec)
            metric = Metrics(fastq)
            isIllOld =  identifier.isIlluminaOld()
            isIllNew =  identifier.isIlluminaNew()
            isSraOld = identifier.isSraOld()
            isSraNew = identifier.isSraNew()
            isPac = identifier.isPacbio()
            isENA = identifier.isENA()
            seqType = ''
            libType = ''
            sample_regex = re.compile('_r1|_r2|_?l001|_?l002|_?l003|_?l004|_R1|_R2|_L001|_?L002|_L003|_L004|_1|_2') #|L001|L002|L003|L004')
            sample = sample_regex.split(os.path.basename(fastq))[0]
            sample_info = self.parseMaRS(sample)
            year = sample_info.Year
            country = sample_info.Country
            site = sample_info.Site
            td = sample_info.TreatmentDay
            treatment = sample_info.Treatment
            sid = sample_info.ID
            gs = sample_info.Genus
            stype = sample_info.Type
            markers = sample_info.Markers
            replicate = sample_info.Replicate
            if isIllOld:
                paired_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isIllNew:
                paired_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isSraOld:
                paired_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isSraNew:
                paired_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isENA:
                paired_regex = re.compile('@[\w\.]+ \d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isPac:
                lib_regex = re.compile('@\w+_\d+_\d+_\w+')
                lib = re.findall(lib_regex, rec.header)[0]
                paired = False
                seqType = 'Pacbio'
                if metric.avgReadLen():
                    libType = 'Long'
            else:
                print(fastq)
                self.logger.warning('Read from {0} with header : {1} does not follow any defined fastq header format.Please correct it'.format(fastq, rec_header))
            try:
                paired = True
                experiment[sample] = Sample(sample, lib, seqType,
                    [experiment[sample].files[0],fastq], libType, paired,
                    year, country, site, td, treatment, sid, gs, stype, markers,
                    replicate)
            except (KeyError, AttributeError):
                experiment[sample] = Sample(sample, lib, seqType, [fastq],
                libType, paired, year, country, site, td, treatment, sid, gs,
                stype, markers, replicate)
        self.logger.info('A total of {0} libraries were identified from the given folder {1}'.format(len(experiment), self.input_path))
        self.logger.debug('The following libraries were detected in the given folder : {0}'.format(self.input_path))
        for sample, values in experiment.items():
            self.logger.debug('Sample : {0}; Library: {1} ; Sequence type: {2} ; Files: {3} ; Library type: {4} ; Paired: {5}'.format(
                    values.sample, values.libname, values.library, ''.join(values.files), values.prep, values.paired))
        for samples, info in experiment.items():
            if not info.paired:
                self.logger.warning('NeST does not currently support single end runs; skipping sample : {0}'.format(samples))
                experiment.pop(samples)
        return(experiment)

if __name__ == '__main__':
    path = os.path.abspath(sys.argv[1])
    prepper = Prepper(path)
    experiment = prepper.prepInputs()
    rone = list()
    rtwo = list()
    for study in experiment:
        for files in experiment[study].files:
            if re.findall('.*_R1.*', files):
                rone.append(files)
            else:
                rtwo.append(files)
    print(rone)
    print(rtwo)
