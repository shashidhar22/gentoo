import os 
import re
import sys
import glob
import math
import logging
import unittest 
import subprocess
import numpy as np 
import pandas as pd 
from multiprocessing import Pool
from collections import namedtuple
from flock.fasta import Fasta
from flock.fastq import Fastq


class Index:

    def __init__(self, inp_path, out_path, threads=1, kmer=31, temp_path=None):
        self.log = logging.getLogger('Gentoo.Index')
        if type(inp_path) is str:
            self.inp_path = os.path.abspath(inp_path)
        else:
            self.inp_path = inp_path
        self.out_path = os.path.abspath(out_path)
        if temp_path is None:
            self.temp_path = self.out_path
        else:
            self.temp_path = temp_path
        if not os.path.exists('{0}/Fastq'.format(self.out_path)):
            os.mkdir('{0}/Fastq'.format(self.out_path))
        if not os.path.exists('{0}/Dkc'.format(self.out_path)):
            os.mkdir('{0}/Dkc'.format(self.out_path))
        self.kmer = kmer
        self.threads = threads
        #self.study = self.prep()

    def prep(self):
        """The prep module prepares all the samples for indexing. It will 
        accept the following formats of input:
        1. Study file: A tab delimited file with sample name, and comma 
                       separated list of files associated with the samples.
        2. SRA accession list: A list of accession numbers, sample name 
                               will default to accession number.
        3. Input path: Path to the directory with all the files, sample
                       will default to file names
        4. File list: Space separated list or wildcard string provide 
                      with input command

        Parameter list:
            None
            
        Return value:
            Dictionary with sample name as key and file list as value
        """
        ext_list=['fa', 'fas', 'fasta', 'fq', 'fastq', 'fq.gz', 'fastq.gz']
        if type(self.inp_path) is list:
            #Group files together, handle R1 and R2
            filenames = self.inp_path
            study = self.groupFiles(filenames)
        elif type(self.inp_path) is str:
            if os.path.isfile(self.inp_path):
                #The file is either a SRA accession list
                #Or a tab delimited file with samples
                study = self.getStudy()
            else:
                #Glob files from input path and 
                #Group files togetherm handle R1 and R2
                filenames = list()
                for subdir, dirname, files in os.walk(self.inp_path):
                    for filename in files:
                        ext = filename.split('.', 1)[1]
                        if ext in ext_list:
                            filepath = subdir + os.sep + filename
                            filenames.append(filepath)
                        else:
                            #log a warning that file is being skipped
                            self.log.warning('{0} is not a valid FASTQ/FASTA file'.format(filename))
                study = self.groupFiles(filenames)
        else:
            #Generate expection as inputs dont match accepted format
            raise SystemExit('{0} is not a valid FASTQ/FASTA file'.format(filename))
        #Return study which is a dictionary with sample name
        # as key and file list as values
        return(study)

    def groupFiles(self, file_list):
        """groupFiles takes a file list as input and groups them into samples.
        Parameter list:
            file_list = Unordered list of files for indexing. Can be FASTQ or 
                        FASTA files
        Return value:
            study = Dictionary with sample name as key and file list as value.
        """
        study = dict()
        fasta_ext = ['fa', 'fas', 'fasta']
        fastq_ext = ['fq', 'fq.gz', 'fastq', 'fastq.gz']
        for files in file_list:
            filename = os.path.basename(files)
            ext = filename.split('.', 1)[1]
            if ext in fasta_ext:
                sample = filename.split('.', 1)[0]
                study[sample] = [files]
            elif ext in fastq_ext:
                sample_regex = re.compile('_r1|_r2|_?l001|_?l002|_?l003|_?l004|_R1|_R2|_L001|_?L002|_L003|_L004|_1|_2') #|L001|L002|L003|L004')
                basename = filename.split('.', 1)[0]
                sample = sample_regex.split(basename)[0]
                try:
                    study[sample].append(files)
                except KeyError:
                    study[sample] = [files]
        return(study)

    def getStudy(self):
        """Given a SRA accession list or Study file as input, return study 
        dictionary for all samples

        Parameter list:
            file_name = Path to SRA accession file or Study file
        Return value:
            study = Dictionary with sample name as key and file list as value.
        """
        study = dict()
        study_file = open(self.inp_path)
        sra_list = list()
        for lines in study_file:
            lines = lines.strip().split('\t')
            sample = lines[0]
            file_list = lines[1].split(',')
            for files in file_list:
                if os.path.isfile(os.path.abspath(files)):
                    files = os.path.abspath(files)
                    try:
                        study[sample].append(files)
                    except KeyError:
                        study[sample] = [files]
                else:
                    if not os.path.exists('{0}/Fastq'.format(self.out_path)):
                        os.mkdir('{0}/Fastq'.format(self.out_path))
                    rone = '{0}/Fastq/{1}_1.fastq.gz'.format(self.out_path, files)
                    rtwo = '{0}/Fastq/{1}_2.fastq.gz'.format(self.out_path, files)
                    rrun = '{0}/Fastq/{1}.fastq.gz'.format(self.out_path, files)
                    sra_list.append(files)
                    try:
                        study[sample].extend((rone, rtwo))
                    except KeyError:
                        study[sample] = [rone, rtwo]
        pool = Pool(self.threads)
        sra_runs = pool.map(self.sraDownload, sra_list)
        for sample, file_list in study.items():
            final_list = list()
            for files in file_list:
                if os.path.exists(files):
                    final_list.append(files)
            study[sample] = final_list
        return(study)

    def sraDownload(self, sra):
        """Given a sra accession number, download file using fastq-dump and return returncode.

        Parameter list:
            sra = SRA accession number
        
        Return value:
            returncode if fastq-dump ran successfully, else will raise SystemExit exception
        """
        self.log.debug('Downloading : {0}'.format(sra))
        out_dir = '{0}/Fastq'.format(self.out_path)
        fqd_cmd = ['fastq-dump', '--gzip', '--split-3', '-O', out_dir, sra]
        fqd_run = subprocess.Popen(fqd_cmd, shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        fqd_run.wait()
        if fqd_run.returncode != 0:
            #self.log.error(' '.join(fqd_cmd))
            raise SystemExit('Could not download {0}, check accession number and try again'.format(sra))
        else:
            self.log.info('Downladed complete: {0}'.format(sra))
        return(fqd_run.returncode)

    def runKanalyze(self, arguments):
        sample = arguments[0]
        files = arguments[1]
        out_file = '{0}/Dkc/{1}.dkc'.format(self.out_path, sample)
        temp_loc = '{0}/Dkc/{1}_tmp'.format(self.temp_path, sample)
        if not os.path.exists(temp_loc):
            os.mkdir(temp_loc)
        fq_ext = ['fq', 'fq.gz', 'fastq', 'fastq.gz']
        if files[0].split('.', 1) in fq_ext:
            kcmd = ['lib/kanalyze/count', '-t', '2', '-m', 'dec', '-k', str(self.kmer), 
                    '-c', 'kmercount:2', '-rcanonical', '--seqfilter' , 'sanger:20', 
                    '--temploc', temp_loc, '-o', out_file] +  files 
        else:
            kcmd = ['lib/kanalyze/count', '-t', '2', '-m', 'dec', '-k', str(self.kmer), 
                    '--temploc', temp_loc, '-o', out_file] + files
        krun = subprocess.Popen(kcmd, stderr=subprocess.DEVNULL, 
                                stdout=subprocess.DEVNULL, shell=False)
        krun.wait()
        if krun.returncode != 0:
            self.log.error('Failed to index {0}, check input files and try again'.format(sample))
            self.log.error(' '.join(kcmd))
            return(out_file, 1)
        else:
            self.log.info('K-mer index generated for {0}'.format(sample))
        return(out_file, krun.returncode)

    def createIndex(self):
        study = self.prep()
        samples = list(study.keys())
        files = list(study.values())
        threads = self.threads //2 if self.threads > 1 else 1
        pool = Pool(threads)
        index_files = pool.map(self.runKanalyze, zip(samples, files))
        return(index_files)

