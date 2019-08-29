import os 
import re
import sys
import glob
import math
import time
import shutil
import logging
import argparse
import unittest
import itertools 
import subprocess
import numpy as np 
import pandas as pd 
from skbio.tree import nj
from itertools import repeat
from skbio import DistanceMatrix
from multiprocessing import Pool
from collections import namedtuple
from flock.fasta import Fasta
from flock.fastq import Fastq

def sraDownload(arguments):
    """Given a sra accession number, download file using fastq-dump and return returncode.

    Parameter list:
        sra = SRA accession number
    
    Return value:
        returncode if fastq-dump ran successfully, else will raise SystemExit exception
    """
    sra = arguments[0]
    out_path = arguments[1]
    logger = logging.getLogger('Gentooo.sra')
    logger.debug('Downloading : {0}'.format(sra))
    out_dir = '{0}/Fastq'.format(out_path)
    out_files = ['{0}/{1}_1.fastq.gz'.format(out_dir, sra), '{0}/{1}_2.fastq.gz'.format(out_dir, sra)]
    if os.path.exists(out_files[0]) and os.path.exists(out_files[1]):
       logger.info('Samples exists {0}'.format(sra))
       return(0)
    fqd_cmd = ['fastq-dump', '--gzip', '--split-3', '-O', out_dir, sra]
    fqd_run = subprocess.Popen(fqd_cmd, shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    fqd_run.wait()
    if fqd_run.returncode != 0:
        #self.log.error(' '.join(fqd_cmd))
        raise SystemExit('Could not download {0}, check accession number and try again'.format(sra))
    else:
        logger.info('Downladed complete: {0}'.format(sra))
    return(fqd_run.returncode)

def runKanalyze(arguments):
    sample = arguments[0]
    files = arguments[1]
    out_path = arguments[2]
    temp_path = arguments[3]
    ksize = arguments[4]
    out_file = '{0}/Dkc/{1}.dkc'.format(out_path, sample)
    temp_loc = '{0}/Dkc/{1}_tmp'.format(temp_path, sample)
    logger = logging.getLogger('Gentoo.kanalyze')
    if not os.path.exists(temp_loc):
        os.mkdir(temp_loc)
    fq_ext = ['fq', 'fq.gz', 'fastq', 'fastq.gz']
    if files[0].split('.', 1)[1] in fq_ext:
        kcmd = ['java', '-jar', 'lib/kanalyze/kanalyze.jar', 'count', '-t', '2', '-m', 'dec', '-k', str(ksize), 
                '-c', 'kmercount:5', '-rcanonical', '--seqfilter' , 'sanger:20', 
                '--temploc', temp_loc, '-o', out_file] +  files 
    else:
        kcmd = ['java', '-jar', 'lib/kanalyze/kanalyze.jar', 'count', '-t', '2', '-m', 'dec', '-k', str(ksize), 
                '--temploc', temp_loc, '-o', out_file] + files
    krun = subprocess.Popen(kcmd, stderr=subprocess.DEVNULL, 
                            stdout=subprocess.DEVNULL, shell=False)
    krun.wait()
    if krun.returncode != 0:
        logger.error('Failed to index {0}, check input files and try again'.format(sample))
        logger.error(' '.join(kcmd))
        return(out_file, 1)
    else:
        shutil.rmtree(temp_loc)
        logger.info('K-mer index generated for {0}'.format(sample))
    return(out_file, krun.returncode)

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
        ext_list=['fa', 'fna', 'fas', 'fasta', 'fq', 'fastq', 'fq.gz', 'fastq.gz']
        if type(self.inp_path) is list:
            #Group files together, handle R1 and R2
            filenames = self.inp_path
            study = self.groupFiles(sorted(filenames))
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
                study = self.groupFiles(sorted(filenames))
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
                sample_regex = re.compile('_r1|_r2|_?l001|_?l002|_?l003|_?l004|_R1|_R2|_L001|_L002|_L003|_L004|_1|_2') #|L001|L002|L003|L004')
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
        sra_runs = pool.map(sraDownload, zip(sra_list, repeat(self.out_path)))
        for sample, file_list in study.items():
            final_list = list()
            for files in file_list:
                if os.path.exists(files):
                    final_list.append(files)
            study[sample] = final_list
        return(study)

    def createIndex(self):
        study = self.prep()
        samples = list(study.keys())
        files = list(study.values())
        threads = self.threads //2 if self.threads > 1 else 1
        pool = Pool(threads)
        index_files = pool.map(runKanalyze, zip(samples, files, repeat(self.out_path), repeat(self.temp_path), repeat(self.kmer)))
        return(index_files)

def stream(kmer_file):
    kmer_reader = open(kmer_file)
    for lines in kmer_reader:
        lines = lines.strip().split('\t')
        kmer = int(lines[0])
        count = int(lines[1])
        yield(kmer, count)

def merge(file_list):
    fone = file_list[0]
    ftwo = file_list[1]
    oname = os.path.splitext(os.path.basename(fone))[0]
    tname = os.path.splitext(os.path.basename(ftwo))[0]
    merge_logger = logging.getLogger('Gentoo')
    #Creating logger for nest
    merge_logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # Create formatter and add it to the handlers
    formatter = logging.Formatter('{asctime} - {name} - {levelname} - {message}', style="{")
    ch.setFormatter(formatter)
    # Add the handlers to the logger
    merge_logger.addHandler(ch)
    ostream = stream(fone)
    tstream = stream(ftwo)
    omer = next(ostream, None)
    tmer = next(tstream, None)
    #Change jaccard to similarity ratio that accounts for abundance
    #Chaning to jaccard index
    intersection = 0
    union = 0
    while omer and tmer:
        if omer[0] == tmer[0]:
            intersection += min(omer[1], tmer[1])
            union += max(omer[1], tmer[1])
            omer = next(ostream, None)
            tmer = next(tstream, None)
        elif omer[0] < tmer[0]:
            union += omer[1]
            omer = next(ostream, None)
        elif omer[0] > tmer[0]:
            union += tmer[1]
            tmer = next(tstream, None)
    while omer:
        union += omer[1]
        omer = next(ostream, None)
    while tmer:
        union += tmer[1]
        tmer = next(tstream, None)
    ## Since a jaccard index depends on just binary presence or absence scenarios,
    ## the similarity of samples can be determine by using a similarity ratio
    ## given by
    ## SRij =  kykiykj / ( kyki2 +  kykj2 -  kykiykj), where
    ##  yki = abundance of kth species in quadrat i
    similarity = intersection/union
    distance = 1 - similarity
    #print('Distance between {0} and {1} = {2}'.format(
    #                oname, tname, distance))
    return(oname, tname, distance)

class Cluster:

    def __init__(self, dkc_path, out_path, kmer=31, threads=4, study=None):
        self.dkc_path = os.path.abspath(dkc_path)
        self.out_path = os.path.abspath(out_path)
        self.kmer = kmer 
        self.threads = threads
        if study is None:
            self.study = 'Gentoo_{0}'.format(int(time.time()))
        else:
            self.study = study
        self.log = logging.getLogger('Gentoo.Cluster')

    def splitter(self):
        file_list = glob.glob('{0}/*.dkc'.format(self.dkc_path))
        products = list()
        combinations = list(itertools.combinations_with_replacement(file_list, 2))
        self.log.debug('Performing {0} pairwise comparisons'.format(len(combinations)))
        print('Performing {0} pairwise comparisons'.format(len(combinations)))
        pools = Pool(self.threads)
        jaccard_list = pools.map(merge, combinations)
        print('Pairwise comparisons completed')
        for groups in jaccard_list:
            products.append(list(groups))
            if [groups[1], groups[0], groups[2]] not in products:
                products.append([groups[1], groups[0], groups[2]])
        return(products)

    def pairwiseToDist(self, outgroup=None):
        jaccard_list = self.splitter()
        dist_file = '{0}/{1}.dist'.format(self.out_path, self.study)
        newick_file = '{0}/{1}.nwk'.format(self.out_path, self.study)
        jaccard_table = pd.DataFrame(jaccard_list, columns=['SampleA', 'SampleB', 'Distance'])
        jaccard_table.to_csv(dist_file, header=True, index=False)
        jaccard_matrix = jaccard_table.pivot(index='SampleA', columns='SampleB', values='Distance')
        samples = jaccard_matrix.index
        jaccard_matrix = jaccard_matrix.to_numpy()
        #jaccard_matrix = np.asmatrix(jaccard_matrix)
        jaccard_matrix = DistanceMatrix(jaccard_matrix, samples)
        newick_str = nj(jaccard_matrix, result_constructor=str)
        newick_file = open(newick_file, 'w')
        newick_file.write('{0}\n'.format(newick_str))
        newick_file.close()


if __name__ == '__main__':

    #Get arguments
    parsers = argparse.ArgumentParser(prog='Gentoo')
    parsers.add_argument('-i', '--inp_path', type=str,
                         help='Path to input file or directory')
    parsers.add_argument('-o', '--out_path', type=str,
                         help='Path to output directory')
    parsers.add_argument('-t', '--tmp_path', type=str, default=None,
                         help='Path to temp directory')
    parsers.add_argument('-p', '--threads', type=int, default=5,
                         help='Number of threads used for analysis')
    parsers.add_argument('-k', '--ksize', type=int, default=31,
                         help='K-mer size used for analysis')
    parsers.add_argument('-s', '--study', type=str, default=None,
                         help='Study name')
    parsers.add_argument('--index', action='store_true', 
                         help='Create KC index files')
    parsers.add_argument('--cluster', action='store_true',
                         help='Create distance matrix from indexed files')
    args = parsers.parse_args()

    #Index code
    if args.index:
       indexer = Index(args.inp_path, args.out_path, args.threads, args.ksize, args.tmp_path)
       index_files = indexer.createIndex()

    if args.cluster:
       cluster = Cluster(args.inp_path, args.out_path, args.ksize, args.threads, args.study)
       cluster.pairwiseToDist()
