import os
import sys
import glob
import math
import logging
import itertools
import subprocess
import pandas as pd
import numpy as np
from pprint import pprint
from multiprocessing import Pool
from flock.prepinputs import Prepper
from flock.fasta import Fasta
from flock.ikc import IKC
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree, faces, TreeStyle

def fqindex(sam_path):
    sam_files = Prepper('{0}/fastq'.format(sam_path), 'fasterq-dump').prepInputs()
    out_files = list()
    if not os.path.exists('{0}/dkc'.format(ref_path)):
        os.mkdir('{0}/dkc'.format(ref_path))
    for files, sample in sam_files.items():
        out_file = '{0}/dkc/{1}_{2}.dkc'.format(sam_path, 
                                                sample.sample, sample.country)
        out_files.append(out_file)
        kcmd = ['lib/kanalyze/count', '-t', '4', '-m', 'dec', '-k', '31', '-c', 
                'kmercount:2', '-rcanonical', '--seqfilter' , 'sanger:20', 
                '-o', out_file] +  sample.files
        krun = subprocess.Popen(kcmd, shell=False)
        krun.wait()
    return(out_files)

def index(ref_path):
    ref_files = glob.glob('{0}/fasta/*.fasta'.format(ref_path))
    out_files = list()
    if not os.path.exists('{0}/dkc'.format(ref_path)):
        os.mkdir('{0}/dkc'.format(ref_path))
    for files in ref_files:
        out_file = '{0}/dkc/{1}.dkc'.format(ref_path, os.path.splitext(os.path.basename(files))[0])
        out_files.append(out_file)
        kcmd = ['lib/kanalyze/count', '-t', '4', '-m', 'dec', '-k', '31', '-o', out_file,
                files]
        krun = subprocess.Popen(kcmd, shell=False)
        krun.wait()
    return(out_files)

def stream(kmer_file):
    kmer_reader = open(kmer_file, buffering=13107200)
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
    merge_logger = logging.getLogger('PyFinch.{0}.{1}'.format(oname, tname))
    ostream = stream(fone)
    tstream = stream(ftwo)
    omer = next(ostream, None)
    tmer = next(tstream, None)
    #Change jaccard to similarity ratio that accounts for abundance
    #Chaning to jaccard index
    #fones = 0
    #ftwos = 0
    #fonetwos = 0
    intersection = 0
    union = 0
    while omer and tmer:
        if omer[0] == tmer[0]:
            #fones += omer[1]**2
            #ftwos += tmer[1]**2
            #fonetwos += omer[1] * tmer[1]
            intersection += min(omer[1], tmer[1])
            union += max(omer[1], tmer[1])
            omer = next(ostream, None)
            tmer = next(tstream, None)
        elif omer[0] < tmer[0]:
            #fones += omer[1] **2
            #ftwos += 0 **2
            #fonetwos += omer[1] * 0
            union += omer[1]
            omer = next(ostream, None)
        elif omer[0] > tmer[0]:
            #fones += 0 **2
            #ftwos += tmer[1]**2
            #fonetwos += 0 * tmer[1]
            union += tmer[1]
            tmer = next(tstream, None)

    while omer:
        #fones += omer[1] **2
        #ftwos += 0 **2
        #fonetwos += omer[1] * 0
        #print(omer, fone, ftwo)
        union += omer[1]
        omer = next(ostream, None)

    while tmer:
        #fones += 0 **2
        #ftwos += tmer[1]**2
        #fonetwos += 0 * tmer[1]
        union += tmer[1]
        tmer = next(tstream, None)
    ## Since a jaccard index depends on just binary presence or absence scenarios,
    ## the similarity of samples can be determine by using a similarity ratio
    ## given by
    ## SRij =  kykiykj / ( kyki2 +  kykj2 -  kykiykj), where
    ##  yki = abundance of kth species in quadrat i
    #similarity = fonetwos / (fones + ftwos - fonetwos)
    #distance = 1 - similarity
    similarity = intersection/union
    distance = 1 - similarity
    merge_logger.debug('Distance between {0} and {1} = {2}'.format(
                       oname, tname, distance))
    return(oname, tname, distance)

def splitter(ref_path):
    file_list = glob.glob('{0}/dkc/*.dkc'.format(ref_path))
    products = list()
    combinations = list(itertools.combinations_with_replacement(file_list, 2))
    splitter_logger = logging.getLogger('PyFinch')
    splitter_logger.debug('Performing {0} pairwise comparisons'.format(len(combinations)))
    pools = Pool(4)
    jaccard_list = pools.map(merge, combinations)
    for groups in jaccard_list:
        products.append(list(groups))
        if [groups[1], groups[0], groups[2]] not in products:
            products.append([groups[1], groups[0], groups[2]])
    return(products)

def pairwise_to_distance(jaccard_list):
    jaccard_table = pd.DataFrame(jaccard_list, columns=['SampleA', 'SampleB', 'Distance'])
    jaccard_table.to_csv('Plasmodium.dist', header=True, index=False)
    jaccard_matrix = jaccard_table.pivot(index='SampleA', columns='SampleB', values='Distance')
    #jaccard_matrix = pd.pivot_table(jaccard_table, values='Distance', index=['SampleA', 'SampleB'], fill_value=1.0)
    #print(jaccard_matrix.shape)
    samples = jaccard_matrix.index
    jaccard_matrix = jaccard_matrix.to_numpy()
    jaccard_matrix = DistanceMatrix(jaccard_matrix, samples)
    newick_str = nj(jaccard_matrix, result_constructor=str)
    print(newick_str)
    tree = Tree(newick_str)
    tree.set_outgroup(tree&"PlasmoDB-41_Pgallinaceum8A_Genome")
    print(tree)
    #np.savetxt('Plasmodium.dist', jaccard_matrix, delimiter=',')
    #print(jaccard_matrix.shape)

if __name__ == '__main__':
    ref_path = os.path.abspath(sys.argv[1])
    #fone_path = sys.argv[2]
    #ftwo_path = sys.argv[3]
    #Creating logger for nest
    logger = logging.getLogger('PyFinch')
    logger.setLevel(logging.DEBUG)
    # Creating a console handler to log info messages
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # Create formatter and add it to the handlers
    formatter = logging.Formatter('{asctime} - {name} - {levelname} - {message}', style="{")
    ch.setFormatter(formatter)
    # Add the handlers to the logger
    logger.addHandler(ch)
    #out_file = fqindex(ref_path)
    #out_file = index(ref_path)
    #intersection, union, jaccard  = merge((fone_path, ftwo_path))
    jaccard_list = splitter(ref_path)
    pairwise_to_distance(jaccard_list)
    #print(intersection, union, jaccard)
