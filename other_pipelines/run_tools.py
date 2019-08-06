'''The tools executes four different methods to calculate ANI or 
pairwise distance. 
Author: Shashidhar Ravishankar
Date: 07/15/2019
'''
import os
import glob
import time
import shutil
import argparse
import subprocess
from itertools import repeat
from multiprocessing import Pool

def run_pyani(input_path, out_path):
    '''Run PyANI under Blastn mode to calculate ANIb'''
    start = time.time()
    pycmd = ['average_nucleotide_identity.py', '--workers', '5', '-i',
             input_path, '-o', out_path, '-m', 'ANIm', '-g', '-f', '-l',
             '{0}/PyANI.log'.format(out_path)]
    pyrun = subprocess.Popen(pycmd, shell=False)
    pyrun.wait()
    stop = time.time() - start
    print('Mode: PyAni; Input: {0}; Runtime: {1}'.format(input_path, stop))

def run_fastani(input_path, out_path):
    '''Run fastANI to calculate probabilistic ANI'''
    start = time.time()
    files = glob.glob('{0}/*.fasta'.format(input_path))
    ref_list = '{0}/Reference_list.txt'.format(out_path)
    qry_list = '{0}/Query_list.txt'.format(out_path)
    out_file = '{0}/fastANI_output.tsv'.format(out_path)
    ref_file = open(ref_list, 'w')
    qry_file = open(qry_list, 'w')
    for paths in files:
        ref_file.write('{0}\n'.format(os.path.abspath(paths)))
        qry_file.write('{0}\n'.format(os.path.abspath(paths)))
    ref_file.close()
    qry_file.close()
    facmd = ['fastANI', '-t', '10', '--ql', qry_list, '--rl', ref_list,
             '--matrix', '-o', out_file]
    farun = subprocess.Popen(facmd, shell=False)
    farun.wait()
    stop = time.time() - start
    print('Mode: fastANI; Input: {0}; Runtime: {1}'.format(input_path, stop))
 
def run_fastqni(input_path, out_path):
    '''Run fastANI to calculate probabilistic ANI'''
    start = time.time()
    files_one = sorted(glob.glob('{0}/Fastq/*_1.fastq'.format(input_path)))
    files_two = sorted(glob.glob('{0}/Fastq/*_2.fastq'.format(input_path)))
    fq_out = '{0}/MergedFq'.format(out_path)
    if not os.path.exists(fq_out):
        os.mkdir(fq_out)
    merged_files = list()
    for rone, rtwo in zip(files_one, files_two):
        basename = os.path.basename(rone).split('_1')[0]
        ccmd = 'cat {0} {1} > {2}/{3}_merged.fq'.format(rone, rtwo, fq_out, basename)
        crun = subprocess.Popen(ccmd, shell=True)
        crun.wait()
        merged_files.append('{0}/{1}_merged.fq'.format(fq_out, basename))
    ref_list = '{0}/Reference_list.txt'.format(out_path)
    qry_list = '{0}/Query_list.txt'.format(out_path)
    out_file = '{0}/fastANI_output.tsv'.format(out_path)
    ref_file = open(ref_list, 'w')
    qry_file = open(qry_list, 'w')
    for paths in merged_files:
        ref_file.write('{0}\n'.format(os.path.abspath(paths)))
        qry_file.write('{0}\n'.format(os.path.abspath(paths)))
    ref_file.close()
    qry_file.close()
    facmd = ['fastANI', '-t', '2', '--ql', qry_list, '--rl', ref_list,
             '--matrix', '-o', out_file]
    farun = subprocess.Popen(facmd, shell=False)
    farun.wait()
    stop = time.time() - start
    print('Mode: fastANI; Input: {0}; Runtime: {1}'.format(input_path, stop))

def run_mash(input_path, out_path):
    '''Run MaSH to get pairwise distance'''
    start = time.time()
    files = glob.glob('{0}/*.fasta'.format(input_path))
    ref_list = '{0}/Reference_list.txt'.format(out_path)
    qry_list = '{0}/Query_list.txt'.format(out_path)
    out_file = '{0}/Mast_output.tsv'.format(out_path)
    skh_file = '{0}/Mash_sketch.msh'.format(out_path)
    skh_base = '{0}/Mash_sketch'.format(out_path)
    #Ref and query list
    ref_file = open(ref_list, 'w')
    qry_file = open(qry_list, 'w')
    for paths in files:
        ref_file.write('{0}\n'.format(os.path.abspath(paths)))
        qry_file.write('{0}\n'.format(os.path.abspath(paths)))
    ref_file.close()
    qry_file.close()
    #Sketch command
    skcmd = ['mash', 'sketch', '-o', skh_base, '-p', '10', '-l', ref_list]
    skrun = subprocess.Popen(skcmd, shell=False)
    skrun.wait()
    #Dist command
    mhcmd = 'mash dist -l -p 10 {1} {0} > {2}'.format(qry_list, skh_file, out_file)
    mhrun = subprocess.Popen(mhcmd, shell=True)
    mhrun.wait()
    stop = time.time() - start
    print('Mode: Mash; Input: {0}; Runtime: {1}'.format(input_path, stop))

def run_mashfq(input_path, out_path):
    '''Run MaSH to get pairwise distance'''
    start = time.time()
    files_one = sorted(glob.glob('{0}/Fastq/*_r1.fastq'.format(input_path)))
    files_two = sorted(glob.glob('{0}/Fastq/*_r2.fastq'.format(input_path)))
    fq_out = '{0}/MergedFq'.format(out_path)
    if not os.path.exists(fq_out):
        os.mkdir(fq_out)
    merged_files = list()
    for rone, rtwo in zip(files_one, files_two):
        basename = os.path.basename(rone).split('_1')[0]
        ccmd = 'cat {0} {1} > {2}/{3}_merged.fq'.format(rone, rtwo, fq_out, basename)
        crun = subprocess.Popen(ccmd, shell=True)
        crun.wait()
        print(ccmd)
        merged_files.append('{0}/{1}_merged.fq'.format(fq_out, basename))
    ref_list = '{0}/Reference_list.txt'.format(out_path)
    qry_list = '{0}/Query_list.txt'.format(out_path)
    out_file = '{0}/Mast_output.tsv'.format(out_path)
    skh_file = '{0}/Mash_sketch.msh'.format(out_path)
    skh_base = '{0}/Mash_sketch'.format(out_path)
    #Ref and query list
    ref_file = open(ref_list, 'w')
    qry_file = open(qry_list, 'w')
    for paths in merged_files:
        print(paths)
        ref_file.write('{0}\n'.format(os.path.abspath(paths)))
        qry_file.write('{0}\n'.format(os.path.abspath(paths)))
    ref_file.close()
    qry_file.close()
    #Sketch command
    skcmd = ['mash', 'sketch', '-o', skh_base, '-p', '10', '-l', ref_list]
    skrun = subprocess.Popen(skcmd, shell=False)
    skrun.wait()
    #Dist command
    mhcmd = 'mash dist -l -p 10 {1} {0} > {2}'.format(qry_list, skh_file, out_file)
    mhrun = subprocess.Popen(mhcmd, shell=True)
    mhrun.wait()
    stop = time.time() - start
    print('Mode: Mash; Input: {0}; Runtime: {1}'.format(input_path, stop))

def run_finch(input_path, out_path):
    '''Run Finch to calculate the pairwise distance'''
    start = time.time()
    #Create index
    dkc_path = '{0}/Dkc'.format(input_path) #create_index(input_path, out_path)
    #Run finch script
    out_file = '{0}/Finch_output.tsv'.format(out_path)
    fcmd = ['Rscript', 'finch.rscript', '-i', dkc_path, '-o', out_file]
    frun = subprocess.Popen(fcmd, shell=False)
    frun.wait()
    stop = time.time() - start
    print('Mode: Finch; Input: {0}; Runtime: {1}'.format(input_path, stop))
    
def run_kanalyze(varsets):
    '''Run k-analyze to create decimel KC files'''
    files = varsets[0]
    out_path = varsets[1]
    temp_path = varsets[2]
    ksize = varsets[3]
    sample = os.path.splitext(os.path.basename(files))
    out_file = '{0}/Dkc/{1}.dkc'.format(out_path, sample)
    temp_loc = '{0}/Dkc/{1}_tmp'.format(temp_path, sample)
    if not os.path.exists(temp_loc):
        os.mkdir(temp_loc)
    fq_ext = ['fq', 'fq.gz', 'fastq', 'fastq.gz']
    if files[0].split('.', 1) in fq_ext:
        kcmd = ['java', '-jar', '../lib/kanalyze/kanalyze.jar', 'count', '-t',
                '2', '-m', 'dec', '-k', str(ksize), '-c', 'kmercount:2',
                '-rcanonical', '--seqfilter', 'sanger:20', '--temploc', temp_loc,
                '-o', out_file] +  files
    else:
        kcmd = ['java', '-jar', '../lib/kanalyze/kanalyze.jar', 'count', '-t',
                '2', '-m', 'dec', '-k', str(ksize), '--temploc',
                temp_loc, '-o', out_file] + files
    krun = subprocess.Popen(kcmd, stderr=subprocess.DEVNULL,
                            stdout=subprocess.DEVNULL, shell=False)
    krun.wait()
    if krun.returncode != 0:
        return
    else:
        shutil.rmtree(temp_loc)
    return

def create_index(input_path, out_path):
    '''Parallely schedule k-analyze analysis'''
    files = glob.glob('{0}/*.fna'.format(input_path))
    dkc_path = '{0}/Dkc'.format(out_path)
    if not os.path.exists(dkc_path):
        os.mkdir(dkc_path)
    threads = 15
    pool = Pool(threads)
    pool.map(run_kanalyze, zip(files, repeat(out_path), repeat(out_path), repeat(31)))
    return dkc_path

if __name__ == '__main__':
    arguments = argparse.ArgumentParser(prog='runTools')
    arguments.add_argument('-i', '--input_path', type=str,
                           help='Path to input data',
                           default='~/scratch/thesis/aimtwo/ref_seq/FastaAni')
    arguments.add_argument('-o', '--output_path', type=str,
                           help='Path to output folder',
                           default='~/scratch/thesis/aimtwo/ref_seq')
    arguments.add_argument('-m', '--mode', type=str,
                           choices=['pyani', 'fastani', 'mash', 'finch', 'mashfq', 'fastqni'],
                           help='Mode of analysis')
    valued = arguments.parse_args()
    if valued.mode == 'pyani':
        run_pyani(valued.input_path, valued.output_path)
    elif valued.mode == 'fastani':
        run_fastani(valued.input_path, valued.output_path)
    elif valued.mode == 'mash':
        run_mash(valued.input_path, valued.output_path)
    elif valued.mode == 'finch':
        run_finch(valued.input_path, valued.output_path)
    elif valued.mode == 'mashfq':
        run_mashfq(valued.input_path, valued.output_path)
    elif valued.mode == 'fastqni':
        run_fastqni(os.path.abspath(valued.input_path), os.path.abspath(valued.output_path))
