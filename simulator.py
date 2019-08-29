import os
import sys
import glob
import itertools
import subprocess
import random
from itertools import repeat
from itertools import product
from multiprocessing import Pool

def get_length(sample):
    contigs = open(sample).read().split('>')[1:]
    length = 0
    for contig in contigs:
        seq = contig.split('\n')
        length += len(''.join(seq))
    return length

def dwgsim(arguments):
    error_rates = arguments[0]
    genomes = arguments[1]
    error_dict = {0.005 : '05p', 0.01: '1p', 0.05: '5p', 0.001: '01p'} 
    out_path = '/nv/hp10/sravishankar9/scratch/thesis/aimtwo/Plasmodium{0}/Fastq'.format(error_dict[error_rates])
    names = os.path.basename(genomes).split('.fasta')[0]
    size = get_length(genomes)
    reads = int(size*30/500)
    seed = random.randint(1,10000000)
    dcmd = ['dwgsim', '-1', '250', '-2', '250', '-d', '1000', '-e', str(error_rates),  
            '-E', str(error_rates), '-z', str(seed), '-c', '0', '-C', '-1', '-N', str(reads), '-I', '2', 
            genomes, '{2}/{0}{1}'.format(names, seed, out_path)]
    drun = subprocess.Popen(dcmd, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    drun.wait()
    if drun.returncode != 0:
        print(' '.join(dcmd))
    else:
        return([names, genomes, reads, error_rates, seed])

def dgparallel(inpath, outpath):
    reference = [os.path.abspath(files) for files in glob.glob('{0}/*.fasta'.format(inpath))]
    error_rates = [0.001, 0.01, 0.05, 0.005]
    iterations = list(itertools.product(*[error_rates, reference]))
    pools = Pool(20)
    details = pools.map(dwgsim, iterations)
    out_file = open('Plasmodium_insilico_study_details.tsv', 'w')
    out_file.write('Number\tSample_name\tGenome_file\tRead_count\tError_rate\tFastq_basename\n')
    for counter, exp in enumerate(details):
        out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{1}{5}\n'.format(counter, exp[0], exp[1], exp[2], exp[3], exp[4]))

    out_file.close()
if __name__ == '__main__':
    inpath = sys.argv[1]
    outpath = sys.argv[2]
    dgparallel(inpath, outpath)

