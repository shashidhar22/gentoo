import os
import sys
import glob
import subprocess
from itertools import repeat
from multiprocessing import Pool

def run_spades(values):
    rone_path, rtwo_path, out_path = values
    basename = os.path.basename(rone_path).split('_r1')[0]
    out_dir = '{0}/{1}'.format(out_path, basename)
    scmd = ['spades.py', '-t', '10', '-k', '21,33,55,77', '--careful', '-1', rone_path,
            '-2', rtwo_path, '-o', out_dir]
    srun = subprocess.Popen(scmd, shell=False)
    srun.wait()
    print(' '.join(scmd))
    out_file = '{0}/{1}_spades.fna'.format(out_path, basename)
    shutil.copy2('{0}/scaffolds.fasta', out_file)
    return

def spawn_spades(input_path, out_path):
    out_dir = '{0}/assemblies'.format(out_path)
    rones = sorted(glob.glob('{0}/*_r1.fastq'.format(input_path)))
    rtwos = sorted(glob.glob('{0}/*_r2.fastq'.format(input_path)))
    print(out_dir)
    print(rones)
    pools = Pool(1)
    pools.map(run_spades, zip(rones, rtwos, repeat(out_dir)))
    
if __name__ == '__main__':
    input_path = os.path.abspath(sys.argv[1])
    out_path = os.path.abspath(sys.argv[2])
    spawn_spades(input_path, out_path)    

