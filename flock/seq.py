import os
import sys
import glob
import time
import logging
from flock.fasta import Fasta
from flock.fastq import Fastq
from collections import namedtuple

#Create logger
FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)


class Seq:

    def __init__(self, input, meta=None):
        self.input = input
        self.read = self.createSeqObj()
        self.meta = meta

    def isFasta(self):
        if ('.fna' in self.input or '.fa' in self.input or '.fasta' in self.input
            or '.faa' in self.input) and '.fastq' not in self.input:
            return True
        else:
            return False

    def isFastq(self):
        if ('.fastq' in self.input or '.fastq.gz' in self.input or
            'fq' in self.input or 'fq.gz' in self.input):
            return True
        else:
            return False

    def isRaw(self):
        if not os.path.isfile(self.input) and isinstance(self.input, str):
            return True

    def isMer(self):
        if isinstance(self.input, int):
            return True

    def createSeqObj(self):
        if self.isFasta():
            seq_obj = Fasta(self.input).read()
            ftype = 'fasta'
        elif self.isFastq():
            seq_obj = Fastq(self.input).read()
            ftype = 'fastq'
        elif self.isRaw():
            if isinstance(self.input, str):
                seq_obj = iter(self.input.split())
            elif isinstance(self.input, list):
                seq_obj = iter(self.input)
            ftype = 'raw'
        elif self.isMer():
            seq_obj = iter([self.input])
            ftype = 'dec'

        Seq = namedtuple('Seq', ['header', 'seq', 'qual', 'length', 'type',
            'source'])
        for sequence in seq_obj:
            if ftype == 'fasta':
                header = sequence.header
                seq = sequence.seq
                qual = None
                length = sequence.length
                source = self.input
            elif ftype == 'fastq':
                header = sequence.header
                seq = ''.join(sequence.seq)
                qual = sequence.quals
                length = len(seq)
                source = self.input
            elif ftype == 'raw':
                header = None
                seq = sequence
                qual = None
                length = len(seq)
                source = None
            else:
                header = None
                seq = sequence
                qual = None
                length = self.meta
                source = None
            seq = Seq(header, seq, qual, length, ftype, source)
            yield(seq)
