#!/usr/bin/env python3

# bootstrap_test_data.py
#
# Copyright (C) 2019 University of Glasgow
#
# Author: Dario Beraldi <dario.beraldi@glasgow.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import argparse
import os
import sys
import subprocess
import tempfile
import pandas
import re

parser = argparse.ArgumentParser(description= """
DESCRIPTION
Create test fastq files. Only for testing.
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--manifest', '-m',
    required= True,
    help='''\
Get name of fastq files from this manifest file. File names ending in .gz will
be gzipped. NB: existing files will be overwritten!''')

parser.add_argument('--ref', '-r',
    required= True,               
    help='''\
Use this fasta reference file for creating reads.''')

parser.add_argument('--snv_spacing', '-vs',
    default= 100000,
    type= int,
    help='''\
Add a Single Nucleotide Variant every so many bases on each chrom in reference.
Default: %(default)s''')

parser.add_argument('--depth', '-d',
                   default= 30,
                   type= int,
                   help='''\
Target read depth. Default %(default)s''')

parser.add_argument('--read_length', '-rl',
    default= 75,
    type= int,
    help='''\
Read length. Default: %(default)s''')

parser.add_argument('--insert_size', '-is',
    default= 300,
    type= int,
    help='''\
Insert size: distance between start of R1 and end of R2. Default: %(default)s''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

args= parser.parse_args()

def text_wrap(x, width):
    """Wrap string x to given line width
    """
    wrapped= []
    xfrom= 0
    xto= width
    while xfrom < len(x):
        wrapped.append(x[xfrom:xto])
        xfrom += width
        xto += width
    return '\n'.join(wrapped)

# Read reference
sys.stderr.write('Reading reference %s...\n' % args.ref)
ref= dict()
chrom= None
with open(args.ref) as fa:
    for line in fa:
        line= line.strip()
        if line.startswith('>'):
            if chrom is not None:
                ref[chrom]= ''.join(seq)
            chrom= line.strip('>').split()[0]
            seq= []
        else:
            seq.append(line)
ref[chrom]= ''.join(seq)

# Add some variants to the reference
sys.stderr.write('Mutating reference...\n')
ref_var= dict()
fromto= {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'A'}
for chrom in ref:
    seq= list(ref[chrom])
    at= range(0, len(seq), args.snv_spacing)    
    sys.stderr.write('Variants on %s at position:\n' % chrom)
    sys.stderr.write('%s\n' % ' '.join([str(x+1) for x in at]))
    for pos in at:
        v= fromto[seq[pos].upper()]
        seq[pos]= v
    ref_var[chrom]= ''.join(seq)
    # Add some SVs
    seq= ref_var[chrom]
    seq_del= seq[0:int(len(seq)*0.2)] + seq[int(len(seq)*0.4):]
    ref_var[chrom]= seq_del

# Write the mutant reference
ref_var_fa= tempfile.NamedTemporaryFile(prefix= os.path.abspath(args.ref) + '.', suffix= '.var.fa', delete= False)
sys.stderr.write('Writing mutated reference to %s\n' % ref_var_fa.name)
with open(ref_var_fa.name, 'w') as fa:
    for chrom in ref_var:
        fa.write('>' + chrom + '\n')
        seq= text_wrap(ref_var[chrom], 50)
        fa.write(seq + '\n')

gsize= 0
for chrom in ref:
    gsize += len(ref[chrom])
N= int(gsize * args.depth / (args.read_length * 2))

mnf= pandas.read_csv(args.manifest, sep= '\t')
libs= set(mnf['library'])
for lib in libs:
    fq1= sorted(mnf[(mnf['library'] == lib) & (mnf['rmate'] == 'R1')]['fastq'])
    fq2= sorted(mnf[(mnf['library'] == lib) & (mnf['rmate'] == 'R2')]['fastq'])
    for i in range(0, len(fq1)):
        os.makedirs(os.path.dirname(fq1[i]), exist_ok=True)
        os.makedirs(os.path.dirname(fq2[i]), exist_ok=True)
        gzip1= False
        if fq1[i].endswith('.gz'):
            gzip1= True
            fq1[i]= re.sub('\.gz$', '', fq1[i])
        gzip2= False
        if fq2[i].endswith('.gz'):
            gzip2= True
            fq2[i]= re.sub('\.gz$', '', fq2[i])
        
        if list(mnf[mnf['library'] == lib]['sample_type'])[0] in ['Normal', 'PON']:
            cmd= ['wgsim', '-e 0.0001', '-N', str(N), '-1', str(args.read_length), '-2', str(args.read_length), '-d', str(args.insert_size), '-S', fq1[i], args.ref, fq1[i], fq2[i]]
            sys.stderr.write(' '.join(cmd) + '\n')
            subprocess.check_call(cmd, stdout= open(os.devnull, 'w'))
        elif list(mnf[mnf['library'] == lib]['sample_type'])[0] in ['Tumour']:
            # Normal reads
            tmp1= tempfile.NamedTemporaryFile(prefix= os.path.abspath(fq1[i]) + '.', suffix= '.fq', delete= False)
            tmp2= tempfile.NamedTemporaryFile(prefix= os.path.abspath(fq2[i]) + '.', suffix= '.fq', delete= False)
            cmd= ['wgsim', '-e 0.0001', '-N', str(int(N * 0.2)), '-1', str(args.read_length), '-2', str(args.read_length), '-d', str(args.insert_size), '-S', fq1[i]+'a', args.ref, tmp1.name, tmp2.name]
            sys.stderr.write(' '.join(cmd) + '\n')
            subprocess.check_call(cmd, stdout= open(os.devnull, 'w'))

            # Tumour reads
            tmp11= tempfile.NamedTemporaryFile(prefix= os.path.abspath(fq1[i]) + '.', suffix= '.fq', delete= False)
            tmp22= tempfile.NamedTemporaryFile(prefix= os.path.abspath(fq2[i]) + '.', suffix= '.fq', delete= False)
            cmd= ['wgsim', '-e 0.0001', '-N', str(int(N * 0.8)), '-1', str(args.read_length), '-2', str(args.read_length), '-d', str(args.insert_size), '-S', fq1[i]+'b', ref_var_fa.name, tmp11.name, tmp22.name]
            sys.stderr.write(' '.join(cmd) + '\n')
            subprocess.check_call(cmd, stdout= open(os.devnull, 'w'))

            subprocess.check_call('cat %s %s > %s' %(tmp1.name, tmp11.name, fq1[i]), shell= True)
            subprocess.check_call('cat %s %s > %s' %(tmp2.name, tmp22.name, fq2[i]), shell= True)
            os.remove(tmp1.name)
            os.remove(tmp11.name)
            os.remove(tmp2.name)
            os.remove(tmp22.name)
        else:   
            raise Exception()

        if gzip1:
            subprocess.check_call(['gzip', '-f', fq1[i]])
        if gzip2:
            subprocess.check_call(['gzip', '-f', fq2[i]])

os.remove(ref_var_fa.name)
