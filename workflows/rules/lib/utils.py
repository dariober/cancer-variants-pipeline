# utils.py  
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

import re
import os
import subprocess
import sys
import datetime
import copy
import pandas
import tempfile
import distutils.spawn
import textwrap
from lib import errors
import pyfaidx

vep_script= r"""
base=`echo '{output.vcf}' | sed 's/.vcf.gz$//'`

vep --format vcf --fasta {input.ref} --hgvs --dont_skip --flag_pick --offline \
    --cache --dir_cache {params.dir_cache} --fork 4 --vcf --force_overwrite \
    -i {input.vcf} -o STDOUT --stats_file ${{base}}.summary.html \
    --warning_file ${{base}}.vep_warnings.txt \
| bgzip > {output.vcf}
tabix {output.vcf}
        """

def get_major_chroms(fasta):
    """Return a list of the "major" chromsosomes from the given fasta genome.
    Major chromsomes are those identified by digits, optionally prefixed by the
    string 'chr'. See code for details.
    """
    pyfaidx.Fasta(fasta) # To create the fasta index
    with open(fasta + '.fai') as f:
        chroms= []
        for line in f:
            chrom= line.strip().split('\t')[0]
            nchrom= re.sub('^chr', '', chrom)
            if nchrom.isnumeric() or nchrom in ['X', 'Y']:
                chroms.append(chrom)
    return chroms

def find_on_path(x):
    """Search user PATH to find file x. 
    """
    path= distutils.spawn.find_executable(x)
    if path is None:
        sys.stderr.write('\nFile "%s" not found on PATH\n\n' % x)
        sys.exit(1)
    else:
        return path

def get_bam_index(bam):
    """Return the index file associated to this bam file. Tipically this is
    `x.bam.bai` but it could be also `x.bai`
    """
    if os.path.exists(bam + '.bai'):
        return bam + '.bai'
    obam= re.sub('\.bam$', '', bam)
    if os.path.exists(obam + '.bai'):
        return obam + '.bai'
    raise FileNotFoundError('No index found for %s' % bam)

def symLink(source, dest):
    """Convenience function to create a symlink. If dest exists, we check
    whether the symlink would be redundant (ok) or if it points to a different
    source (error).
    """
    source= os.path.abspath(source)
    dest= os.path.abspath(dest)
    if source == dest:
        return
    try:
        os.symlink(source, dest)
    except FileExistsError:
        # We need to check if the existing file is the symlink we
        # tried to create. If it is not, we exit otherwise you
        # think you are analyising the bam in manifest and instead
        # you are analysing the one already there.
        if os.path.realpath(dest) != source:
            raise Exception('A file already exists in %s.' % dest)

def get_log_dir(root_dir, sub_dir= 'cluster_log', create= True):
    """Create and return the path to the log directory for this root_directory
    """
    ld= os.path.join(root_dir, sub_dir)
    if create:
        os.makedirs(ld, exist_ok= True)
    return ld

