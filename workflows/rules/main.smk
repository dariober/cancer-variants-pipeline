# main.smk 
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
import pandas 
import sys
import yaml
from lib import utils
from lib import Manifest

mnf= Manifest.Manifest(config['manifest'])
mnf.symlink_bam()

config['chroms']= []
if 'ref' in config:
    config['chroms']= utils.get_major_chroms(config['ref'])

print('----')
print(mnf)

wildcard_constraints:
    sample= '|'.join([re.escape(x) for x in mnf.sample]),
    tumour_smp= '|'.join([re.escape(x) for x in mnf.contrasts['tumour_smp']]),
    normal_smp= '|'.join([re.escape(x) for x in mnf.contrasts['normal_smp']]),
    chrom= '|'.join([re.escape(x) for x in config['chroms']]),
    fastq_md5= '|'.join([re.escape(x) for x in mnf.fastqc['fastq_md5']]),
    pair_id= '|'.join([re.escape(x) for x in mnf.fastq['pair_id']]),

known_tasks= yaml.load(open(os.path.join(config['medusa_dir'], 'tasks.yaml')))
tasks= {} 
for task in config['tasks'].split(' '):
    if task in known_tasks:
        tasks[task]= known_tasks[task]

print('----')
print('Output per task:')
[print(x, tasks[x]['output']) for x in tasks]
print('----')

rule all_task_output:
    input:
        [eval(tasks[x]['output']) for x in tasks]

include: 'mutect.smk'
include: 'facets.smk'
include: 'proc_bam.smk'
include: 'accessories.smk'
include: 'calculateContamination.smk'
if not mnf.fastq.empty:
    include: 'proc_fastq.smk'

onsuccess:
    for x in mnf.sample:
        try:
            os.removedirs(utils.get_log_dir(x, create= False))
        except OSError:
            pass
