# mutect_pon.smk
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

from collections import OrderedDict
import re
import os
import pandas 
import sys
import utils

config= utils.set_default_config(config)

manifest= utils.read_manifest(config['manifest'])

if 'sample_type' in manifest:
    sample= list(OrderedDict.fromkeys(manifest[(manifest['sample_type'] == 'pon')]['sample']))
else:
    # If the sample_type is absent, every sample is taken as PON
    sample= list(set(manifest['sample']))

wildcard_constraints:
    sample= '|'.join([re.escape(x) for x in sample]),
    chrom= '|'.join([re.escape(x) for x in config['chroms']])

tasks= {
    'bam':        expand('{sample}/bam/{sample}.bam.bai', sample= sample),
    'mutect_pon': config['mutect_pon'] + '.tbi', 
}
tasks= utils.requested_tasks(config, tasks)

mutect_pon_log_dir= os.path.dirname(config['mutect_pon'])
onstart:
    """Create the directories where cluster logs will go. This should not be
    necessary as snakemake should prepare them. However, see
    https://stackoverflow.com/questions/50042226/output-log-file-to-cluster-option
    """
    for smp in sample:
        os.makedirs(os.path.join(os.getcwd(), smp, 'cluster_log'), exist_ok= True)
    
    os.makedirs(mutect_pon_log_dir, exist_ok= True)
    
onsuccess:
    try:
        os.removedirs(mutect_pon_log_dir)
    except:
        pass

    for smp in set(manifest['sample']):
        try:
            os.removedirs(os.path.join(os.getcwd(), smp, 'cluster_log'))
        except OSError:
            pass
    

rule all_mutect_pon:
    input:
        tasks.values()

normal_variants= expand('mutect_pon/{sample}.vcf.gz', sample= sample)
rule createSomaticPanelOfNormals:
    input:
        normal_variants
    output:
        vcf= config['mutect_pon'],
        idx= config['mutect_pon'] + '.tbi',
    params:
        log_dir= mutect_pon_log_dir,
        vcfs= ' '.join(['-vcfs ' + re.sub('.tbi$', '', x) for x in normal_variants]),
        jobname= 'panelOfNormals',
        cpus_per_task= 2,
        mem= 8000,
        tmp_dir= 'mutect_pon/',
    shell:
        """
        gatk --java-options '-Xmx7500m' CreateSomaticPanelOfNormals \\
            {params.vcfs} \\
            --output {output.vcf} \\
            --TMP_DIR {params.tmp_dir}
        """

rule mergeMutectNormalsChroms:
    input:
        vcf= lambda wildcards: expand("mutect_pon/{sample}.{chrom}.vcf.gz", sample= wildcards.sample, chrom= config['chroms']),
        tbi= lambda wildcards: expand("mutect_pon/{sample}.{chrom}.vcf.gz.tbi", sample= wildcards.sample, chrom= config['chroms'])
    output:
        vcf= 'mutect_pon/{sample}.vcf.gz',
        tbi= 'mutect_pon/{sample}.vcf.gz.tbi',
    params:
        log_dir= CLUSTER_LOG_DIR,
        jobname= '{sample}',
        cpus_per_task= 1,
        mem= 2000,
    shell:
        """
        bcftools concat {input.vcf} | bgzip > {output.vcf}
        tabix -f {output.vcf}
        """

rule mutectNormalChrom:
    input:
        bam= '{sample}/bam/{sample}.bam',
        bai= '{sample}/bam/{sample}.bam.bai',
        ref_dict= os.path.splitext(config['ref'])[0] + '.dict',
    output:
        vcf= temp('mutect_pon/{sample}.{chrom}.vcf.gz'),
        tbi= temp('mutect_pon/{sample}.{chrom}.vcf.gz.tbi'),
    params:
        log_dir= CLUSTER_LOG_DIR,
        chrom= '{chrom}',
        ref= config['ref'],
        sample= '{sample}',
        jobname= '{sample}.{chrom}',
        cpus_per_task= 1,
        mem= 3000,
        tmp_dir= 'mutect_pon/',
    shell:
        """
        gatk --java-options '-Xmx2500m' Mutect2 \\
            -R {params.ref} \\
            -I {input.bam} \\
            -tumor {params.sample} \\
            -L {params.chrom} \\
            -O {output.vcf} \\
            --TMP_DIR {params.tmp_dir}
        """

include: 
    os.path.join(os.path.abspath(config['workflows']), 'zzz_align.smk')
