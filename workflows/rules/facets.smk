# facets.smk 
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

rule cnv_facets:
    input:
        tumour= '{tumour_smp}/{aligner}/{tumour_smp}.bam',
        tbai= '{tumour_smp}/{aligner}/{tumour_smp}.bam.bai',
        normal= '{normal_smp}/{aligner}/{normal_smp}.bam',
        nbai= '{normal_smp}/{aligner}/{normal_smp}.bam.bai',
        snp_vcf= config['cnv_facets_snp'] if 'cnv_facets_snp' in config else '',
    output:
        vcf= '{tumour_smp}/{aligner}/facets/{normal_smp}/{tumour_smp}.vcf.gz',
        tbi= '{tumour_smp}/{aligner}/facets/{normal_smp}/{tumour_smp}.vcf.gz.tbi',
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.tumour_smp),
        jobname= '{tumour_smp}',
        cpus_per_task= 8,
        mem= 30000,
        gbuild= config['genome_build']
    shell:
        r"""
        out=`echo {output.vcf} | sed 's/\.vcf\.gz$//'`
        cnv_facets.R \
            --snp-tumour {input.tumour} \
            --snp-normal {input.normal} \
            --snp-vcf {input.snp_vcf} \
            --snp-mapq 5 \
            --snp-baq 10 \
            --snp-nprocs {params.cpus_per_task} \
            --depth 25 8000 \
            --cval 25 400 \
            --nbhd-snp auto \
            --gbuild {params.gbuild} \
            --out ${{out}}
        """

