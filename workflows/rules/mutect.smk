# mutect.smk 
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
 
rule mutect:
    input:
        tumour= '{tumour_smp}/{aligner}/{tumour_smp}.bam',
        tbai= '{tumour_smp}/{aligner}/{tumour_smp}.bam.bai',
        normal= '{normal_smp}/{aligner}/{normal_smp}.bam',
        nbai= '{normal_smp}/{aligner}/{normal_smp}.bam.bai',
        ref_dict= os.path.splitext(config['ref'])[0] + '.dict',
    output:
        vcf= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.{chrom}.tmp.vcf.gz'),
        tbi= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.{chrom}.tmp.vcf.gz.tbi'),
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.tumour_smp),
        pon= '--panel-of-normals ' + config['mutect_pon']  if 'mutect_pon' in config else '',
        ref= config['ref'],
        chrom= '{chrom}',
        tumour_smp= '{tumour_smp}',
        normal_smp= '{normal_smp}',
        jobname= '{tumour_smp}.{normal_smp}.{chrom}.{aligner}',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        gatk --java-options '-Xmx3500m' Mutect2 {params.pon} \\
            -R {params.ref} \\
            -I {input.tumour} \\
            -tumor {params.tumour_smp} \\
            -I {input.normal} \\
            -normal {params.normal_smp} \\
            -L {params.chrom} \\
            -O {output.vcf}
        """

rule mergeMutect:
    input:
        vcf= lambda wildcards: expand("{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.{chrom}.tmp.vcf.gz", 
                                   tumour_smp= wildcards.tumour_smp, 
                                   normal_smp= wildcards.normal_smp, 
                                   chrom= config['chroms'], 
                                   aligner= wildcards.aligner),
        tbi= lambda wildcards: expand("{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.{chrom}.tmp.vcf.gz.tbi", 
                                   tumour_smp= wildcards.tumour_smp, 
                                   normal_smp= wildcards.normal_smp, 
                                   chrom= config['chroms'], 
                                   aligner= wildcards.aligner),
    output:
        vcf= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp.vcf.gz'),
        tbi= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp.vcf.gz.tbi'),
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.tumour_smp),
        jobname= '{tumour_smp}',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        bcftools concat {input.vcf} | bgzip > {output.vcf}
        tabix -f {output.vcf}
        """

rule filterMutectCalls:
    input:
        vcf= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp.vcf.gz',
        tbi= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp.vcf.gz.tbi',
    output:
        vcf= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp2.vcf.gz'),
        tbi= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp2.vcf.gz.tbi'),
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.tumour_smp),
        jobname= '{tumour_smp}',
        cpus_per_task= 2,
        mem= 4000,
    shell:
        """
        gatk --java-options '-Xmx3500m' FilterMutectCalls --variant {input.vcf} --output {output.vcf}
        """

if 'vep_dir_cache' in config and config['vep_dir_cache'] != '':
    normalized_vcf_out= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.noVep.vcf.gz')
    normalized_tbi_out= temp('{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.noVep.vcf.gz.tbi')
else:
    normalized_vcf_out= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.vcf.gz'
    normalized_tbi_out= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.vcf.gz.tbi'

rule normalizeMutect:
    input:
        vcf= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp2.vcf.gz',
        tbi= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.tmp2.vcf.gz.tbi',
        ref= config['ref'],
    output:
        vcf= normalized_vcf_out,
        tbi= normalized_tbi_out,
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.tumour_smp),
        jobname= '{tumour_smp}',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        bcftools norm --output-type z --fasta-ref {input.ref} --multiallelics - --output {output.vcf} {input.vcf}
        tabix {output.vcf}
        """

if 'vep_dir_cache' in config:
    from lib import utils
    rule vepMutect:
        input:
            vcf= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.noVep.vcf.gz',
            tbi= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.noVep.vcf.gz.tbi',
            ref= config['ref'],
        output:
            vcf= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.vcf.gz',
            tbi= '{tumour_smp}/{aligner}/mutect/{normal_smp}/{tumour_smp}.vcf.gz.tbi',
        params:
            log_dir= lambda wildcards: utils.get_log_dir(wildcards.tumour_smp),
            jobname= '{tumour_smp}',
            cpus_per_task= 4,
            mem= 4000,
            dir_cache= config['vep_dir_cache'],
        shell:
            utils.vep_script
