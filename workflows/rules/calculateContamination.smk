# calculateContamination.smk 
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
#
rule getPileupSummaries:
    input:
        bam= '{sample}/{aligner}/{sample}.bam',
        bai= '{sample}/{aligner}/{sample}.bam.bai',
        snps= config['calculateContamination_vcf'] if 'calculateContamination_vcf' in config else '',
    output:
        '{sample}/{aligner}/calculateContamination/{sample}.pst',
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.sample),
        jobname= '{sample}',
        cpus_per_task= 2,
        mem= 4000,
    shell:
        r"""
        gatk --java-options '-Xmx3500m' GetPileupSummaries \
            -I {input.bam} \
            -V {input.snps} \
            -O {output} \
            --intervals {input.snps}
        """

rule calculateContamination:
    input:
        '{sample}/{aligner}/calculateContamination/{sample}.pst',
    output:
        '{sample}/{aligner}/calculateContamination/{sample}.cntm'
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.sample),
        jobname= '{sample}',
        cpus_per_task= 2,
        mem= 4000,
    shell:
        r"""
        gatk --java-options '-Xmx3500m' CalculateContamination \
            -I {input} \
            -O {output}
        """
