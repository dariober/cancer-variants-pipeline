# proc_fastq.smk 
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

rule fastqc: 
    input: 
        lambda wc: mnf.fastqc[mnf.fastqc['fastq_md5'] == wc.fastq_md5]['fastq']
    output: 
        md5= '{fastq_md5}',
    params: 
        jobname= lambda wc: os.path.basename(wc.fastq_md5), 
        log_dir= lambda wc: utils.get_log_dir(list(mnf.fastqc[mnf.fastqc['fastq_md5'] == wc.fastq_md5]['sample'])[0]),
        fastqc_dir= lambda wc: os.path.dirname(wc.fastq_md5),
        cpus_per_task= 1, 
        mem= 2000, 
    shell: 
        """ 
        md5sum {input} > {output.md5}
        fastqc --noextract --outdir {params.fastqc_dir} {input} 
        """ 

rule bwa:
    priority: -10,
    input:
        config['ref'] + '.bwt',
        fq1= lambda wc: mnf.get_fq1(wc.pair_id),
        fq2= lambda wc: mnf.get_fq2(wc.pair_id),
    output:
        bam= temp('{sample}/bwa/{pair_id}.tmp.bam'),
    log:
        trim= '{sample}/bwa/{pair_id}.trim.log',
    params:
        log_dir= lambda wc: utils.get_log_dir(wc.sample),
        ref= config['ref'],
        RG= lambda wc: r'@RG\tID:{id}\tSM:{sm}\tLB:{lb}\tPL:ILLUMINA\tPU:NA'.format(
            id= wc.pair_id,
            lb= list(mnf.fastq[(mnf.fastq['pair_id'] == wc.pair_id)]['library'])[0],
            sm= wc.sample),
        jobname= '{sample}.{pair_id}',
        cpus_per_task= 16,
        mem= 48000,
        trimq= config['bbduk_trimq'],
        minlength= config['bbduk_minlength'],
    shell:
        r"""
        bbduk.sh -Xmx500m in={input.fq1} in2={input.fq2} out=stdout.fq qtrim=rl minlength={params.minlength} trimq={params.trimq} 2> {log.trim} \
        | bwa mem -R '{params.RG}' -t 12 -p {params.ref} - \
        | samtools sort -T {output.bam} -@ 4 > {output.bam}
        """

rule bowtie:
    priority: -10,
    input:
        config['ref'] + '.rev.2.bt2',
        fq1= lambda wc: mnf.get_fq1(wc.pair_id),
        fq2= lambda wc: mnf.get_fq2(wc.pair_id),
    output:
        bam= temp('{sample}/bowtie/{pair_id}.tmp.bam'),
    log:
        trim= '{sample}/bowtie/{pair_id}.trim.log',
    params:
        log_dir= lambda wc: utils.get_log_dir(wc.sample),
        ref= config['ref'],
        RG= lambda wc: '--rg-id {id} --rg SM:{sm} --rg LB:{lb} --rg PL:ILLUMINA --rg PU:NA'.format(
            id= wc.pair_id,
            lb= list(mnf.fastq[(mnf.fastq['pair_id'] == wc.pair_id)]['library'])[0],
            sm= wc.sample),
        jobname= '{sample}.{pair_id}',
        cpus_per_task= 16,
        mem= 48000,
        trimq= config['bbduk_trimq'],
        minlength= config['bbduk_minlength'],
    shell:
        """
        bbduk.sh -Xmx500m in={input.fq1} in2={input.fq2} out=stdout.fq qtrim=rl minlength={params.minlength} trimq={params.trimq} 2> {log.trim} \\
        | bowtie2 '{params.RG}' --maxins 1000 --threads 16 --local --interleaved - -x {params.ref} \\
        | samtools sort -T {output.bam} -@ 8 > {output.bam}
        """

