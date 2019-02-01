# proc_bam.smk
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

bam_merge= "samtools merge -@ 8 {output.bam} {input.bam}"
bam_rename= "mv {input.bam} {output.bam}"
rule mergeOrRenameBam:
    # There may be multiple pairs of fastq files assiged to the same sample,
    # for example if the library was sequenced across multiple lanes. In this
    # case, we merge bam files. If instead there is just one pair we go for a
    # simple renaming.
    priority: -10,
    input:
        bam= lambda wc: expand('{sample}/{aligner}/{pair_id}.tmp.bam', sample= wc.sample, aligner= wc.aligner, pair_id= mnf.pair_ids_for_sample(wc.sample)),
    output:
        bam= temp('{sample}/{aligner}/{sample}.tmp.bam'),
    params:
        log_dir= lambda wc: utils.get_log_dir(wc.sample),
        jobname= '{sample}.{aligner}',
        cpus_per_task= 4,
        mem= 2000,
    message:
        "If the number of input files is > 1, merge input with:\n%s\n\nOtherwise rename input with\n%s" % (bam_merge, bam_rename) ,
    run:
        shell("samtools quickcheck {input.bam}")
        if len(input.bam) > 1:
            print(bam_merge)
            shell(bam_merge)
        else:
            print(bam_rename)
            shell(bam_rename)

rule bamToTdf:
    input:
        bam= '{sample}/{aligner}/{sample}.bam',
        genome= re.sub('\.fa$|\.fasta$', '.chrom.sizes', config['ref'])
    output:
        '{sample}/{aligner}/{sample}.tdf',
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.sample),
        jobname= '{sample}',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        igvtools count --minMapQuality 5 {input.bam} {output} {input.genome}
        """

rule samstats:
    input:
        '{sample}/{aligner}/{sample}.bam',
    output:
        '{sample}/{aligner}/{sample}.stats',
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.sample),
        jobname= '{sample}',
        cpus_per_task= 2,
        mem= 2000,
    shell:
        """
        samtools stats -@ {params.cpus_per_task} {input} > {output}
        """

rule markDuplicates:
    priority: -10,
    input:
        bam= '{sample}/{aligner}/{sample}.tmp.bam',
    output:
        bam= '{sample}/{aligner}/{sample}.bam',
        bai= '{sample}/{aligner}/{sample}.bam.bai',
        dupmetrics= '{sample}/{aligner}/{sample}.dupmetrics.txt',
        md5= '{sample}/{aligner}/{sample}.bam.md5',
    log:
        trimlog= '{sample}/{aligner}/{sample}.trim.log'
    params:
        log_dir= lambda wildcards: utils.get_log_dir(wildcards.sample),
        jobname= '{sample}.{aligner}',
        cpus_per_task= 4,
        mem= 8000,
        tmp_dir= '{sample}/{aligner}',
    shell:
        r"""
        gatk --java-options '-XX:ParallelGCThreads=8 -Xmx6g' MarkDuplicates \
            -I={input.bam} \
            -O={output.bam} \
            -M={output.dupmetrics} \
            --TMP_DIR={params.tmp_dir} \
            --MAX_RECORDS_IN_RAM=5000000 \
            --ADD_PG_TAG_TO_READS=false \
            --CREATE_INDEX=true
        
        # Rename index from x.bai to x.bam.bai
        bai=`echo {output.bam} | sed s/\.bam$/.bai/`
        mv ${{bai}} {output.bai}

        md5sum {output.bam} > {output.md5}
        """

