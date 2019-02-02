#!/usr/bin/env python3

# test.py
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

import unittest
import shutil
import os
import subprocess as sp
import gzip
import argparse
import yaml
import sys

wdir= os.path.dirname(os.path.realpath(__file__))
base_path= os.path.dirname(
           os.path.dirname(wdir))
module_path= os.path.join(base_path, 'workflows')
sys.path.insert(0, module_path)
import medusa

def read_bam(bam):
    """Read bam file and return it as a list of lines
    """
    p= sp.Popen("samtools view -h %s" % bam, shell= True, stdout= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0
    sam= stdout.decode().split('\n')
    return sam 

def read_vcf(vcf):
    with gzip.open(vcf) as gz:
        recs= gz.read().decode().split('\n')
    return recs

class Medusa(unittest.TestCase):

    def setUp(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    #def testCompileCmd(self):
    #    opts_yaml= yaml.load(open("../../workflows/tasks_opts.yaml", "r"))
    #    parser= argparse.ArgumentParser()
    #    parser.add_argument('--manifest')
    #    parser.add_argument('--genome_build')
    #    parser.add_argument('--config', default= [], nargs= '*')
    #    args= parser.parse_args(['--manifest', 'mnf.txt', '--genome_build', 'hg38'])    
    #    medusa.compile_snakemake_cmd(args)

    def testNotRunOption(self):
        p= sp.Popen(r"""../../workflows/medusa.py --not_run -m manifest.txt -gb hg38""",
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertEqual(0, len(stdout.decode()))
        self.assertTrue('snakemake ' in stderr.decode())

    def testDryRun(self):
        p= sp.Popen(r"""
            ../../workflows/medusa.py -m manifest.txt \
                    --dryrun \
                    --printshellcmds \
                    -gb hg38 \
                    --directory test_out \
                    --mutect \
                    --cnv_facets \
                    --samstats \
                    --tdf \
                    --fastqc \
                    --calculateContamination \
                    --ref ../../test_data/ref/hg38.chroms.fa \
                    --bbduk_trimq 15 \
                    --bbduk_minlength 20 \
                    --vep_dir_cache ../../test_data/vep \
                    --cnv_facets_snp ../../test_data/ref/gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz \
                    --calculateContamination_vcf ../../test_data/ref/gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz
        """, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(len(stdout.decode()) > 1000)

    def testBwa(self):
        p= sp.Popen(r"""
            ../../workflows/medusa.py -m manifest.txt \
                    -p \
                    -gb hg38 \
                    --directory test_out \
                    --bwa \
                    --bbduk_trimq 15 \
                    --bbduk_minlength 20 \
                    --ref ../../test_data/ref/hg38.chroms.fa
        """, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('minlength=20' in stderr.decode() and 'trimq=15' in stderr.decode())
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/TCRBOA6-T.bam.bai'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-N/bwa/TCRBOA6-N.bam.bai'))
        
        sam= read_bam("test_out/TCRBOA6-T/bwa/TCRBOA6-T.bam")

        # Two bam files have been merged
        self.assertTrue(len(sam) > 130000 and len(sam) < 150000)
        rg= [x for x in sam if x.startswith('@RG\t')]
        self.assertEqual(2, len(rg))
        self.assertTrue('@RG\tID:TCRBOA6-T.L1.' in '\n'.join(rg)) 
        self.assertTrue('@RG\tID:TCRBOA6-T.L2.' in '\n'.join(rg)) 
        
        self.assertTrue('@PG\tID:MarkDuplicates' in '\n'.join(sam))

        # Test normal sample
        sam= read_bam("test_out/TCRBOA6-N/bwa/TCRBOA6-N.bam")
        self.assertTrue(len(sam) > 150000 and len(sam) < 160000)
        self.assertTrue(len([x for x in sam if x.startswith('@RG\tID:TCRBOA6-N.L1')]) == 1) # Only one read group

    def testRunAllTasks(self):
        p= sp.Popen(r"""
            ../../workflows/medusa.py -m manifest.txt \
                    -p \
                    -gb hg38 \
                    --directory test_out \
                    --mutect \
                    --cnv_facets \
                    --samstats \
                    --tdf \
                    --fastqc \
                    --calculateContamination \
                    --ref ../../test_data/ref/hg38.chroms.fa \
                    --vep_dir_cache ../../test_data/vep \
                    --cnv_facets_snp ../../test_data/ref/gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz \
                    --calculateContamination_vcf ../../test_data/ref/gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz
        """, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        print(stderr)
        print(stdout)
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/mutect/TCRBOA6-N/TCRBOA6-T.vcf.gz.tbi'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/TCRBOA6-T.stats'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/TCRBOA6-T.tdf'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/fastqc/TCRBOA6-T.L1.R1_fastqc.html'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/calculateContamination/TCRBOA6-T.cntm'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/facets/TCRBOA6-N/TCRBOA6-T.vcf.gz.tbi'))
        
        # Test mutect vcf has been VEP'd
        vcf= read_vcf('test_out/TCRBOA6-T/bwa/mutect/TCRBOA6-N/TCRBOA6-T.vcf.gz')
        self.assertTrue(';CSQ=' in '\n'.join(vcf))

    def testRunWithoutTumourNormalPairs(self):
        """If samples are not in tumour-normal pairs, we can still run tasks
        that do not depend on pairing such as align, fastqc etc.
        """
        p= sp.Popen(r"""
            ../../workflows/medusa.py -m manifest_wo_pairs.txt \
                    -p \
                    -gb hg38 \
                    --directory test_out \
                    --fastqc \
                    --bwa \
                    --samstats \
                    --tdf \
                    --calculateContamination \
                    --ref ../../test_data/ref/hg38.chroms.fa \
                    --calculateContamination_vcf ../../test_data/ref/gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz
        """, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/TCRBOA6-T.bam.bai'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/TCRBOA6-T.stats'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/TCRBOA6-T.tdf'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/fastqc/TCRBOA6-T.L1.R1_fastqc.html'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/calculateContamination/TCRBOA6-T.cntm'))

    def testSampleSwap(self):
        """In addition to the regular tumour vs normal, swap normal and
        tumour samples  
        """
        p= sp.Popen(r"""
            ../../workflows/medusa.py -m manifest_swap.txt \
                    -p \
                    -gb hg38 \
                    --directory test_out \
                    --mutect \
                    --ref ../../test_data/ref/hg38.chroms.fa
        """, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/mutect/TCRBOA6-N/TCRBOA6-T.vcf.gz.tbi'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-N/bwa/mutect/TCRBOA6-T/TCRBOA6-N.vcf.gz.tbi'))
        
        # Check sample pairing is as expected. This relies on mutect adding the
        # sample entries 
        vcf= read_vcf('test_out/TCRBOA6-T/bwa/mutect/TCRBOA6-N/TCRBOA6-T.vcf.gz') 
        self.assertTrue('##normal_sample=TCRBOA6-N' in '\n'.join(vcf))
        self.assertTrue('##tumor_sample=TCRBOA6-T' in '\n'.join(vcf))

        vcf= read_vcf('test_out/TCRBOA6-N/bwa/mutect/TCRBOA6-T/TCRBOA6-N.vcf.gz') 
        self.assertTrue('##normal_sample=TCRBOA6-T' in '\n'.join(vcf)) # Note swap
        self.assertTrue('##tumor_sample=TCRBOA6-N' in '\n'.join(vcf))
    
    def testBamInput(self):
        p= sp.check_call("""
            mkdir -p test_out &&
            cp ../../test_data/ref/hg38.chroms.fa test_out/
            bwa index test_out/hg38.chroms.fa
            for x in TCRBOA6-N TCRBOA6-T
            do
            bwa mem test_out/hg38.chroms.fa ../../test_data/fastq/${x}.L1.R1.fq.gz \
                | samtools sort > test_out/${x}.bam && 
                samtools index test_out/${x}.bam
            done
            cp test_out/TCRBOA6-T.bam test_out/TCRBOA6-T2.bam
            cp test_out/TCRBOA6-T.bam.bai test_out/TCRBOA6-T2.bam.bai
            """, shell= True)

        # NB: To run mutect you would need RG groups in BAMs
        p= sp.Popen(r"""
            ../../workflows/medusa.py -m manifest_bam.txt \
                    -p \
                    -gb hg38 \
                    --directory test_out \
                    --samstats \
                    --cnv_facets \
                    --cnv_facets_snp ../../test_data/ref/gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz \
                    --ref test_out/hg38.chroms.fa
        """, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()

        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/TCRBOA6-T.stats'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T2/bwa/TCRBOA6-T2.stats'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-N/bwa/TCRBOA6-N.stats'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T/bwa/facets/TCRBOA6-N/TCRBOA6-T.vcf.gz.tbi'))
        self.assertTrue(os.path.exists('test_out/TCRBOA6-T2/bwa/facets/TCRBOA6-N/TCRBOA6-T2.vcf.gz.tbi'))

if __name__ == '__main__':
    unittest.main()
