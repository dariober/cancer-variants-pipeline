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


"""Execute tests and check results with:
    ./tests.py
"""

import unittest
import os
import sys
import shutil

# Import module. This is a test in itself:
wdir= os.path.dirname(os.path.realpath(__file__))
base_path= os.path.dirname(
           os.path.dirname(wdir))
module_path= os.path.join(base_path, 'workflows/rules')
sys.path.insert(0, module_path)

from lib import Manifest
from lib import errors

os.chdir(wdir)

# ------------------------------------------------------------------------------

class ManifestTest(unittest.TestCase):
    
    def testCanGetLibraries(self):
        mnf= Manifest.Manifest(os.path.join(wdir, 'mnf.txt'))
        self.assertEqual(3, len(mnf.sample))

    def testLibraryNameSameAsSampleName(self):
        # Since manifest does not have 'library' we use the sample name as library id.
        mnf= Manifest.Manifest(os.path.join(wdir, 'mnf.txt'))
        self.assertEqual(list(mnf.df['sample']), list(mnf.df['library']))

    def testCanGetMatchedNormalContrasts(self):
        mnf= Manifest.Manifest(os.path.join(wdir, 'mnf.txt'))
        self.assertEqual(1, len(mnf.contrasts))
        self.assertEqual('N1', mnf.contrasts['normal_smp'][0])
        self.assertEqual('T1', mnf.contrasts['tumour_smp'][0])

    def testCanMatchedManyToMany(self):
        mnf= Manifest.Manifest(os.path.join(wdir, 'manyToMany.txt'))
        self.assertEqual(4, len(mnf.contrasts))

    def testCanGetFastqcTable(self):
        mnf= Manifest.Manifest(os.path.join(wdir, 'mnf.txt'))
        self.assertEqual(6, len(mnf.fastqc)) # mnf.txt has 4 rows  with a duplicate fastq pair
        self.assertEqual('N1/fastqc/N1_S4_R1_001.fastq.gz.md5', list(mnf.fastqc['fastq_md5'])[0])

        mnf= Manifest.Manifest(os.path.join(wdir, 'fastqc.txt'))
        self.assertEqual(2, len(mnf.fastqc)) 
        self.assertEqual('L1/fastqc/R1/reads.fastq.gz.md5', list(mnf.fastqc['fastq_md5'])[0])
        self.assertEqual('L1/fastqc/R2/reads.fastq.gz.md5', list(mnf.fastqc['fastq_md5'])[1])

    def testCanSymlinkBamFiles(self):
        shutil.rmtree('N1', ignore_errors= True)
        shutil.rmtree('P1', ignore_errors= True)
        shutil.rmtree('T1', ignore_errors= True)
        shutil.rmtree('T3', ignore_errors= True)
        mnf= Manifest.Manifest('manifest_bam.txt')
        mnf.symlink_bam()
        # Note Output tree is created in cwd
        self.assertTrue(os.path.islink('N1/bwa/N1.bam'))
        self.assertTrue(os.path.islink('N1/bwa/N1.bam.bai'))
        self.assertTrue(os.path.islink('T1/bwa/T1.bam'))
        self.assertTrue(os.path.islink('T1/bwa/T1.bam.bai'))
        shutil.rmtree('N1', ignore_errors= True)
        shutil.rmtree('P1', ignore_errors= True)
        shutil.rmtree('T1', ignore_errors= True)
        shutil.rmtree('T3', ignore_errors= True)

    def testFailMultipleBamsToSameSample(self):
        shutil.rmtree('N1', ignore_errors= True)
        shutil.rmtree('T1', ignore_errors= True)
        mnf= Manifest.Manifest('manifest_bam_dup.txt')
        with self.assertRaises(errors.InvalidManifestError):
            mnf.symlink_bam()

    def testFailOnManifestWithDuplicateFq(self):
        with self.assertRaises(errors.InvalidManifestError):
            mnf= Manifest.Manifest('invalid.txt')

    def testFailTooManySamplesToOneContrast(self):
        with self.assertRaises(errors.InvalidManifestError):
            mnf= Manifest.Manifest('contrast.txt')

    def test_tn_pair_can_be_null(self):
        mnf= Manifest.Manifest('contrast2.txt')
        self.assertEqual(0, len(mnf.contrasts)) 

    def testGetTumourSamples(self):
        mnf= Manifest.Manifest(os.path.join(wdir, 'contrast2.txt'))
        ts= mnf.get_tumour_samples()
        self.assertTrue(type(ts) == set)
        self.assertEqual(['B', 'C'], sorted(list(ts)))
        
        # W/o tumour samples
        mnf= Manifest.Manifest(os.path.join(wdir, 'fastqc.txt'))
        self.assertEqual(set(), mnf.get_tumour_samples())

if __name__ == '__main__':
    unittest.main()
