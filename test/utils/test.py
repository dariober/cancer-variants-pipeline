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
import pathlib
import shutil

# Import module. This is a test in itself:
# WGS_pipeline/test/utils
wdir= os.path.dirname(os.path.realpath(__file__))
base_path= os.path.dirname(
           os.path.dirname(wdir))
module_path= os.path.join(base_path, 'workflows/rules')
sys.path.insert(0, module_path)

from lib import utils

# ------------------------------------------------------------------------------

class TestUtils(unittest.TestCase):
    
    def testSymlink(self):
        shutil.rmtree('test_out_link', ignore_errors= True)
        os.makedirs('test_out_link/L1/bwa', exist_ok= True) # Output dir
        os.makedirs('test_out_link/source', exist_ok= True)
        pathlib.Path('test_out_link/source/L1.bam').touch() # Source file
        
        expected= 'test_out_link/L1/bwa/L1.bam'
        utils.symLink('test_out_link/source/L1.bam', expected)
        self.assertTrue(os.path.islink(expected))
        
        # symlink again: Ok
        utils.symLink('test_out_link/source/L1.bam', expected)

        # Source and dest are the same: Ok
        shutil.rmtree('test_out_link') 
        os.makedirs('test_out_link', exist_ok= True)
        pathlib.Path('test_out_link/L1.bam').touch()
        utils.symLink('test_out_link/L1.bam', './test_out_link/../test_out_link/L1.bam')
        self.assertTrue(os.path.islink('test_out_link/L1.bam') == False) # symlink not done
        self.assertTrue(os.path.exists('test_out_link/L1.bam'))

        # Symlink exists but it points to a different file: Fail
        shutil.rmtree('test_out_link') 
        os.makedirs('test_out_link/L1/bwa', exist_ok= True)
        os.makedirs('test_out_link/source', exist_ok= True)
        pathlib.Path('test_out_link/source/L1.bam').touch()
        pathlib.Path('test_out_link/source/L2.bam').touch()

        dest= 'test_out_link/L1/bwa/L1.bam'
        os.symlink('test_out_link/source/L2.bam', dest) 
        with self.assertRaises(Exception):
            utils.symLink('test_out_link/source/L1.bam', dest)

        # Destintation file exists: Fail
        shutil.rmtree('test_out_link') 
        dest= 'test_out_link/L1/bwa/L1.bam'
        os.makedirs('test_out_link/L1/bwa', exist_ok= True)
        pathlib.Path(dest).touch()
        
        os.makedirs('test_out_link/source', exist_ok= True)
        pathlib.Path('test_out_link/source/L1.bam').touch()
        with self.assertRaises(Exception):
            utils.symLink('test_out_link/source/L1.bam', dest)

        shutil.rmtree('test_out_link')

if __name__ == '__main__':
    unittest.main()
