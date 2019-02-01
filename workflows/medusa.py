#!/usr/bin/env python3

# medusa.py 
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

import subprocess
import argparse
import os
import sys
import yaml
import textwrap
import shlex
from colorama import Fore
from colorama import Style

def validate_arg_name(x):
    for c in x:
        if not c.isalnum() and c not in ['_']:
            raise ValueError('%s is not a valid parameter name' % x)

parser = argparse.ArgumentParser(description= """
DESCRIPTION

medusa is a pipeline for detecting cancer somatic variants in tumour-normal
sample pairs using next-generation-sequencing data. The typical input is raw
fastq files and the main output is VCF files for different types of somatic
variants together with quality control files.

{p} is the main entry point for the user to execute the pipeline.

- Use the manifest file to pass files and sample information to the pipeline

- Set which tasks should be executed via the TASKS option

Note that {p} only collects the user's option and submits a snakemake command.
To add or modify analyses, see the snakemake rules.
""".format(p= os.path.basename(__file__)), formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--manifest', '-m', required= True, help='''\
Manifest file linking sample libraries to input fastq or bam files and pairing
tumour and normal samples. See online docs for details''')

parser.add_argument('--directory', '-d', help='''\
Specify working directory. Relative paths will use this as their origin''')

parser.add_argument('--genome_build', '-gb', required= True, help='''\
Identifier of the genome build. E.g. hg38 ''')

parser.add_argument('--not_run', '-N', action= 'store_true', help= '''\
Do not run snakemake, only print the command that would be executed''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

# -----------------------------------------------------------------------------
# Create the parser dynamically on the basis of the yaml file(s)

task_arg= parser.add_argument_group('TASKS')

tasks_yaml= yaml.load(open(os.path.join(os.path.dirname(__file__), "tasks.yaml"), "r"))
for x in tasks_yaml:
    validate_arg_name(x)
    xhelp= tasks_yaml[x]['description']
    task_arg.add_argument('--%s' % x, 
            help= '\n'.join(textwrap.wrap(xhelp, width= 70)),
            action= 'store_true')
# -----------------------------------------------------------------------------
opts_arg= parser.add_argument_group('OPTIONS CONFIGURING TASKS')

opts_yaml= yaml.load(open(os.path.join(os.path.dirname(__file__), "tasks_opts.yaml"), "r"))
for x in opts_yaml:
    validate_arg_name(x)
    if 'type' not in opts_yaml[x]:
        opts_yaml[x]['type']= 'str'
    if opts_yaml[x]['type'] != 'file':
        xtype= eval(opts_yaml[x]['type'])
    else:
        xtype= str

    default= None
    if 'default' in opts_yaml[x]:
        default= opts_yaml[x]['default']
    
    xhelp= opts_yaml[x]['description']
    if default is not None:
        xhelp= xhelp.rstrip('.')
        xhelp= xhelp + '. Default %s' % default
    xhelp= '\n'.join(textwrap.wrap(xhelp, width= 70))

    opts_arg.add_argument('--%s' % x, help= xhelp, default= default,
            type= xtype)

# -----------------------------------------------------------------------------
# Options passed to snakemake verbatim. If adding more options, remeber to
# update the snakemake string passed to subprocess

smk_arg= parser.add_argument_group('OPTIONS TO SNAKEMAKE')

smk_arg.add_argument('--config', '-C', nargs= '*', metavar="KEY=VALUE", default= [], help='''\
Further configuration options passed to the pipeline''')

smk_arg.add_argument('--dryrun', '-n', action= 'store_true', help= '''\
Do not execute jobs, just print the execution plan''')

smk_arg.add_argument('--printshellcmds', '-p', action= 'store_true', help= '''\
Print out the shell commands that will be executed''')

smk_arg.add_argument('--snakemake', '-S', metavar= "SNAKEMAKE='OPTS'", help='''\
Further snakemake options not included above. This is a single string passed
as-is to the snakemake command.  e.g. -S="--unlock -q"''')

# -----------------------------------------------------------------------------

def validate_args(args):
    """Validate the arguments in the Namespace args. This validation gives a
    useful message to the user if some inconsistencies are found
    """
    if args.bwa:
        if args.ref is None:
            sys.stderr.write(f'{Fore.RED}A reference file is required to execute task bwa\n{Style.RESET_ALL}')
            sys.exit(1)

    if args.cnv_facets:
        if args.cnv_facets_snp is None:
            sys.stderr.write(f'{Fore.RED}Please provide a VCF of SNPs to execute task cnv_facets\n{Style.RESET_ALL}')
            sys.exit(1)

    if args.calculateContamination:
        if args.calculateContamination_vcf is None:
            sys.stderr.write(f'{Fore.RED}Please provide a VCF of SNPs to execute task calculateContamination\n{Style.RESET_ALL}')
            sys.exit(1)

def compile_snakemake_cmd(args):
    """Return a string ready to be executed as a subprocess.
    args: Namespace object as produced by argparse.
    """
    args_dict= vars(args)
    requested_tasks= []
    task_options= []
    for opt in args_dict:
        if opt in tasks_yaml and args_dict[opt]:
            requested_tasks.append(opt)
        if opt in opts_yaml and args_dict[opt] is not None:
            if opts_yaml[opt]['type'] == 'file' and os.path.exists(args_dict[opt]):
                args_dict[opt]= os.path.abspath(args_dict[opt])
            try:
                task_options.append(opt + '=' + shlex.quote(args_dict[opt]))
            except TypeError:
                task_options.append(opt + '=' + str(args_dict[opt]))
    smk_config= ['medusa_dir=%s' % os.path.abspath(os.path.dirname(__file__))]
    smk_config.append('genome_build=%s' % args.genome_build)
    smk_config.append('manifest=%s' % os.path.abspath(args.manifest))
    smk_config.append("tasks='%s'" % ' '.join(requested_tasks))
    smk_config.extend(task_options)
    smk_config.extend(args.config)
    smk_config= ' \\\n        '.join(smk_config)

    # Start compiling the snakemake command
    smk_cmd= r"""
snakemake -d %(xdir)s \
    --config \
        %(smk_config)s \
    -s %(main)s""" % {'xdir': xdir, 'smk_config': smk_config, 'main': os.path.join(os.path.dirname(__file__), 'rules/main.smk')}

    if args.dryrun:
        smk_cmd += ' \\\n   --dryrun'
    if args.printshellcmds:
        smk_cmd += ' \\\n   --printshellcmds'

    if args.snakemake is not None and args.snakemake.strip() != '':
        smk_cmd += ' \\\n    ' + args.snakemake
    
    return smk_cmd


if __name__ == '__main__':
    args= parser.parse_args()
    validate_args(args)
    if args.directory is None:
        xdir= os.getcwd()
    else:
        xdir= args.directory
    
    # Gather configuration options that will be passed to snakemake via --config
    args_dict= vars(args)
    requested_tasks= []
    task_options= []
    for opt in args_dict:
        if opt in tasks_yaml and args_dict[opt]:
            requested_tasks.append(opt)
        if opt in opts_yaml and args_dict[opt] is not None:
            if opts_yaml[opt]['type'] == 'file' and os.path.exists(args_dict[opt]):
                args_dict[opt]= os.path.abspath(args_dict[opt])
            try:
                task_options.append(opt + '=' + shlex.quote(args_dict[opt]))
            except TypeError:
                task_options.append(opt + '=' + str(args_dict[opt]))
    smk_config= ['medusa_dir=%s' % os.path.abspath(os.path.dirname(__file__))]
    smk_config.append('genome_build=%s' % args.genome_build)
    smk_config.append('manifest=%s' % os.path.abspath(args.manifest))
    smk_config.append("tasks='%s'" % ' '.join(requested_tasks))
    smk_config.extend(task_options)
    smk_config.extend(args.config)
    smk_config= ' \\\n        '.join(smk_config)

    # Start compiling the snakemake command
    smk_cmd= r"""
snakemake -d %(xdir)s \
    --config \
        %(smk_config)s \
    -s %(main)s""" % {'xdir': xdir, 'smk_config': smk_config, 'main': os.path.join(os.path.dirname(__file__), 'rules/main.smk')}

    if args.dryrun:
        smk_cmd += ' \\\n   --dryrun'
    if args.printshellcmds:
        smk_cmd += ' \\\n   --printshellcmds'

    if args.snakemake is not None and args.snakemake.strip() != '':
        smk_cmd += ' \\\n    ' + args.snakemake

    sys.stderr.write(f'{Fore.MAGENTA}%s\n{Style.RESET_ALL}' % smk_cmd)

    if not args.not_run:
        sp= subprocess.Popen(smk_cmd, stderr= subprocess.PIPE, shell= True)
        stdout, stderr= sp.communicate()
        stderr= stderr.decode().split('\n')
        stderr= '\n'.join([x for x in stderr])
        if stderr != '':
            sys.stderr.write(stderr)

        if sp.returncode != 0:
            sys.exit(sp.returncode)
