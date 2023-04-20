#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to prepare and run recent HGT detection automatically.
# Created by galaxy on 2016/10/19 11:19
#
# DEPENDENCIES
# =============================================================================
# Python (>=2.7):
# o Biopython (http://www.biopython.org)
# o Pyani (https://github.com/widdowquinn/pyani)
# R (>=3.0.1):
# o fitdistplus, ggplot2
# Software:
# o EMBOSS Needle (http://emboss.sourceforge.net/download/)
# o MUMmer (http://mummer.sourceforge.net/)
# =============================================================================
#
# USAGE
# =============================================================================
# calculate_ani.py [options]
#
# Options:
#   -h, --help            show this help message and exit
#   -o OUTDIRNAME, --outdir=OUTDIRNAME
#                         Output directory
#   -i INDIRNAME, --indir=INDIRNAME
#                         Input directory name
#   -v, --verbose         Give verbose output
#   -f, --force           Force file overwriting
#   --noclobber           Don't nuke existing files
#   -p, --part            Which part will be run? [0|1|2|3|4]
#   -t, --threads         How many threads will be used? [default all]
#   -g GFORMAT            Graphics output format(s) [pdf|png|jpg|svg]
#   -l, --logfile         Logfile location
# ==============================================================================
#
# The MIT License
#
# Copyright (c) 2017-2018 Northwest A&U University
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
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ===============================================================================


import os
import re
import sys
import time
import shutil
import subprocess
import logging.handlers
import traceback
import tarfile
import itertools
import numpy as np
import pandas as pd
from os import path
from argparse import ArgumentParser
from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from subprocess import *
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool



def parse_cmdline():
    """
    Parse command-line arguments for script.
    :return: Input command-line arguments
    """
    parser = ArgumentParser(prog="recentHGT.py")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None,
                        help="Input directory name")
    parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None,
                        help="Output directory")
    parser.add_argument("-t", "--threads", type=int, dest="threads", default=cpu_count(),
                        help="How many threads will be used? [default all]")
    parser.add_argument("-p", "--part", type=int, dest="part", default=0,
                        help="Which part will be run? [0|1|2|3|4|5]")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile", action="store", default=None,
                        help="Logfile location")
    parser.add_argument("-f", "--force", dest="force", action="store_true", default=False,
                        help="Force file overwriting")
    parser.add_argument("--noclobber", dest="noclobber", action="store_true", default=False,
                        help="Don't nuke existing files")
    parser.add_argument("-g", dest="gformat", action="store", default="pdf",
                        help="Graphics output format(s) [pdf|png|jpg|svg]")
    parser.add_argument("-d", dest="drawing", action="store", default="t",
                        help="Do we need to draw some distribution pictures? "
                             "[default 'f' or 'F', you can set 't' or 'T']")
    return parser.parse_args()


def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))


def make_outdir():
    """Make the output directory, if required.

    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:

    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error("Output directory %s would overwrite existing " +
                         "files (exiting)", args.outdirname)
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything below it",
                        args.outdirname)
            if args.noclobber:
                logger.warning("NOCLOBBER: not actually deleting directory")
            else:
                shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s", args.outdirname)
    try:
        os.makedirs(args.outdirname)  # We make the directory recursively
        # Depending on the choice of method, a subdirectory will be made for
        # alignment output files
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)


def load_strains_info(strain_file, form=1):
    """
    This function is used to load strains information.
    :param form: using different form to return useful information to the functions latter.
    :param strain_file: the path of the file contains the ids and names of strains.
    :return: a python dict
    """
    logger.info('Loading strain information, form {0}'.format(str(form)))
    strain_dict = defaultdict()
    try:
        with open(strain_file, 'r') as f:
            for a_line in f.readlines()[1:]:
                a_list = a_line.strip().split('\t')
                #strain_name = a_list[1].split(' ')[-1]
                strain_name = a_list[1]
                strain_id = a_list[2]
                if len(a_list) > 3:
                    chr_id = a_list[3]
                    psym_id = a_list[4]
                if form == 1:
                    strain_dict[strain_id] = [strain_name, chr_id, psym_id]
                elif form == 2:
                    strain_dict[strain_name] = [strain_id, chr_id, psym_id]
                elif form == 3:
                    strain_dict[strain_id] = [strain_name]
    except IOError:
        logger.error("There is no file contains strain information or the file is locked, please check.")
        logger.error(last_exception())
        sys.exit(1)
    return strain_dict


def compress_delete_outdir(outdir):
    """Compress the contents of the passed directory to .tar.gz and delete."""
    # Compress output in .tar.gz file and remove raw output
    tar_fname = outdir + '.tar.gz'
    logger.info("\tCompressing output from %s to %s", outdir, tar_fname)
    with tarfile.open(tar_fname, "w:gz") as fh:
        fh.add(outdir)
    logger.info("\tRemoving output directory %s", outdir)
    shutil.rmtree(outdir)


def load_genome(genbank, strain_chr_id, strain_psym_id):
    pattern_1 = re.compile(r'LOCUS\s+(.*?)\s+')
    pattern_2 = re.compile(r'/db_xref="SEED:fig\|(.*)"')
    genome_dict = defaultdict()
    locus_value = 0
    with open(genbank, 'r') as f1:
        for each_line in f1.readlines():
            if 'LOCUS' in each_line:
                m = re.search(pattern_1, each_line.strip())
                locus = m.group(1)
                if locus == strain_chr_id:
                    locus_value = '*'
                elif locus == strain_psym_id:
                    locus_value = '#'
                else:
                    locus_value = '+'
            elif '/db_xref="SEED:' in each_line:
                n = re.search(pattern_2, each_line.strip())
                gene_id = n.group(1)
                genome_dict[gene_id] = locus_value  # key: *.peg.*, value: +,#,*,
    return genome_dict


def og_location(og_id, all_gene_dict, pair_gene_dir):
    if not os.path.exists(pair_gene_dir):
        logger.error("There is no directory contains gene file, please check.")
        logger.error(last_exception())
        sys.exit(1)
    tmp_gene_fasta = os.path.join(pair_gene_dir, og_id + '.fa')
    re_pattern = re.compile(r'(\d+\.\d+)_(\d+)')
    og_loc_value = ''  # 0: Both on Chromosome, 1: One on Plasmid, the other on Chromosome, 2: Both on Plasmid.
    for record in SeqIO.parse(tmp_gene_fasta, 'fasta'):
        m = re.search(re_pattern, record.description)
        gene_id = m.group(1)
        location_value = all_gene_dict[gene_id]
        og_loc_value += location_value
    return og_loc_value


def load_ani():
    ani_out_dir = os.path.join(args.outdirname, 'ANIm')
    if not os.path.exists(ani_out_dir):
        logger.info('IOError: try to open ANIm directory but failed, please check if it exists or renamed.')
        logger.error(last_exception())
        sys.exit(1)
    strain_ani_matrix_file = os.path.join(ani_out_dir, 'ANIm_percentage_identity.tab')
    matrix_strain_dict = defaultdict()
    strain_label_file = os.path.join(args.indirname, 'strain_info.txt')
    strain_dict = load_strains_info(strain_label_file, form=3)
    tmp_matrix_lines = ''
    with open(strain_ani_matrix_file, 'r') as f1:
        all_lines = f1.readlines()
        strain_list = all_lines[0].strip().split('\t')
        for strain_id in strain_list:
            strain_name = strain_dict[strain_id][0]
            strain_index = strain_list.index(strain_id)
            matrix_strain_dict[strain_name] = int(strain_index)  # key: strain name, e.g. CFN42, value: index, e.g. 0,1
        for each_line in all_lines[1:]:
            tmp_list = each_line.strip().split('\t')
            tmp_line = ''
            for i in tmp_list[1:]:
                ani = float('%.3f' % float(i)) * 100
                tmp_line += '{0}\t'.format(str(ani))
            tmp_matrix_lines += tmp_line.strip('\t') + '\n'
    tmp_matrix_file = os.path.join(args.outdirname, 'tmp.txt')
    with open(tmp_matrix_file, 'w') as f2:
        f2.write(tmp_matrix_lines)
    ani_matrix = np.loadtxt(tmp_matrix_file, delimiter='\t')
    os.remove(tmp_matrix_file)
    return ani_matrix, matrix_strain_dict, strain_dict


def draw_distribution(item, result_file, r):
    devnull = open(os.devnull, 'w')
    try:
        subprocess.call(['Rscript', r, item[2], result_file, item[0]],
                        stdout=devnull, stderr=devnull)
    except OSError:
        logger.info('Try to run {0} but failed, please check.'.format(r))
        logger.error(last_exception())
        sys.exit(1)


def each_needle_run(pair_gene_dir, tmp_gene_converted_dir, pair_gene_alignment_dir, og_id, strain_dict):
    """
    This function is used to call Needle program to do pairwise sequence alignment
    :param pair_gene_dir: each homologous gene directory
    :param tmp_gene_converted_dir: used to put some temporary files and will be deleted in the end
    :param pair_gene_alignment_dir: each orthologous gene pair-wised alignment directory
    :param og_id: each orthologous gene id
    :param strain_dict: inherit from load_strains_label function with strain information
    :return: the alignment result of each gene
    """
    if not os.path.exists(pair_gene_dir):
        logger.error("There is no directory contains gene file, please check.")
        logger.error(last_exception())
        sys.exit(1)
    tmp_gene_fasta = os.path.join(pair_gene_dir, og_id + '.fa')
    converted_records = []
    #记得调整
    #re_pattern = re.compile(r'(GCA_\d+\.\d+)_(.*)_(.*)_(.*)')
    #in_pattern = re.compile(r'Identity.*\((\d+\.\d+)%\)')
    annotation = ''
    og_list = []
    for record in SeqIO.parse(tmp_gene_fasta, 'fasta'):
        recordid = record.id
        strain_id = recordid.split('-')[0]
        #locustag = recordid.split('-')[3]
        locustag = recordid.split('-')[2]
        #logger.info(strain_id)
        #gene_id = '{0}.peg.{1}'.format(strain_id, locustag)
        gene_id = locustag
        og_list.append(gene_id)
        #annotation = m.group(3)
        record.id = strain_dict[strain_id][0]
        final_record = SeqRecord(record.seq, record.id, description='')
        converted_records.append(final_record)
    the_strain_fasta = os.path.join(tmp_gene_converted_dir, 'a.fasta')
    other_strain_fasta = os.path.join(tmp_gene_converted_dir, 'b.fasta')
    SeqIO.write(converted_records[0], the_strain_fasta, 'fasta')
    SeqIO.write(converted_records[1], other_strain_fasta, 'fasta')
    result_file = os.path.join(pair_gene_alignment_dir, "{0}.txt".format(og_id))
    needle_cline = NeedleCommandline()
    needle_cline.asequence = the_strain_fasta
    needle_cline.bsequence = other_strain_fasta
    needle_cline.gapopen = 10
    needle_cline.gapextend = 0.5
    needle_cline.outfile = result_file
    devnull = open(os.devnull, 'w')
    try:
        subprocess.call(str(needle_cline), shell=True, stdout=devnull, stderr=devnull)
    except OSError:
        logger.info('Try to call Needle program failed, please check if Needle has been installed successfully.')
        logger.error(last_exception())
        sys.exit(1)
    os.remove(the_strain_fasta)
    os.remove(other_strain_fasta)
    gene_alignment_result = ''
    with open(result_file, 'r') as f:
        for a_line in f.readlines():
            if 'Identity' in a_line:
                #m = re.search(in_pattern, a_line.strip())
                #identity = m.group(1)
                identity = a_line.strip().split('(')[1].split('%')[0]
                gene_alignment_result = '{0}\t[{1}|{2}]\t{3}\t{4}\n'.format(og_id, og_list[0],
                                                                            og_list[1], identity,
                                                                            annotation)
    return gene_alignment_result


def setPairOrthofinder(fasta_folder, strain_pair_dir, CDS_Dict):
    retval = os.getcwd()
    outdirname = path.join(retval, args.outdirname)
    # 10.27
    (ani_matrix, matrix_strain_dict, strain_dict) = load_ani()
    if not os.path.exists(strain_pair_dir):
        os.makedirs(strain_pair_dir)
    filenames = os.listdir(fasta_folder)
    fastapath = [os.path.join(fasta_folder, filename) for filename in filenames]
    strain_pair_fa = os.path.join(args.indirname, 'strain_pair')
    if not os.path.exists(strain_pair_fa):
        os.makedirs(strain_pair_fa)
    with open(args.outdirname + '/pair_ani_result.txt', 'a')as ani_file:
        for pair in combinations(fastapath, 2):
            strain_name_pair = [str(os.path.splitext(strain)[0]).split('/')[-1] for strain in pair]
            strain_pair_genes_dir = os.path.join(strain_pair_dir, str(strain_name_pair[0] + '-' + strain_name_pair[1]))
            pair1name = strain_dict[strain_name_pair[0]][0]
            pair2name = strain_dict[strain_name_pair[1]][0]
            pair_ani = ani_matrix[matrix_strain_dict[pair1name]][matrix_strain_dict[pair2name]]
            ani_file.writelines(pair1name+'\t'+pair2name+'\t'+str(pair_ani)+'\n')
            if pair_ani <= 95.0:
                each_strain_pair_path = os.path.join(strain_pair_fa, str(strain_name_pair[0] + '-' + strain_name_pair[1]))
                if not os.path.exists(each_strain_pair_path):
                    subprocess.call(['mkdir', each_strain_pair_path])
                for strain in pair:
                    shutil.copy(strain, each_strain_pair_path)
                devnull = open(os.devnull, 'w')
                if not os.path.exists(os.path.join(each_strain_pair_path, "OrthoFinder")):
                    try:
                        os.system("orthofinder -f " + each_strain_pair_path + " -t " + str(args.threads))
                        subprocess.call(['orthofinder', '-f', each_strain_pair_path, '-t', str(args.threads)], shell=True, stdout=devnull, stderr=devnull)
                    except OSError:
                        logger.info('Try to run' + str(strain_name_pair[0]) + '-' + strain_name_pair[1] + 'failed, please check.')
                        logger.error(last_exception())
                        sys.exit(1)
                    #continue
                ortho_path = os.path.join(each_strain_pair_path, "OrthoFinder")
                ortho_result = os.listdir(ortho_path)
                ortho_result_pathname = ortho_result[-1]
                ortho_result_path = os.path.join(ortho_path, ortho_result_pathname)
                singlecopy_OG_seq_path = path.join(ortho_result_path, "Single_Copy_Orthologue_Sequences")
                prokka_res_dir_name = path.join(outdirname, "prodigal_res")
                cds_folder = path.join(args.indirname, 'cds')
                p = Pool(args.threads)
                if not os.path.exists(strain_pair_genes_dir):
                    os.makedirs(strain_pair_genes_dir)
                    p.apply_async(GetCDS, args=(singlecopy_OG_seq_path, strain_pair_genes_dir, CDS_Dict))
                p.close()
                p.join()
            #GetCDS(singlecopy_OG_seq_path, prokka_res_dir_name, strain_pair_genes_dir)
            #for root, dirs, files in os.walk(singlecopy_OG_seq_path):
            #    for each_files in files:
            #        og_file_path = os.path.join(root, each_files)
            #        #logger.info(og_file_path)
            #        shutil.copy(og_file_path, strain_pair_genes_dir)
    message = 'Orthofinder-running have been finished.'
    return message


def runProdigal(fasta_folder):
    retval = os.getcwd()
    fastaFileList = os.listdir(fasta_folder)
    outdirname = path.join(retval, args.outdirname)
    prokka_res_dir_name = path.join(outdirname, "prodigal_res")

    os.chdir(path.join(fasta_folder, "../"))
    if path.exists(prokka_res_dir_name):
        logger.warning("prodigal res dir exists")
    else:
        logger.info(f"create prodigal res dir: {prokka_res_dir_name}")
        os.makedirs(prokka_res_dir_name)
    # 检查有没有prodigal命令
    try:
        # stderr=PIPE, stdout=PIPE 不要输出内容 把内容导导管道里面去
        ret = run("prodigal", shell=False, stderr=PIPE, stdout=PIPE)
        # print(ret)
        if ret.returncode < 0:
            print("Child was terminated by signal", ret.returncode, file=sys.stderr)
        else:
            print("Child returned", ret.returncode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)
        return
    print(f"{'#' * 10}Running prodigal{'#' * 10}\n{'-' * 60}")

    p = Pool(args.threads)
    for fna in fastaFileList:
        p.apply_async(prodigalCMD, args=(fasta_folder, fna, prokka_res_dir_name))
    p.close()
    p.join()
    os.chdir(retval)

def prodigalCMD(fasta_folder, fna, prokka_res_dir_name):
    fna_path = path.join(fasta_folder, fna)

    fna_name = ".".join(fna.split(".")[:-1])

    logger.info(f"Running prodigal for: {fna}")

    run(["prodigal", "-i", fna_path, "-a", f"{prokka_res_dir_name}/{fna_name}.faa", "-d",
         f"{prokka_res_dir_name}/{fna_name}.cds", "-o", f"{prokka_res_dir_name}/{fna_name}.gff"], stderr=PIPE,
        stdout=PIPE)
    # 重新命名
    seqtk_res_faa = run(["seqtk", "rename", f"{prokka_res_dir_name}/{fna_name}.faa", f"{fna_name}_"], stdout=PIPE)
    # 把seqtk的输出覆盖掉原始文件
    with open(f"{prokka_res_dir_name}/{fna_name}.faa", "wb") as f:
        f.write(seqtk_res_faa.stdout)
    # cds
    seqtk_res_cds = run(["seqtk", "rename", f"{prokka_res_dir_name}/{fna_name}.cds", f"{fna_name}_"], stdout=PIPE)
    with open(f"{prokka_res_dir_name}/{fna_name}.cds", "wb") as f:
        f.write(seqtk_res_cds.stdout)


def GetCDS(singlecopy_OG_seq_path, strain_pair_genes_dir, CDS_Dict):
    og_list = os.listdir(singlecopy_OG_seq_path)
    for filename in og_list:
        anum = str(filename)
        OG_number = anum[0:9]
        with open(strain_pair_genes_dir + '/' + OG_number + '.fa',
                  "a") as cds_result:
            pr_records = list(SeqIO.parse(
                singlecopy_OG_seq_path + '/' + filename,
                "fasta"))
            for pr_seq_record in pr_records:
                pr_id = str(pr_seq_record.id)
                #strain = pr_id.split('-')[0]
                #protein_id = pr_id.split('-')[1]
                #locus_tag = pr_id.split('-')[2]
                #s = pr_id.rfind("_")
                #pr_strain = pr_id[:s]
                cds_record = CDS_Dict[pr_id]
                cds_result.write('>' + str(pr_id) + '\n' + str(cds_record) + '\n')
    logger.info(singlecopy_OG_seq_path+' cds sequence done')


def getCDSDict():
    CDS_Dict = {}
    cds_folder = path.join(args.indirname, 'cds')
    filenames = os.listdir(cds_folder)
    for each_file in filenames:
        each_strain = os.path.splitext(each_file)[0]
        each_strain_dict = {}
        for record in SeqIO.parse(cds_folder +'/'+each_file,'fasta'):
            #根据cds文件进行调整
            cds_id = str(record.id)
            cds_seq = str(record.seq)
            #each_strain_dict[cds_id] = [cds_seq]
            #CDS_Dict[each_strain] = each_strain_dict
            CDS_Dict[cds_id] = [cds_seq]
    return CDS_Dict


def convertProteinandCDS0(genbank_folder, faa_folder):
    filenames = os.listdir(genbank_folder)
    cds_folder = os.path.join(args.indirname, 'cds')
    for each_file in filenames:
        logger.info(each_file)
        each_strain = os.path.splitext(each_file)[0]
        with open(faa_folder + '/' + each_strain + '.faa', 'a') as faa_out_file:
            #re_pattern = re.compile(r'fig\|(\d+\.\d+)\.peg\.(\d+)(.*)')
            for seq_record in SeqIO.parse(genbank_folder + '/' + each_file, "genbank"):
                cds_num = 1
                cds_fasta = ""
                complete_seq = str(seq_record.seq)
                for seq_feature in seq_record.features:
                    if seq_feature.type == 'CDS':
                        seq = seq_feature.qualifiers.get('translation')
                        if seq != None:
                            #logger.info(type(seq))
                            pseq = str(seq[0])
                            pr_seq = ''.join(pseq)
                            cds_seq = ""
                            f = seq_feature.qualifiers
                            print(f)
                            print(f['locus_tag'])

                            #seq_ana = '>' +each_strain + '-' + seq_record.id + '-' + seq_feature.qualifiers['protein_id'][0] + '-' + str(cds_num) + '\n'
                            # (10.5)
                            #seq_ana = '>' + each_strain + '-' + seq_record.id + '-' + '-' + str(cds_num) + '\n'
                            #seq_ana = '>' + each_strain + '-' + seq_record.id + '-' + seq_feature.qualifiers.get('protein_id') + '-' + str(cds_num) + '\n'
                            seq_ana = '>' + each_strain + '-' + seq_record.id + '-' + seq_feature.qualifiers['locus_tag'][0]+ '-' + str(cds_num) + '\n'
                            faa_out_file.write(seq_ana + pr_seq + '\n')
                            cds_num += 1
                            for ele1 in seq_feature.location.parts:
                                cds_seq += complete_seq[ele1.start:ele1.end]
                            cds_fasta += format_fasta(seq_ana, cds_seq, 70)
                cds_file_obj = open(cds_folder + '/' + each_strain + '.cds', "a")
                cds_file_obj.write(cds_fasta)
    message = 'protein file and cds file have been created.'
    return message

def convertProteinandCDS(genbank_folder, faa_folder):
    filenames = os.listdir(genbank_folder)
    cds_folder = os.path.join(args.indirname, 'cds')
    for each_file in filenames:
        logger.info(each_file)
        each_strain = os.path.splitext(each_file)[0]
        #with open(faa_folder + '/' + each_strain + '.faa', 'a') as faa_out_file:
        faa_file_obj = open(faa_folder + '/' + each_strain + '.faa', "a")
        #re_pattern = re.compile(r'fig\|(\d+\.\d+)\.peg\.(\d+)(.*)')
        for seq_record in SeqIO.parse(genbank_folder + '/' + each_file, "genbank"):
            cds_num = 1
            cds_fasta = ""
            com_seq = seq_record.seq
            complete_seq = str(seq_record.seq)
            for seq_feature in seq_record.features:
                if seq_feature.type == 'CDS':
                    seq = seq_feature.qualifiers.get('translation')
                    if seq != None:
                        #logger.info(type(seq))
                        pseq = str(seq[0])
                        pr_seq = ''.join(pseq)
                        cds_seq = ""
                        #seq_ana = '>' +each_strain + '-' + seq_record.id + '-' + seq_feature.qualifiers['protein_id'][0] + '-' + str(cds_num) + '\n'
                        seq_ana = '>' + each_strain + '-' + seq_record.id + '-' + seq_feature.qualifiers['locus_tag'][0]+ '-' + str(cds_num) + '\n'
                        #seq_ana = '>' + each_strain + '-' + seq_record.id + '-' + seq_feature.qualifiers.get('protein_id') + '-' + str(cds_num) + '\n'
                        faa_file_obj.write(seq_ana + pr_seq + '\n')
                        cds_num += 1
                        for ele1 in seq_feature.location.parts:
                            if seq_feature.strand == -1:
                                the_seq = com_seq[ele1.start:ele1.end]
                                complement_seq = the_seq.reverse_complement()
                                cds_seq += str(complement_seq)
                            if seq_feature.strand == 1:
                                cds_seq += complete_seq[ele1.start:ele1.end]
                        cds_fasta += format_fasta(seq_ana, cds_seq, 70)
            cds_file_obj = open(cds_folder + '/' + each_strain + '.cds', "a")
            cds_file_obj.write(cds_fasta)
    message = 'protein file and cds file have been created.'
    return message

def format_fasta(ana, seq, num):
    """
    格式化文本为 fasta格式
    :param ana: 注释信息
    :param seq: 序列
    :param num: 序列换行时的字符个数
    :return: fasta格式文本
    """
    format_seq = ""
    for i, char in enumerate(seq):
        format_seq += char
        if (i + 1) % num == 0:
            format_seq += "\n"
    return ana + format_seq + "\n"


def convertCDS():
    genbank_folder = os.path.join(args.indirname, 'genbank')
    cds_folder = os.path.join(args.indirname, 'cds')
    gbk_files = os.listdir(genbank_folder)
    for gbk_file in gbk_files:
        each_strain = os.path.splitext(gbk_file)[0]
        for gb_seq in SeqIO.parse(genbank_folder + '/' + gbk_file, "genbank"):
            cds_num = 1
            cds_fasta = ""
            complete_seq = str(gb_seq.seq)
            for ele in gb_seq.features:
                if ele.type == "CDS":
                    cds_seq = ""
                    cds_ana = '>' + each_strain + '-' + gb_seq.id + "-" + ele.qualifiers['protein_id'][0] + "-" + ele.qualifiers['locus_tag'][0] + " \n"
                    cds_num += 1
                    for ele1 in ele.location.parts:
                        cds_seq += complete_seq[ele1.start:ele1.end]
                    cds_fasta += format_fasta(cds_ana, cds_seq, 70)
            cds_file_obj = open(cds_folder + '/' + each_strain + '.cds', "a")
            cds_file_obj.write(cds_fasta)
    message = 'cds file have been created.'
    return message


def each_strain_pair_run(strain_pair, all_genes_dir, result_dir, strain_dict, strain_results_dir):
    """
    This function is used to call the functions above to do each strain pair genes alignment.
    :param strain_pair: each strain pair name, such as A_B or B_A
    :param all_genes_dir: a directory with all genes of each strain pair
    :param result_dir: result directory
    :param strain_dict: inherit from load_strains_label function with strain information
    :param strain_results_dir: a directory with all strain results
    :return: ......
    """
    gene_list = []
    #strain_results_dir = strain_results_collection_dir = os.path.join(args.outdirname, 'strain_pair_result')
    #all_genes_dir = strain_pair_dir = os.path.join(args.indirname, 'strain_pair_OG')
    strain_pair_genes_dir = os.path.join(all_genes_dir, strain_pair)
    for rs, ds, fs in os.walk(strain_pair_genes_dir):
        for f in fs:
            fname = os.path.splitext(f)[0]
            gene_list.append(fname)
    strain_pair_result_dir = os.path.join(result_dir, strain_pair)
    if os.path.exists(strain_pair_result_dir):
        shutil.rmtree(strain_pair_result_dir)
    os.makedirs(strain_pair_result_dir)
    tmp_gene_converted_dir = os.path.join(strain_pair_result_dir, 'tmp')
    os.makedirs(tmp_gene_converted_dir)
    pair_results = 'Clusters\tProteins\tSimilarity\tAnnotation\n'
    for gene in gene_list:
        gene_alignment_result = each_needle_run(strain_pair_genes_dir, tmp_gene_converted_dir,
                                                strain_pair_result_dir, gene, strain_dict)
        pair_results += gene_alignment_result
    shutil.rmtree(tmp_gene_converted_dir)
    pair_result_collection_file = os.path.join(strain_results_dir, '{0}.txt'.format(strain_pair))
    with open(pair_result_collection_file, 'w') as f:
        f.write(pair_results)
    logger.info('{0} is over.'.format(strain_pair))
    # time.sleep(1)   # keep process safe


def getStrainInfo():
    logger.info('strain_info file does not exist,creating from gbff file…………')
    genbank_folder = os.path.join(args.indirname, 'genbank')
    gbk_files = os.listdir(genbank_folder)
    with open(args.indirname + '/strain_info.txt', 'a')as info_file:
        info_file.write('No.' + '\t' + 'Strain' + '\t' + 'RastID')
        num = 1
        for gbk_file in gbk_files:
            id = gbk_file[0:15]
            first_record = next(SeqIO.parse(genbank_folder+'/'+gbk_file, "genbank"))
            organisam = first_record.features[0].qualifiers.get('organism')
            strain = first_record.features[0].qualifiers.get('strain')
            logger.info(gbk_file)
            if strain != None:
                if strain[0] in organisam[0]:
                    strain_name = organisam[0]
                else:
                    strain_name = organisam[0] + '_' + strain[0]
            else:
                with open(genbank_folder + '/' + gbk_file, "r", encoding='utf-8') as f:
                    definition = f.readlines()[1]
                strain_name = definition.strip().split(',')[0].split('DEFINITION  ')[1]
            info_file.write('\n' + str(num) + '\t' + strain_name + '\t' + id)
            num += 1
    message = 'strain_info.txt have been created.'
    return message

def get_ani_db_result():
    (ani_matrix, matrix_strain_dict, strain_dict) = load_ani()
    with open(args.outdirname + '/pair_HGT_result.txt', 'a')as ani_file:
        with open(args.outdirname + '/recent_HGT_results.txt', 'r')as HGT_file:
            a=0
            for line in HGT_file.readlines()[1:]:
                strain_pair_genes_dir = line.strip().split('\t')[0]
                strain_name1 = strain_pair_genes_dir.split('-')[0]
                strain_name2 = strain_pair_genes_dir.split('-')[1]
                pairs = strain_pair_genes_dir.split('-')
                pair1name = strain_dict[pairs[0]][0]
                pair2name = strain_dict[pairs[1]][0]
                pair_ani = ani_matrix[matrix_strain_dict[pair1name]][matrix_strain_dict[pair2name]]
                pair_hgt = line.strip().split('\t')[1]
                a += 1
                ani_file.writelines(str(a)+'\t' + pair1name+'\t'+pair2name+'\t'+str(pair_ani)+'\t'+str(pair_hgt)+'\t'+strain_name1+'\t'+strain_name2+'\n')
                hgt_folder = os.path.join(args.outdirname, 'csv')
                hgt_result_folder = os.path.join(args.outdirname, 'strain_pair_result')
                if not os.path.exists(hgt_folder):
                    os.makedirs(hgt_folder)
                each_distri_file = hgt_folder+'/'+strain_pair_genes_dir+'.csv'
                pair_file = hgt_result_folder+'/'+strain_pair_genes_dir+'.txt'
                data = pd.read_table(pair_file,sep='\t')
                data1 = data.iloc[:, 2]
                b = pd.cut(data1, [0, 5, 10, 15, 20, 25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100],
                           labels=[u"1", u"2", u"3", u"4", u"5",u"6", u"7", u"8", u"9", u"10",u"11", u"12", u"13", u"14", u"15",u"16", u"17", u"18", u"19", u"20"])
                c = b.value_counts()
                c2 = c.sort_index()
                d = {'Run': c2.index, 'Speed': c2.values}
                e = pd.DataFrame(d)
                e.to_csv(each_distri_file,index = False,header=True)






def first_part():
    """
    It is the first part. It is used to call pyani program to calculate ANI of each strain pair.
    :return: success message
    """
    logger.info('Part 1: Calculating average nucleotide identity (ANI) of each strain pair...')
    genbank_folder = os.path.join(args.indirname, 'genbank')
    fasta_folder = os.path.join(args.outdirname, 'fasta')
    if not os.path.exists(fasta_folder):
        os.makedirs(fasta_folder)
    i = 0
    for root, dirs, files in os.walk(genbank_folder):
        for each_file in files:
            strain_id = os.path.splitext(each_file)[0]
            gbk = os.path.join(genbank_folder, each_file)
            fasta = os.path.join(fasta_folder, strain_id + '.fasta')
            SeqIO.convert(gbk, 'genbank', fasta, 'fasta')
            i += 1
    logger.info('{0} genbank files have been converted to fasta.'.format(str(i)))
    ani_folder = os.path.join(args.outdirname, 'ANIm')
    subprocess.call(['average_nucleotide_identity.py', '-i', fasta_folder,
                     '-o', ani_folder, '-m', 'ANIm', '-g'])
    message = 'Average nucleotide identity analyses have been done.'
    retval = os.getcwd()
    indirname = os.path.join(retval, args.indirname)
    outdirname = os.path.join(retval, args.outdirname)
    genbank_folder = os.path.join(indirname, 'genbank')
    faa_folder = os.path.join(outdirname, 'faa')
    fna_folder = os.path.join(outdirname, 'fasta')
    cds_folder = os.path.join(indirname, 'cds')
    # 可以考虑把条件判断放到内部
    if not os.path.exists(faa_folder):
        os.makedirs(faa_folder)
    if not os.path.exists(cds_folder):
        os.makedirs(cds_folder)
        convertProteinandCDS(genbank_folder, faa_folder)
    # if not os.path.exists(cds_folder):
    #    os.makedirs(cds_folder)
    #    convertCDS()
    (CDS_Dict) = getCDSDict()
    strain_label_file = os.path.join(indirname, 'strain_info.txt')
    if not os.path.exists(strain_label_file):
        getStrainInfo()
    message = 'strain file convertion have been done.'
    return message


def second_part():
    """
    It is the second part. It is used to call Needle program to do pairwise sequence alignment.
    :return: success message
    """
    logger.info('Part 2: Processing pairwise sequence alignment...')
    retval = os.getcwd()
    indirname = os.path.join(retval, args.indirname)
    outdirname = os.path.join(retval, args.outdirname)
    genbank_folder = os.path.join(indirname, 'genbank')
    faa_folder = os.path.join(outdirname, 'faa')
    fna_folder = os.path.join(outdirname, 'fasta')
    cds_folder = os.path.join(indirname, 'cds')
    strain_pair_dir = os.path.join(indirname, 'strain_pair_OG')
    #if not os.path.exists(os.path.join(outdirname, "prodigal_res")):
    #    runProdigal(fna_folder)
    (CDS_Dict) = getCDSDict()
    setPairOrthofinder(faa_folder, strain_pair_dir, CDS_Dict)
    p = Pool(args.threads)
    strain_label_file = os.path.join(indirname, 'strain_info.txt')
    strain_dict = load_strains_info(strain_label_file, form=3)
    output_dir = os.path.join(outdirname, 'strain_pair_OG_alignment')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    strain_results_collection_dir = os.path.join(outdirname, 'strain_pair_result')
    if not os.path.exists(strain_results_collection_dir):
        os.makedirs(strain_results_collection_dir)
    pairs = []
    for root, dirs, files in os.walk(strain_pair_dir):
        for each_dir in dirs:
            pairs.append(each_dir)
    logger.info('Running Needle sequence alignment in {0} threads, please wait...'.format(str(args.threads)))
    for each_pair in pairs:
        p.apply_async(each_strain_pair_run, args=(each_pair, strain_pair_dir, output_dir,
                                                  strain_dict, strain_results_collection_dir))
    p.close()
    p.join()
    message = 'All alignments have been finished.'
    logger.info("Compressing/deleting %s", output_dir)
    #compress_delete_outdir(output_dir)
    return message


def third_part():
    """
    It is the third part. It is used to call R script to draw similarity distribution pictures.
    :return: success message
    """
    logger.info('Part 3: Drawing alignment distribution pictures...')
    strain_result_dir = os.path.join(args.outdirname, 'strain_pair_result')
    all_result_dir = os.path.join(args.outdirname, 'strain_result')
    if not os.path.exists(all_result_dir):
        os.makedirs(all_result_dir)
    (ani_matrix, matrix_strain_dict, strain_dict) = load_ani()

    # result_dict = defaultdict(list)
    valid_pair_list = []
    for root, dirs, files in os.walk(strain_result_dir):
        for each_file in files:
            file_name = os.path.splitext(each_file)[0]
            file_path = os.path.join(strain_result_dir, each_file)
            pairs = file_name.split('-')
            # if pairs[0] not in result_dict:
            #     result_dict[pairs[0]] = []
            pair1name = strain_dict[pairs[0]][0]
            pair2name = strain_dict[pairs[1]][0]
            pair_ani = ani_matrix[matrix_strain_dict[pair1name]][matrix_strain_dict[pair2name]]
            pair_name = '{0} ~ {1} (ANI={2}%)'.format(str(pair1name), str(pair2name), str(pair_ani))
            if pair_ani <= 100.0:
                valid_pair_list.append([pair_name, file_name, file_path])
                # with open(file_path) as f3:
                #     for line in f3.readlines()[1:]:
                #         result_line = '{0} ({1}%)\t{2}'.format(pair_name, str(pair_ani), line)
                #         result_dict[pairs[0]].append(result_line)
    r_script = os.path.join(src_dir_name, 'draw_distribution.R')
    # header_line = 'Pair\tClusters\tProteins\tSimilarity\tAnnotation\n'
    logger.info('Saving {0} pictures to {1} format.'.format(str(len(valid_pair_list)), args.gformat))
    drawing = args.drawing
    message = 'Do not draw the distributions.'
    # for each_strain, results in result_dict.items():
    #     strain_result_file = os.path.join(all_result_dir, each_strain + '.txt')
    #     with open(strain_result_file, 'w') as f4:
    #         f4.write(header_line)
    #         for each_result in results:
    #             f4.write(each_result)
    if drawing in 'tT':
        p = Pool(args.threads)
        for each_item in valid_pair_list:
            each_result_picture = os.path.join(all_result_dir, '{0}.{1}'.format(each_item[1], args.gformat))
            p.apply_async(draw_distribution, args=(each_item, each_result_picture, r_script))
        p.close()
        p.join()
        message = 'All pictures have been saved in {0}'.format(all_result_dir)
    return message


def fourth_part():
    """
    It is the fourth part. It is used to call R script to infer the number of recent HGT genes.
    :return: success message
    """
    result_dir = args.outdirname
    logger.info('Part 4: Processing recent HGT detection...')
    strain_result_dir = os.path.join(result_dir, 'strain_pair_result')
    (ani_matrix, matrix_strain_dict, strain_dict) = load_ani()
    tmp_result_dir = os.path.join(result_dir, 'tmp')
    if not os.path.exists(tmp_result_dir):
        os.makedirs(tmp_result_dir)
    else:
        shutil.rmtree(tmp_result_dir)
        os.makedirs(tmp_result_dir)
    for root, dirs, files in os.walk(strain_result_dir):
        for each_file in files:
            f_name = os.path.splitext(each_file)[0]
            f_path = os.path.join(strain_result_dir, each_file)
            pairs = f_name.split('-')
            pair1name = strain_dict[pairs[0]][0]
            pair2name = strain_dict[pairs[1]][0]
            # pair_name = str(pairs[0]) + ' ~ ' + str(pairs[1])
            pair_ani = ani_matrix[matrix_strain_dict[pair1name]][matrix_strain_dict[pair2name]]
            #if pair_ani < 94:
            tmp_file = os.path.join(tmp_result_dir, each_file)
            #if not os.path.exists(tmp_result_dir):
            if not os.path.exists(tmp_file):
                shutil.copy(f_path, tmp_file)
    r_script = os.path.join(src_dir_name, 'rHGT_alpha.R')
    param_min = 40.0
    param_max = 98.0
    devnull = open(os.devnull, 'w')
    result_file = os.path.join(result_dir, 'recent_HGT_results.txt')
    if not os.path.exists(result_file):
        try:
            subprocess.call(['Rscript', r_script, tmp_result_dir, result_file,
                             str(param_min), str(param_max)], stdout=devnull, stderr=devnull)
        except OSError:
            logger.info('Try to run {0} but failed, please check.'.format(r_script))
            logger.error(last_exception())
            sys.exit(1)
        #shutil.rmtree(tmp_result_dir)
    message = 'Recent HGT detections have been finished.'
    get_ani_db_result()
    return message


def fifth_part():
    """
    It is the fifth part. It is used to call R script to draw the comparison between the number
    of rHGTs and specific location genes (chromosome and plasmid genes) to show the accuracy.
    :return: success message
    """
    logger.info('Part 5: Drawing the comparison pictures...')
    logger.info('Loading recent HGT results.')
    genome_dir = os.path.join(args.indirname, 'genbank')
    if not os.path.exists(genome_dir):
        logger.error('These is no directory contains genbank files, please check here.')
        sys.exit(1)
    strain_id_list = []
    strain_label_file = os.path.join(args.indirname, 'strain_info.txt')
    strain_dict1 = load_strains_info(strain_label_file, form=1)
    all_gene_dict = defaultdict()
    for root, dirs, files in os.walk(genome_dir):
        for each_file in files:
            fname = os.path.splitext(each_file)
            strain_id = fname[0]                # The fname[0] is a strain id, e.g. 379.140.
            strain_id_list.append(fname[0])
            genbank = os.path.join(genome_dir, each_file)
            strain_chr_id = strain_dict1[strain_id][1]
            strain_psym_id = strain_dict1[strain_id][2]
            each_genome_dict = load_genome(genbank, strain_chr_id, strain_psym_id)
            all_gene_dict.update(each_genome_dict)
    result_dir = args.outdirname
    strain_result_dir = os.path.join(result_dir, 'strain_pair_result')
    hgt_result_file = os.path.join(args.outdirname, 'recent_HGT_results.txt')
    strain_pair_gene_dir = os.path.join(args.indirname, 'strain_pair_OG')
    detect_dict = defaultdict()
    with open(hgt_result_file, 'r') as f1:
        for each_line in f1.readlines()[1:]:
            a_list = each_line.strip().split('\t')
            detect_dict[a_list[0]] = a_list[1]  # key: strain_pair || Value: inferred HGT number
    pair_gene_location_dict = defaultdict()
    for root, dirs, files in os.walk(strain_result_dir):
        for each_file in files:
            fname = os.path.splitext(each_file)
            strain_pair = fname[0]
            pair_gene_dir = os.path.join(strain_pair_gene_dir, strain_pair)
            result_file = os.path.join(strain_result_dir, each_file)
            plasmid_gene = 0
            chr_gene = 0
            psym_gene = 0
            other_gene = 0
            with open(result_file, 'r') as f2:
                for each_line in f2.readlines()[1:]:
                    b_list = each_line.strip().split('\t')
                    identity = float(b_list[2])
                    if identity >= 98.5:
                        og_id = b_list[0]
                        og_loc_value = og_location(og_id, all_gene_dict, pair_gene_dir)
                        if og_loc_value == '**':
                            chr_gene += 1
                        elif og_loc_value == '##':
                            psym_gene += 1
                            plasmid_gene += 1
                        elif og_loc_value == '++':
                            plasmid_gene += 1
                        else:
                            other_gene += 1
            pair_gene_location_dict[strain_pair] = '{0}|{1}|{2}|{3}'.format(str(chr_gene), str(plasmid_gene),
                                                                            str(psym_gene), str(other_gene))
    tmp_result_dict = defaultdict()
    with open(hgt_result_file, 'w') as f3:
        header_line = 'strain.pair\tHGT.number\tgene.number(Chromosome|Plasmids|pSym|Ambiguity)\n'
        result_lines = ''
        for each_pair, loc_result in pair_gene_location_dict.items():
            detect_num = detect_dict.get(each_pair)
            result_lines += '{0}\t{1}\t{2}\n'.format(each_pair, detect_dict.get(each_pair), loc_result)
            loc_list = loc_result.split('|')
            tmp_result_dict[each_pair] = [detect_num, loc_list[0], loc_list[1], loc_list[2], loc_list[3]]
        f3.write(header_line + result_lines)
    tmp_order_list = ['HGT', 'Chromosome', 'Plasmids', 'pSym', 'Ambiguity']
    logger.info('Drawing comparison pictures to show detection accuracy.')
    (ani_matrix, matrix_strain_dict, strain_dict) = load_ani()
    strain_num = list(strain_dict.keys())
    strain_name = list(strain_dict.values())
    combine_result_dir = os.path.join(args.outdirname, 'combined_results')
    if not os.path.exists(combine_result_dir):
        os.makedirs(combine_result_dir)
    strain_list = []
    for each_strain_id in strain_dict1.keys():
        each_strain_name = strain_dict1[each_strain_id][0]
        strain_list.append(each_strain_name)
    r_script = os.path.join(src_dir_name, 'draw_rHGT_loc_genes.R')
    devnull = open(os.devnull, 'w')
    for query_strain in strain_list:
        query_strain_result_file = os.path.join(combine_result_dir, query_strain + '.txt')
        with open(query_strain_result_file, 'w') as f4:
            h = 'query.strain\ttype\tnumber\tani\n'
            l = ''
            for other_strain in strain_list:
                if other_strain != query_strain:
                    query_name_list = []
                    query_name_list.append(query_strain)
                    ind_query = strain_name.index(query_name_list)
                    query_num = strain_num[ind_query]
                    other_name_list = []
                    other_name_list.append(other_strain)
                    ind_other = strain_name.index(other_name_list)
                    other_num = strain_num[ind_other]
                    tmp_result_list = tmp_result_dict.get('{0}_{1}'.format(query_num, other_num))
                    if tmp_result_list is None:
                        tmp_result_list = tmp_result_dict.get('{0}_{1}'.format(other_num, query_num))
                    pair_ani = ani_matrix[matrix_strain_dict[query_strain]][matrix_strain_dict[other_strain]]
                    for i in range(len(tmp_result_list)):
                        l += '{0}\t{1}\t{2}\t{3}\n'.format(other_strain, tmp_order_list[i],
                                                           tmp_result_list[i], pair_ani)
            f4.write(h + l)
        comparison_pictures = os.path.join(combine_result_dir, '{0}.{1}'.format(query_strain, args.gformat))
        try:
            subprocess.call(['Rscript', r_script, query_strain_result_file, comparison_pictures],
                            stdout=devnull, stderr=devnull)
        except OSError:
            logger.info('Try to run {0} but failed, please check.'.format(r_script))
            logger.error(last_exception())
            sys.exit(1)
    message = 'All pictures have been saved in {0}'.format(combine_result_dir)
    return message


def auto_run():
    """
    This function is used to run all parts automatically.
    :return: success message
    """
    logger.info('Automatically run all processes...')
    message_1 = first_part()
    logger.info(message_1)
    message_2 = second_part()
    logger.info(message_2)
    message_3 = third_part()
    logger.info(message_3)
    message_4 = fourth_part()
    logger.info(message_4)
    # message_5 = fifth_part()
    # logger.info(message_5)
    done_message = 'All 4 parts have been done.'
    return done_message


def separate_run(part):
    """
    This function is used to choose one of four part to run by user.
    :param part: The part will be run.
    :return: success message
    """
    message = ''
    if part == 1:
        message = first_part()
    elif part == 2:
        message = second_part()
    elif part == 3:
        message = third_part()
    elif part == 4:
        message = fourth_part()
    elif part == 5:
        message = fifth_part()
    return message


if __name__ == '__main__':
    # Run as script
    # Parse command-line
    args = parse_cmdline()
    # Set up logging
    logger = logging.getLogger('recentHGT.py: %s' % time.asctime())
    t0 = time.time()
    src_dir_name = 'src'
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            log_stream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(log_stream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except IOError:
            logger.error("Could not open %s for logging", args.logfile)
            sys.exit(1)
    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    # Report arguments, if verbose
    # logger.info("pyani version: %s", VERSION)
    logger.info(args)
    logger.info("command-line: %s", ' '.join(sys.argv))
    # Have we got an input and output directory? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s", args.indirname)
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s", args.outdirname)
    run_message = ''
    if args.part == 0:
        try:
            run_message = auto_run()
        except OSError:
            logger.info('Try to run all 4 parts automatically but failed, please check.')
    elif 1 <= args.part <= 5:
        try:
            run_message = separate_run(args.part)
        except OSError:
            logger.info('Try to run part {0} but failed, please check.'.format(args.part))
            logger.error(last_exception())
            sys.exit(1)
    else:
        logger.error('Part {0} is not valid, please choose one of [0|1|2|3|4|5].'.format(args.part))
        sys.exit(1)
    logger.info(run_message)
    # Report that we've finished
    logger.info("All jobs have been done: %s.", time.asctime())
    logger.info("Total time taken: %.2fs", (time.time() - t0))
