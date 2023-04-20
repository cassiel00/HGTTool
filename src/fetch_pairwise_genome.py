#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to fetch pairwise core genomes between each strain pair.
#               It is better to put this script in ITEP directory.
# Created by galaxy on 2016/10/21 15:11


import os
import sys
import shutil
from collections import defaultdict
from itertools import combinations


def replace(cur_dir):
    old_id = "all_I_2.0_c_{0}_m_maxbit_".format(maxbit)
    new_id = ""
    for parent, dirnames, filenames in os.walk(cur_dir):
        for filename in filenames:
            if filename.find(old_id) != -1:
                new_name = filename.replace(old_id, new_id)
                # print(filename, "---->", newName)
                os.rename(os.path.join(parent, filename), os.path.join(parent, new_name))


my_path = os.getcwd()
maxbit = 0.4
# source_file = os.path.join(my_path, 'SourceMe.sh')
# os.system('source {0}'.format(source_file))
strain_information_file = os.path.join(my_path, 'strain_info.txt')
strain_dict = defaultdict()
strain_list = []
with open(strain_information_file, 'r') as f1:
    for each_line in f1.readlines()[1:]:
        a_list = each_line.strip().split('\t')
        strain_dict[a_list[1]] = [a_list[2]]
        strain_list.append(a_list[1])
# fetch_strain_pair = []
# for i in strain_list:
#     for j in strain_list:
#         if i != j:
#             fetch_strain_pair.append((i, j))
fetch_strain_pair = list(combinations(strain_list, 2))
strain_pair_dir = os.path.join(my_path, 'all_strain_pairs_{0}'.format(maxbit))
if not os.path.exists(strain_pair_dir):
    os.makedirs(strain_pair_dir)
else:
    shutil.rmtree(strain_pair_dir)
    os.makedirs(strain_pair_dir)
tmp_organisms = os.path.join(my_path, 'tmp_organisms.txt')
replace_file = os.path.join(my_path, 'replace.py')
for each_pair in fetch_strain_pair:
    the_strain = each_pair[0]
    other_strain = each_pair[1]
    # if strain_dict[the_strain][2] != strain_dict[other_strain][2]:
    # if strain_dict[the_strain][1] != strain_dict[other_strain][1]:
    the_strain_name = the_strain.split(' ')[-1]
    other_strain_name = other_strain.split(' ')[-1]
    strain_gene_dir = os.path.join(strain_pair_dir, '{0}_{1}'.format(the_strain_name, other_strain_name))
    if not os.path.exists(strain_gene_dir):
        the_strain_line = '{0}\t{1}\n'.format(
            the_strain, strain_dict[the_strain][0])
        other_strain_line = '{0}\t{1}\n'.format(
            other_strain, strain_dict[other_strain][0])
        with open(tmp_organisms, 'w') as f2:
            result_line = the_strain_line + other_strain_line
            f2.write(result_line)
        cmd = 'cat {0} | python {1}/src/db_findClustersByOrganismList.py -a -u all_I_2.0_c_{3}_m_maxbit | python ' \
              '{1}/src/db_getClusterGeneInformation.py|grep -F -f {0} |python ' \
              '{1}/src/getClusterFastas.py -n {2}'.format(
                  tmp_organisms, my_path, strain_gene_dir, maxbit)
        os.system(cmd)
        replace(strain_gene_dir)
        os.remove(tmp_organisms)
os.system('mv all_strain_pairs_{0} all_strain_pairs'.format(maxbit))
