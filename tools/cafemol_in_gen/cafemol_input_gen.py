#!/usr/bin/env python

import numpy as np
import sys

pro_resid_charge = []

def print_man():
    """Print out usage information.
    """

    man_str = """ [0;31;1m Usage: ./cafemol_input_gen.py N_BDNA PRO_NAME [0m
    So far only monomer and dimer supported.
    """
    print(man_str)

def gen_cafemol_file_part0(pro_name, dimer_opt):
    cafemol_file_name = pro_name + '.cafein'
    if dimer_opt == 1:
        cafe_template_name = "./template/pro_1_part0.txt"
    else:
        cafe_template_name = "./template/pro_2_part0.txt"
    fin_cafe = open(cafe_template_name, 'r')
    fout_cafe = open(cafemol_file_name, 'w')
    template_contents = fin_cafe.read()
    fout_cafe.write(template_contents)
    fin_cafe.close()
    fout_cafe.close()

def gen_cafemol_file_part2(pro_name, dimer_opt):
    cafemol_file_name = pro_name + '.cafein'
    if dimer_opt == 1:
        cafe_template_name = "./template/pro_1_part2.txt"
    else:
        cafe_template_name = "./template/pro_2_part2.txt"
    fin_cafe = open(cafe_template_name, 'r')
    fout_cafe = open(cafemol_file_name, 'a')
    template_contents = fin_cafe.read()
    fout_cafe.write(template_contents)
    fin_cafe.close()
    fout_cafe.close()

def gen_cafemol_file_part4(pro_name):
    cafemol_file_name = pro_name + '.cafein'
    cafe_template_name = "./template/pro_part4.txt"

    fin_cafe = open(cafe_template_name, 'r')
    fout_cafe = open(cafemol_file_name, 'a')
    template_contents = fin_cafe.read()
    fout_cafe.write(template_contents)
    fin_cafe.close()
    fout_cafe.close()

def gen_charge_modifications(n_dsDNA, pro_name):
    cg_id_DNA_shift = 6 * n_dsDNA - 2

    charge_result_name = '../../results/' + pro_name + '.charge'
    with open(charge_result_name, 'r') as fin:
        for lines in fin:
            words = lines.split()
            if len(words) < 0:
                continue
            id = int(words[0])
            ch = float(words[1])
            pro_resid_charge.append( (id, ch) )

    # -------------------- output for CafeMol input --------------------
    cafemol_file_name = pro_name + '.cafein'
    fout_cafe = open(cafemol_file_name, 'a')

    out_str = "CHARGE_CHANGE  {0:4d}   {1:8.3f}\n"
    for iv in pro_resid_charge:
        i, c = iv[0] + cg_id_DNA_shift, iv[1]
        fout_cafe.write(out_str.format(i, c))
    fout_cafe.close()

def gen_group_information(n_dsDNA, pro_name, pro_aa_num):
    cg_id_DNA_shift = 6 * n_dsDNA - 2

    cafemol_file_name = pro_name + '.cafein'
    fout_cafe = open(cafemol_file_name, 'a')

    out_str = "GROUP({0:1d})   ({1:d}-{2:d})\n"
    fout_cafe.write(out_str.format(1, 1, cg_id_DNA_shift))
    fout_cafe.write(out_str.format(2, cg_id_DNA_shift + 1, cg_id_DNA_shift + pro_aa_num))
    fout_cafe.close()

def read_pdb(pro_name):
    pdb_name = '../../pdb_protein/' + pro_name + '.pdb'
    aa_id = -1
    cg_id = 0
    ch_id = '*'
    ch_num = 0
    with open(pdb_name, 'r') as pdb_in:
        for line in pdb_in:
            if line.startswith('ATOM  '):
                pdb_resid = int(line[22:26])
                chain_id = line[21]
                if chain_id != ch_id:
                    aa_id = pdb_resid
                    ch_id = chain_id
                    cg_id += 1
                    ch_num += 1
                elif pdb_resid == aa_id + 1:
                    aa_id = pdb_resid
                    cg_id += 1
                elif pdb_resid > aa_id:
                    print(" Something like gap in backbone found! ", chain_id, ' ', aa_id)
                    return
    return {'res_num':cg_id, 'chain_num':ch_num}

if __name__ == '__main__':
    n_dsDNA = 100
    dimer_opt = 1
    pro_aa_num = 0

    if len(sys.argv) < 3:
        print_man()
        exit()
    try:
        n_dsDNA = int(sys.argv[1])
        pro_name = sys.argv[2]
    except:
        print_man()
        exit()

    pdb_info = read_pdb(pro_name)
    pro_aa_num = pdb_info['res_num']
    dimer_opt = pdb_info['chain_num']

    if dimer_opt > 2:
        print(" [0;31;1m  More than 2 chains found in protein! [0m")

    gen_cafemol_file_part0(pro_name, dimer_opt)
    gen_charge_modifications(n_dsDNA, pro_name)
    gen_cafemol_file_part2(pro_name, dimer_opt)
    gen_group_information(n_dsDNA, pro_name, pro_aa_num)
    gen_cafemol_file_part4(pro_name)
