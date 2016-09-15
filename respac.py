#!/usr/bin/env python

import os
import re

def main(pro_name):
    '''Perform RESPAC calculations for given structure of a protein.
    pro_name:  name of protein.
    '''

    # -------------------- Commands & Params --------------------
    pqr_command_path = "/home/noinil/Workspace/respac/pdb2pqr-linux-bin64-2.0.0/pdb2pqr "
    apbs_command_path = "/home/noinil/Workspace/respac/APBS-1.4.1-binary/bin/apbs "
    dxmath_command_path = "/home/noinil/Workspace/respac/APBS-1.4.1-binary/share/apbs/tools/bin/dxmath "
    surface_command_path = "/home/noinil/Workspace/respac/surface/bin/surface "
    pdcp_command_path = "/home/noinil/Workspace/respac/pdcp/bin/pdcp"

    shell_path = "export LD_LIBRARY_PATH=~/Workspace/respac/APBS-1.4.1-binary/lib/:/usr/local/lib/:/usr/lib/:$LD_LIBRARY_PATH; "

    # -------------------- Basic Variables --------------------
    ionic_strength  = 0.001
    apbs_box_margin = 20.0
    apbs_grid_size  = 0.45
    apbs_radius_A   = 3.0
    apbs_radius_B   = 12.0

    # don't modify the following...
    length_x, length_y, length_z = 0, 0, 0
    debye_length = 100.0

    # -------------------- File Names --------------------
    pdb_name = './pdb_protein/' + pro_name + '.pdb'

    if not os.path.exists(pdb_name):
        print(" !!! ERROR:", pdb_name, " not found. ")
        return

    os.system('mkdir -p run/pqr run/apbs_in run/apbs_out run/surf_in run/pdc_in run/pdb')
    pdb_tmp_name = './run/pdb/' + pro_name + '.pdb'
    pqr_name = './run/pqr/' + pro_name + '.pqr'
    apbs_name = './run/apbs_in/' + pro_name + ".in"
    apbs_vol_A_name = './run/apbs_in/' + pro_name + "_vol_A.in"
    apbs_vol_B_name = './run/apbs_in/' + pro_name + "_vol_B.in"
    apbs_out_name = './run/apbs_out/' + pro_name + '_apbs_potential.dx'
    volm_out_name = './run/apbs_out/' + pro_name + '_delta_volm.dx'
    surf_name = './run/surf_in/' + pro_name + '.surf'
    pdc_name = './run/pdc_in/' + pro_name + '.pdcin'
    charge_name = './results/' + pro_name + '.charge'

    # -------------------- reading PDB --------------------
    print('============================================================')
    print(' Processing PDB file:', pdb_name)
    aa_id = -1
    ch_id = '*'
    cg_id = 0
    fout_pdb_tmp = open(pdb_tmp_name, 'w')
    with open(pdb_name, 'r') as pdb_in:
        coor_x = []
        coor_y = []
        coor_z = []
        for line in pdb_in:
            if line.startswith('ATOM  ') or line.startswith('HETATM'):
                pdb_resid = int(line[22:26])
                chain_id = line[21]
                res_name = line[17:20]
                if res_name == 'ZN ' or res_name == ' ZN':
                    res_name = 'ZN2'
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coor_x.append(x)
                coor_y.append(y)
                coor_z.append(z)
                if pdb_resid > aa_id or chain_id != ch_id:
                    aa_id = pdb_resid
                    ch_id = chain_id
                    cg_id += 1
                newline = 'ATOM  ' + line[6:17] + res_name + line[20:22] + str(cg_id).rjust(4) + line.rstrip()[26:] + '\n'
                fout_pdb_tmp.write(newline)
        fout_pdb_tmp.write('END')
        x_min, x_max = min(coor_x), max(coor_x)
        y_min, y_max = min(coor_y), max(coor_y)
        z_min, z_max = min(coor_z), max(coor_z)
        length_x = round(x_max - x_min) + 2.0
        length_y = round(y_max - y_min) + 2.0
        length_z = round(z_max - z_min) + 2.0

    print(' Getting box size:', length_x, '*', length_y, '*', length_z)
    fout_pdb_tmp.close()

    # -------------------- Generating PQR file --------------------
    print('============================================================')
    print(' Generating pqr file:', pqr_name)
    pqr_command_params = "--ff=CHARMM " + pdb_tmp_name + " " + pqr_name + " >./run/PDB2PQR.log 2>&1"
    pqr_command = pqr_command_path + " " + pqr_command_params
    try:
        os.system(pqr_command)
    except:
        print(" !!! ERROR: pdb2pqr failed!")
        return
    print(' Done... ')

    # -------------------- Generating APBS files --------------------
    print('============================================================')
    print(' Generating APBS files for:', pqr_name)
    apbs_box_xlen = length_x + apbs_box_margin
    apbs_box_ylen = length_y + apbs_box_margin
    apbs_box_zlen = length_z + apbs_box_margin
    apbs_grid_n_x = int(apbs_box_xlen / apbs_grid_size)
    apbs_grid_n_y = int(apbs_box_ylen / apbs_grid_size)
    apbs_grid_n_z = int(apbs_box_xlen / apbs_grid_size)

    fout_apbs_in = open(apbs_name, 'w')
    with open('./lib/template/apbs_in_template') as fin_apbs_temp:
        for line in fin_apbs_temp:
            new_line = line
            new_line = re.sub('PQRFILE', str(pqr_name), new_line)
            new_line = re.sub('DIMX', str(apbs_grid_n_x), new_line)
            new_line = re.sub('DIMY', str(apbs_grid_n_y), new_line)
            new_line = re.sub('DIMZ', str(apbs_grid_n_z), new_line)
            new_line = re.sub('BOXLX', str(apbs_box_xlen), new_line)
            new_line = re.sub('BOXLY', str(apbs_box_ylen), new_line)
            new_line = re.sub('BOXLZ', str(apbs_box_zlen), new_line)
            new_line = re.sub('IONIC_STRENGTH', str(ionic_strength), new_line)
            fout_apbs_in.write(new_line)
    fout_apbs_in.close()
    print(' APBS input file ', apbs_name, ' created.')

    fout_apbs_vol_A_in = open(apbs_vol_A_name, 'w')
    with open('./lib/template/apbs_vol_in_template') as fin_apbs_vol_temp:
        for line in fin_apbs_vol_temp:
            new_line = line
            new_line = re.sub('PQRFILE', str(pqr_name), new_line)
            new_line = re.sub('DIMX', str(apbs_grid_n_x), new_line)
            new_line = re.sub('DIMY', str(apbs_grid_n_y), new_line)
            new_line = re.sub('DIMZ', str(apbs_grid_n_z), new_line)
            new_line = re.sub('BOXLX', str(apbs_box_xlen), new_line)
            new_line = re.sub('BOXLY', str(apbs_box_ylen), new_line)
            new_line = re.sub('BOXLZ', str(apbs_box_zlen), new_line)
            new_line = re.sub('RADIUS', str(apbs_radius_A), new_line)
            new_line = re.sub('OUTPUT', 'vol_A', new_line)
            new_line = re.sub('IONIC_STRENGTH', str(ionic_strength), new_line)
            fout_apbs_vol_A_in.write(new_line)
    fout_apbs_vol_A_in.close()
    print(' APBS input file ', apbs_vol_A_name, ' created.')

    fout_apbs_vol_B_in = open(apbs_vol_B_name, 'w')
    with open('./lib/template/apbs_vol_in_template') as fin_apbs_vol_temp:
        for line in fin_apbs_vol_temp:
            new_line = line
            new_line = re.sub('PQRFILE', str(pqr_name), new_line)
            new_line = re.sub('DIMX', str(apbs_grid_n_x), new_line)
            new_line = re.sub('DIMY', str(apbs_grid_n_y), new_line)
            new_line = re.sub('DIMZ', str(apbs_grid_n_z), new_line)
            new_line = re.sub('BOXLX', str(apbs_box_xlen), new_line)
            new_line = re.sub('BOXLY', str(apbs_box_ylen), new_line)
            new_line = re.sub('BOXLZ', str(apbs_box_zlen), new_line)
            new_line = re.sub('RADIUS', str(apbs_radius_B), new_line)
            new_line = re.sub('OUTPUT', 'vol_B', new_line)
            new_line = re.sub('IONIC_STRENGTH', str(ionic_strength), new_line)
            fout_apbs_vol_B_in.write(new_line)
    fout_apbs_vol_B_in.close()
    print(' APBS input file ', apbs_vol_B_name, ' created.')

    # -------------------- Calc Ele Potential --------------------
    print('============================================================')
    print(' Calculating electrostatic potential... ')
    print(' Step 1 of 3: APBS potentials...')
    apbs_command = shell_path + apbs_command_path + " " + apbs_name + " >./run/APBS1.log 2>&1"
    try:
        os.system(apbs_command)
    except:
        print(" !!! ERROR: APBS calculation 1 failed!")
        return
    print(' Done... ')

    print(' Step 2 of 3: APBS volume A...')
    apbs_command = shell_path + apbs_command_path + " " + apbs_vol_A_name + " >./run/APBS2.log 2>&1"
    try:
        os.system(apbs_command)
    except:
        print(" !!! ERROR: APBS calculation 2 failed!")
        return
    print(' Done... ')

    print(' Step 3 of 3: APBS volume B...')
    apbs_command = shell_path + apbs_command_path + " " + apbs_vol_B_name + " >./run/APBS3.log 2>&1"
    try:
        os.system(apbs_command)
    except:
        print(" !!! ERROR: APBS calculation 3 failed!")
        return
    print(' Done... ')

    print(' DXMATH calculating...')
    dxmath_command = shell_path + dxmath_command_path + " ./lib/template/dxmath.inp" + " >./run/DXMATH.log 2>&1"
    try:
        os.system(dxmath_command)
    except:
        print(" !!! ERROR: dxmath calculation failed!")
        return
    print(' Done... ')

    mv_apbs_out_command = "mv apbs_potential.dx " + apbs_out_name
    mv_volm_out_command = "mv delta_vol.dx " + volm_out_name
    try:
        os.system(mv_apbs_out_command)
        os.system(mv_volm_out_command)
        os.system("rm -f vol_*.dx")
    except:
        print(" Something wrong with APBS claculations...")
        return

    # -------------------- Determine surface --------------------
    print('============================================================')
    print(' Determining surface residues... ')
    surface_params1 = " --pqr " + pqr_name + " --ofname " + surf_name
    surface_params2 = " --dbox 6.0 --r_probe 4.0 "  + " >./run/SURFACE.log 2>&1"
    surface_command = shell_path + surface_command_path + surface_params1 + surface_params2
    try:
        os.system(surface_command)
    except:
        print(" !!! ERROR: Program surface failed!")
    print(" Done...")

    # -------------------- Calc PDC --------------------
    print('============================================================')
    print(' Computing Potential Derived Charges (PDC)...')

    # reading Debye Length
    with open('io.mc', 'r') as fin_iomc:
        for line in fin_iomc:
            if "Debye length =" in line:
                words = line.split()
                debye_length = float(words[4])
                break
    os.system("rm -f io.mc")
    print(" Detected Debye Length = ", debye_length)

    fout_pdc_in = open(pdc_name, 'w')
    with open('./lib/template/pdc_in_template') as fin_pdc_temp:
        for line in fin_pdc_temp:
            new_line = line
            new_line = re.sub('DEBYE', str(debye_length), new_line)
            fout_pdc_in.write(new_line)
    fout_pdc_in.close()
    print(' PDC input file: ', pdc_name, ' created.')

    pdcp_params1 = ' --ifname ' + pdc_name + ' --pqr ' + pqr_name
    pdcp_params2 = ' --pot ' + apbs_out_name + ' --vol ' + volm_out_name
    pdcp_params3 = ' --site All ' + ' --residue ' + surf_name
    pdcp_params4 = ' --ofname ' + charge_name  + " >./run/RESPAC.log 2>&1"
    pdcp_command = pdcp_command_path + pdcp_params1 + pdcp_params2 + pdcp_params3 + pdcp_params4
    try:
        os.system(pdcp_command)
    except:
        print(" !!! ERROR: Program pdcp failed!")

    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("[1;32m Calculation finished! [0m Please see results in ", charge_name)
    print(" Additional results are also provided in [4;36m tools.[0m")


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
