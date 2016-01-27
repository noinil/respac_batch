#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

def main(filename):

    resid = []
    charge = []

    with open(filename, 'r') as fin:
        for lines in fin:
            words = lines.split()
            if len(words) < 0:
                continue
            id = int(words[0])
            ch = float(words[1])
            resid.append(id)
            charge.append(ch)

    # -------------------- output for CafeMol input --------------------
    def cafe_out(dna_length):
        fout = open('cafemol.insert', 'w')
        out_str = "CHARGE_CHANGE  {0:4d}   {1:8.3f}\n"
        for i, v in enumerate(resid):
            fout.write(out_str.format(v + dna_length, charge[i]))
        fout.close()

    # ---------- for HU specifically ----------
    cafe_out(0)

if __name__ == '__main__':
    main(sys.argv[1])
