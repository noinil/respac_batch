#!/usr/bin/env python

def main(filename):
    import numpy as np
    chg_list = []
    with open(filename, 'r') as fin:
        for lines in fin:
            words = lines.split()
            if len(words) < 1:
                continue
            charge = float(words[1])
            chg_list.append(charge)
    a = np.array(chg_list)
    print(sum(a))

if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    main(filename)
