#!/usr/bin/env python

def main():
    import numpy as np
    charge_list = []
    with open('cafemol.insert', 'r') as fin:
        for lines in fin:
            words = lines.split()
            if len(words) < 1:
                continue
            c = float(words[2])
            charge_list.append(c)
    charge_arr = np.array(charge_list)
    print(sum(charge_arr))

if __name__ == '__main__':
    main()
