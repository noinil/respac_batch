#!/usr/bin/env python

def main(seq_fname, charge_fname):
    import numpy as np
    import matplotlib.pyplot as plt

    resid = []
    charge = []

    with open(charge_fname, 'r') as fin:
        for lines in fin:
            words = lines.split()
            if len(words) < 0:
                continue
            id = int(words[0])
            ch = float(words[1])
            resid.append(id)
            charge.append(ch)

    with open(seq_fname, 'r') as fseq:
        for lines in fseq:
            if not lines.startswith('>'):
                pro_seq = lines
                break

    pro_len = len(pro_seq)
    X = np.array(resid)
    Y = np.array(charge)
    net_charge = sum(Y)
    width=0.45

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(X, Y, width, color='g')
    ax.axhline(y=0, linewidth=0.5, alpha=0.3)
    ax.grid(color='gray', alpha=0.7, ls='--', axis='y')
    ax.set_xlim(0, pro_len + 1)
    ax.set_ylim(-1.8, 1.8)
    ax.set_xlabel('Residue Index', size=20)
    ax.set_ylabel('Charge', size=20)
    ax.set_xticks([i * 20 + 20 for i in range(int(pro_len + 1) // 20)])
    ax.set_xticklabels([str(i*20+20) for i in range(int(pro_len + 1) // 20 )], fontsize=14, family='serif')
    ax.set_yticks([i -1 for i in range(3)])
    ax.set_yticklabels([r'$-1$',r'$0$',r'$1$'], fontsize=16)

    for i, v in enumerate(resid):
        resname = pro_seq[v - 1]
        if resname in ['D', 'E', 'K', 'R']:
            c = 'r' if resname in ['D', 'E'] else 'c'
            x = v - 0.2
            y = charge[i] + 0.05 if charge[i] > 0 else charge[i] - 0.04
            ax.text(x, y, r'$\circ$', ha='center', va='center', fontsize=12, color=c)

    ax.text(0.6*pro_len, -1.6, "Net Charge =" + str(net_charge), ha='center', va='bottom', fontsize=18)

    plt.savefig("_charge_distribution.png", dpi=150)
    # plt.show()


if __name__ == '__main__':
    import sys
    if sys.argv[3] != '-c' or sys.argv[1] != '-s':
        print(' Usage: ', sys.argv[0], ' -s seq.fasta -c charge.dat')
    main(sys.argv[2], sys.argv[4])
