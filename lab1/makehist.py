import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep
import json
import random



def load_data(filename):
    codon_seqs = []
    J_arr = []
    with open(filename, "r") as f:
        lines = f.readlines()
    nr_tries= len(lines)/2
    for n, line in enumerate(lines):
        if n < nr_tries:
            codon_seqs.append(line)
        else:
            J_arr.append(line)
    return codon_seqs, J_arr


def mk_hist(J_arr):
    # histstop = .4
    histstop = max(J_arr)
    boxwidth = .02
    nrboxes = histstop/boxwidth
    bins = np.linspace(0, histstop, int(nrboxes))
    # bins = [0, .02, .04,.06]
    plt.hist(J_arr, bins)
    plt.show()


def conv_J_to_int(J_arr):
    for i, J in enumerate(J_arr):
        J_arr[i] =float(J.strip())
        # print(type(J_arr[i]))
        # print(J)
# nr_tries = 100
# codon_seq_arr,J_arr  = load_data("simdata_task4_v4.txt")
codon_seq_arr,J_arr  = load_data("simdata_task4_unshuffled_v1.txt")
conv_J_to_int(J_arr)
mk_hist(J_arr)
print(f"codons : {codon_seq_arr}\nJs : {J_arr}")
