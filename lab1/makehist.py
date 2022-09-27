import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep
import json
import random

# load json data
with open("rates_from_codon.json","r") as f:
    rates_from_codon = json.load(f)


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


# Codon class
class triplet:
    def __init__(self,codon, rate = 0):
        self.codon= codon
        self.rate = rate

def read_codons(bp_seq):
    codon_seq = []
    for i in range(0,len(bp_seq),3):
        codon_seq.append(triplet(bp_seq[i:i+3]))
    return codon_seq

def get_rates(codon_seq):
    for codon in codon_seq:
        codon.rate = rates_from_codon.get(codon.codon)





def conv_J_to_int(J_arr):
    for i, J in enumerate(J_arr):
        J_arr[i] =float(J.strip())
        # print(type(J_arr[i]))
        # print(J)

def mk_hist_int(J_arr):
    # histstop = .4
    # plt.figure()
    histstop = .25
    boxwidth = .02
    nrboxes = histstop/boxwidth
    bins = np.linspace(0, histstop, int(nrboxes))
    # bins = [0, .02, .04,.06]
    plt.hist(J_arr, bins, label="Original codon sequence.")
    plt.legend(loc="best")
    plt.xlabel("Protein current J [1/s]")
    plt.ylabel("Number of simulations [1]")
    # plt.show()

def mk_hist(J_arr):
    # histstop = .4
    plt.figure()
    histstop = max(J_arr)
    boxwidth = .02
    nrboxes = histstop/boxwidth
    bins = np.linspace(0, histstop, int(nrboxes))
    # bins = [0, .02, .04,.06]
    plt.hist(J_arr, bins, label="Randomized codons making the same protein.")
    # plt.show()


def get_top_rates(J_arr, codon_seqs):
    nr_runs = len(J_arr)
    inds = sorted(range(nr_runs), key=lambda x: J_arr[x])[-int(nr_runs/4):]
    codon_seq = []
    mean_rate = []
    for i in inds:
        bp_seq = codon_seqs[i]
        codon_seq.append(read_codons(bp_seq))
        get_rates(codon_seq[-1])
        # mean_rate.append(np.mean(codon_seq))
    # print(mean_rate)
    return codon_seq

def get_bottom_rates(J_arr, codon_seqs):
    nr_runs = len(J_arr)
    inds = sorted(range(nr_runs), key=lambda x: J_arr[x])[:int(nr_runs/4)]
    codon_seq = []
    mean_rate = []
    for i in inds:
        bp_seq = codon_seqs[i]
        codon_seq.append(read_codons(bp_seq))
        get_rates(codon_seq[-1])
        # mean_rate.append(np.mean(codon_seq))
    # print(mean_rate)
    return codon_seq

def get_mean_rates(rate_codon_seqs):
    rates  = np.zeros([len(rate_codon_seqs), len(rate_codon_seqs[0])])
    # print(len(rate_codon_seqs[0]))
    nr_codons = len(rate_codon_seqs[0][:-1])
    mean_rate = np.zeros([nr_codons,1])
    for i in range(nr_codons):
        for j, codon_seq in enumerate(rate_codon_seqs):
            rates[j,i] = codon_seq[i].rate
        mean_rate[i] = np.mean(rates[:,i])
    # print(mean_rate)
    return mean_rate


def plot_rates(top_mean_rates, bottom_mean_rates, og_rate):
    plt.figure()

    plt.plot(top_mean_rates, "-", label="Top performer")
    plt.plot(bottom_mean_rates,"--", label="Bottom performer")
    plt.plot(og_rate, "-.", label="Original codon sequence")
    plt.legend(loc="best")
    plt.xlabel("Index on codon sequence [1]")
    plt.ylabel("Average translation rate [1/s]")


def linfit_mean_rates(rate_codon_seqs):
    rates  = np.zeros([len(rate_codon_seqs), len(rate_codon_seqs[0])])
    # print(len(rate_codon_seqs[0]))
    nr_codons = len(rate_codon_seqs[0][:-1])
    mean_rate = np.zeros([len(rate_codon_seqs),1])
    for i, codon_seq in enumerate(rate_codon_seqs):
        for j, codon in enumerate(codon_seq[:-1]):
            print(codon.rate)
            rates[i, j] = codon.rate
        mean_rate[i] = np.mean(rates[i, :])
    # print(mean_rate)
    return mean_rate


def lin_fit(J_arr, codon_seq_arr):
    codon_seq = []
    for bp_seq in codon_seq_arr:
        codon_seq.append(read_codons(bp_seq))
        get_rates(codon_seq[-1])
    mean_rate = linfit_mean_rates(codon_seq)
    print((mean_rate))
    print(len(J_arr))
    P = np.polyfit(J_arr, mean_rate, 1)
    plt.figure()
    plt.plot(J_arr, np.polyval(P,J_arr))
    plt.xlabel("Protein current J [1/s]")
    plt.ylabel("Average translation rate for that codon sequence [1/s]")
    plt.show()



def print_top_rates(codon_seq):
    for codon in codon_seq:
        print(type(codon))
# nr_tries = 100
# codon_seq_arr,J_arr  = load_data("simdata_task4_v4.txt")
codon_seq_arr,J_arr  = load_data("simdata_task4_unshuffled_v1.txt")
codon_seq_arr_shuff,J_arr_shuff  = load_data("simdata_task4_with_alpha_v2.txt")

top_rate_seq = get_top_rates(J_arr_shuff, codon_seq_arr_shuff)
top_mean_rates = get_mean_rates(top_rate_seq)
bottom_rate_seq = get_bottom_rates(J_arr_shuff, codon_seq_arr_shuff)
bottom_mean_rates=  get_mean_rates(bottom_rate_seq)

# print(codon_seq_arr[0])
og_rate = read_codons(codon_seq_arr[0])
get_rates(og_rate)
# print(og_rate[0].rate)
og_rate = get_mean_rates([og_rate])
plot_rates(top_mean_rates, bottom_mean_rates, og_rate)


conv_J_to_int(J_arr)
conv_J_to_int(J_arr_shuff)
lin_fit(J_arr_shuff, codon_seq_arr_shuff)
mk_hist(J_arr_shuff)
mk_hist_int(np.mean(J_arr))
plt.show()
