import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep
import json
import random

bp_seq ="AGCGCGCGGUCACAACGUUACUGUUAUCGAUCCGGUCGAAAAACUGCUGGCAGUGGGGCAUUACCUCGAAUCUACCGUCGAUAUUGCUGA"

qi = 1



L  = 30
alpha = .9
# beta = 2
# q = .3

Nr = 4
tmax = 3600 * 100
dt = 1/2

kp = 1/600
kd = 1/1800
Np0 = 0
kd_rib = 1/1800
kd_mrna = 1/300
kp_mrna = 1/600

meanlength = int(tmax/2)
# mean = np.empty(10)
tarr = np.linspace(0, tmax-1, int((tmax-1)/dt))

#---------------------------------------------------------------------------------------------------
# task 3 classes and funcs
class mRNA:
    def __init__(self, lattice = np.zeros(L), ribs = 0):
        self.lattice = lattice
        self.ribs = ribs

    # one time step of the ribosomes
    def rib_step(self, proteins_produced, Nr, codon_seq):
        # alpha = codon_seq[0].rate
        beta = codon_seq[-1].rate

        rib_index = np.where(self.lattice>0)[0]
        np.random.shuffle(rib_index)
        occupied = 0
        for i in rib_index:
            if i != L-1:
                if ((rand() < dt * codon_seq[i].rate) and valid_move(self.lattice,i)):
                    self.lattice[i+1] += 1
                    self.lattice[i] -= 1
        if ((self.lattice[0]==0) and (rand() < dt * alpha * Nr)):
            self.lattice[0] += 1
            Nr -=1
        if (self.lattice[L-1]>0):
            occupied += 1
            if (rand() < dt * beta):
                self.lattice[L-1] -= 1
                proteins_produced += 1
                # print("MADE PROTEIN")
                Nr +=1
        if (rand() < 1*dt * kd * proteins_produced):
            # print("REMOVED PROTEIN")
            proteins_produced -= 1

        return proteins_produced, Nr, occupied


class LM:
    def __init__(self, mRNA = mRNA(), matr = [], nr_mRNAs = np.zeros(len(tarr))):
        self.mRNA = mRNA
        self.matr = matr
        self.nr_mRNAs = nr_mRNAs

    def create_mRNA(self, t):
        self.matr.append(mRNA(np.zeros(L)))
        self.nr_mRNAs[t:] += 1

    def remove_mRNA(self, index, t):
        ribs = np.sum(self.matr[index].lattice)
        del(self.matr[index])
        self.nr_mRNAs[t:] -= 1
        return ribs


def valid_move(lattice, i):
    return (lattice[i+1] == 0) and (lattice[i] > 0)



#---------------------------------------------------------------------------------------------------
# task 4 classes and funcs

# Codon class
class triplet:
    def __init__(self,codon, rate = 0):
        self.codon= codon
        self.rate = rate

# load json data
with open("rates_from_codon.json","r") as f:
    rates_from_codon = json.load(f)
with open("aminoacid_from_codon.json","r") as f:
    aminoacid_from_codon = json.load(f)
with open("codons_from_aminoacid.json","r") as f:
    codons_from_aminoacid= json.load(f)

def read_codons(bp_seq):
    codon_seq = []
    for i in range(0,len(bp_seq),3):
        codon_seq.append(triplet(bp_seq[i:i+3]))
    return codon_seq


def randomize_codons(sequence):
    alt_codons = []
    for index, codon in enumerate(sequence):
        AA = aminoacid_from_codon.get(codon.codon)
        cods = codons_from_aminoacid.get(AA)
        np.random.shuffle(cods)
        alt_codons.append(triplet(cods[0]))
    return alt_codons


def get_AA_seq(codon_seq):
    AA_seq = []
    for i ,codon in enumerate(codon_seq):
        AA_seq.append(aminoacid_from_codon.get(codon.codon))
    return AA_seq


def check_shuffle(bp_seq):
    codon_seq = read_codons(bp_seq)
    unshuffled_AA_seq = get_AA_seq(codon_seq)
    shuffled_codons = randomize_codons(codon_seq)
    shuffled_AA_seq = get_AA_seq(shuffled_codons)

    for i in range(len(codon_seq)):
        print(f"""Unshuffled codseq: {codon_seq[i].codon}\n
                    shuffled codseq: {shuffled_codons[i].codon}\n
                    Are equal      : {codon_seq[i].codon==shuffled_codons[i].codon}""")

    print(f"""Unshuffled AA seq: {unshuffled_AA_seq}\n
                shuffled AA seq: {shuffled_AA_seq}\n
                Are equal      : {unshuffled_AA_seq == shuffled_AA_seq}""")

def get_rates(codon_seq):
    for codon in codon_seq:
        codon.rate = rates_from_codon.get(codon.codon)

def check_rates(codon_seq):
    for codon in codon_seq:
        print(codon.rate)

def write_codonseq(codon_seq, filename):
    codstr = ""
    ratestr  = ""
    for codon in codon_seq:
        codstr = codstr + str(codon.codon)
        ratestr = ratestr +","+ str(codon.rate)

    # np.savetxt("simdata_task4.txt", codstr)
    with open(filename, "a+") as f:
        f.write(f"{codstr}\n")
        # f.write(f"{ratestr}\n")




codon_seq = read_codons(bp_seq)
# shuffled_codons = randomize_codons(codon_seq)
# get_rates(shuffled_codons)

# check_rates(shuffled_codons)

#--------------------------------------------------------------------------------------------------
# main loop part
nr_tries = 100
filename = "simdata_task4_with_alpha_v4.txt"
occupied = np.zeros(nr_tries)
proteins_produced = np.zeros([nr_tries, len(tarr)])
mrna_prod_decay = 1
J = np.zeros(nr_tries)
# lm.create_mRNA(0)
lm = []
# main loop
for n in range(nr_tries):
    get_rates(codon_seq)
    # initiate lattice matrix
    # print(lm[n].nr_mRNAs)
    # print(lm[n].lattice)
    lm.append(LM())
    ind = 0
    mrnas = np.zeros(len(tarr))
    for t in range(len(tarr)):
        if t != 0:
            mrnas[t] = mrnas[t-1]
        proteins_produced[n, ind] = proteins_produced[n, ind - 1]
        # make new mrna
        if (mrna_prod_decay*kp_mrna * dt > rand())  :
            lm[n].create_mRNA(t)
            mrnas[t] += 1
        # Do a time step with ribosomes for each mRNA
        for index, mrna in enumerate(lm[n].matr):
            proteins_produced[n ,ind], Nr, end_occupied = mrna.rib_step(proteins_produced[n, ind], Nr,codon_seq)
            occupied[n] += end_occupied
            # remove mrna
            if ((mrna_prod_decay*kd_mrna * dt * mrnas[t] > rand()) ):
                Nr += lm[n].remove_mRNA(index, t)
                mrnas[t] -= 1

        ind +=1
    # plt.figure()
    # # plt.ylim([0,130])
    # plt.plot(tarr, proteins_produced[n,:])

    J[n] = occupied[n]/len(tarr)*codon_seq[-1].rate
    print(np.mean(mrnas))
    print(np.mean(lm[n].nr_mRNAs))
    print(J)
    write_codonseq(codon_seq, filename)
    codon_seq= randomize_codons(codon_seq)

with open(filename, "a+") as f:
    np.savetxt(f,J)
# np.savetxt("simdata_task4.txt",proteins_produced)

# print(f"nr_mRNAs : {lm[n].nr_mRNAs}\nAverage nr of mRNAs : {np.mean(lm[n].nr_mRNAs)}")
# np.savetxt("nr_mrnas.txt", lm.nr_mRNAs)
# initiate lattice matrix
# lm = LM()
# plt.plot(tarr, proteins_produced[0,:])
# plt.plot(tarr, proteins_produced[1,:])
plt.show()
