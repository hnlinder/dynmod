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
beta = 2
q = .3

Nr = 14
tmax = 3600 * 100
dt = 1/2

kp = 1/60
kd = 1/1800
Np0 = 0
kd_rib = 1/1800
kd_mrna = 1/300
kp_mrna = 1/600

meanlength = int(tmax/2)
mean = np.empty(10)
tarr = np.linspace(0, tmax-1, int((tmax-1)/dt))

#---------------------------------------------------------------------------------------------------
# task 3 classes and funcs
class mRNA:
    def __init__(self, lattice = np.zeros(L), ribs = 0):
        self.lattice = lattice
        self.ribs = ribs

    # one time step of the ribosomes
    def rib_step(self, proteins_produced, Nr, codon_seq):
        rib_index = np.where(self.lattice>0)[0]
        np.random.shuffle(rib_index)
        for i in rib_index:
            if i != L-1:
                if ((rand() < dt * q) and valid_move(self.lattice,i)):
                    self.lattice[i+1] += 1
                    self.lattice[i] -= 1
        if ((self.lattice[0]==0) and (rand() < dt * alpha * Nr)):
            self.lattice[0] += 1
            Nr -=1
        if ((self.lattice[L-1]>0) and (rand() < dt * beta)):
            self.lattice[L-1] -= 1
            proteins_produced += 1
            # print("MADE PROTEIN")
            Nr +=1
        if (rand() < 1*dt * kd * proteins_produced):
            # print("REMOVED PROTEIN")
            proteins_produced -= 1

        return proteins_produced, Nr


class LM:
    def __init__(self, mRNA = mRNA(), matr = []):
        self.mRNA = mRNA
        self.matr = matr

    def create_mRNA(self):
        self.matr.append(mRNA(np.zeros(L)))

    def remove_mRNA(self, index):
        ribs = np.sum(self.matr[index].lattice)
        del(self.matr[index])
        return ribs


def valid_move(lattice, i):
    return (lattice[i+1] == 0) and (lattice[i] > 0)



#---------------------------------------------------------------------------------------------------
# task 3 classes and funcs

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

codon_seq = read_codons(bp_seq)
shuffled_codons = randomize_codons(codon_seq)
get_rates(shuffled_codons)

# check_rates(shuffled_codons)

#--------------------------------------------------------------------------------------------------
# main loop part
ind = 0
proteins_produced = np.zeros(len(tarr))
lm = LM()
lm.create_mRNA()
# main loop
for t in range(len(tarr)):
    proteins_produced[ind] = proteins_produced[ind - 1]
    # make new mrna
    if (1*kp_mrna * dt > rand())  :
        lm.create_mRNA()
    # Do a time step with ribosomes for each mRNA
    for index, mrna in enumerate(lm.matr):
        proteins_produced[ind], Nr = mrna.rib_step(proteins_produced[ind], Nr, shuffled_codons)
        # remove mrna
        if (1*kd_mrna * dt > rand()):
            Nr += lm.remove_mRNA(index)
        # print(mrna.lattice)
    # sleep(.1)

    # print(Nr)
    ind +=1

plt.plot(tarr, proteins_produced)
plt.show()
