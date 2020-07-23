#Aviran Lab for patteRNA software using RNAstructure
#UC Davis Department of Biomedical Engineering

import RNAstructure as r
import matplotlib.pyplot as plt
import statistics as stats
import time
from numpy import log 

# Time how long it takes to execute the program

start = time.time()

# This code reads all sequences from score file and appends data to sequences array

sequences = []
with open("score.txt","r") as f:
    line = f.readlines()
    for l in line[1:]:
        sequences.append(l.split()[-1])

# This code reads all motifs from file and appends data to motifs array

motifs = []
with open("score.txt","r") as file:
    lines = file.readlines()
    for i in lines[1:]:
        motifs.append(i.split()[-3])

def find_parentheses(s: str):
    stack = []
    parentheses_locs = {}
    for i, c in enumerate(s):
        if c == '(':
            stack.append(i)
        elif c == ')':
            parentheses_locs[stack.pop()] = i

    return parentheses_locs

new = []
for i, s in enumerate(motifs):
    parentheses_locs = find_parentheses(s)
    new.append(parentheses_locs)

# This code returns the integers that represents the index of pairs

gibbs_energy_scores = []
probs = []

for key,value in zip(sequences, new):
    sample_rna = key
    sample_rna = r.RNA.fromString(key)
    sample_rna.FoldSingleStrand()
    a = value.items()

    for i,j in a:
        sample_rna.PartitionFunction()
        sample_rna.SpecifyPair(i+1,j+1)
        prob = sample_rna.GetPairProbability(i+1,j+1)
        probs.append(prob)

    score = sample_rna.CalculateFreeEnergy(1)
    gibbs_energy_scores.append(score)

 #Computing Pseudo Enegry scores    

# def pseudo_energy(num):
#     pseudos = []
#     for i in whatever:
#         pseudos.append(pseudo_score)
#     shape_reactivity = .98
#     y_int = -.1212
#     slope = .1212
#     pseudo_score =  slope * [log(shape_reactivity) + 1] + y_int

combined = []
for i,j in zip(gibbs_energy_scores,pseudos):
    combined.append(i+j)    


end = time.time()

if __name__ == "__main__":

    print("\nScores of the given sequences are: " +"\n" + str(gibbs_energy_scores))
    print("\nProbabilities of the base pairs are : " + "\n" + str(probs))
    print("\n"+"Time Taken: "+str(end - start) + " seconds")
