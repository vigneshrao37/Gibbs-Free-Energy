#Aviran Lab Gibbs Free Energy Computation for patteRNA software
#UC Davis Department of Biomedical Engineering

import RNAstructure as r
import matplotlib.pyplot as plt
import statistics as stats
import time

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
for key,value in zip(sequences, new):
    sample_rna = key
    sample_rna = r.RNA.fromString(key)
    a = value.items()
    for i,j in a:
        sample_rna.SpecifyPair(i+1,j+1)

    score = sample_rna.CalculateFreeEnergy(1)
    gibbs_energy_scores.append(score)

print("\nScores of the given sequences are: " +"\n" + str(gibbs_energy_scores))

end = time.time()
print("\n"+"Time Taken: "+str(end - start) + " seconds")