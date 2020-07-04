#Aviran Lab Gibbs Free Energy Computation


#import RNAstructure as r
import time
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
    #print(parentheses_locs)
    return parentheses_locs

new = []
for i, s in enumerate(motifs):
    #print('\nmotif {}:\n{}'.format(i, s))
    parentheses_locs = find_parentheses(s)
    new.append(parentheses_locs)

# This code returns the integers that represents the index of pairs

for key,value in zip(sequences, new):
    print(key)
    a = value.items()
    for i,j in a:
        print(i+1,j+1)

end = time.time()
print("\n"+"Time Taken: "+str(end - start) + " seconds")