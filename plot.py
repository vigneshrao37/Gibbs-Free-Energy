#Aviran Lab Gibbs Free Energy Computation for patteRNA software
#UC Davis Department of Biomedical Engineering

#from AviranJuly2 import sequences
import RNAstructure as r
import time
import matplotlib.pyplot as plt
import numpy as np

sequences = []
with open("score.txt","r") as f:
    line = f.readlines()
    for l in line[1:]:
        sequences.append(l.split()[-1])

times = []
lenghts  = []
for i in sequences:
    sample = r.RNA.fromString(i)
    lenghts.append(len(i))

    #timing Partition Function
    start = time.time()
    sample.PartitionFunction()
    end = time.time()
    #timing Partiton Function

    times.append(end - start)

times1 = []
for i in sequences:
    sample1 = r.RNA.fromString(i)

    #timing Linear Fold
    start1 = time.time()
    sample1.FoldSingleStrand()
    end1 = time.time()
    #timing Linear Fold

    times1.append(end1-start1)

if __name__ == "__main__":
    
    plt.plot(lenghts,times,"-r",label = "Partition Function")
    plt.plot(lenghts,times1,"-b",label = "Linear Fold Function")
    plt.xlabel("Lenghts of RNA samples (nucleotides)")
    plt.ylabel("Time Needed for Execution (seconds)")
    plt.title("Lenghts of RNA versus Time of Execution")
    plt.legend()
    plt.show()
    