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

if __name__ == "__main__":
    print(lenghts)
    print(times)
    plt.plot(lenghts,times,"-r",label = "Partition Function")
    #plt.plot(times,"-b",label = "Linear Fold")
    plt.xlabel("Lenghts of RNA samples (nucleotides)")
    plt.ylabel("Time Needed for Execution (seconds)")
    plt.title("Lenghts of RNA versus Time of Execution")
    plt.legend()
    plt.show()