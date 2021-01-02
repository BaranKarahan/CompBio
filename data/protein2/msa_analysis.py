import copy
from ens210_project_lib import *


sequences, names = fastareader("sequences_1000.fas")
distances = []
distances2 = []
for i in range(len(sequences[0])):
    p1 = sequences[0][i]
    dist1 = 0
    dist2 = 0
    if p1 != '-':
        for j in range(1,len(sequences)):
            p2 = sequences[j][i]
            if p1 != p2:
                if p2 != '-':
                    dist1 = dist1 + 1
            m1,m2 = min(prot_dict[p1],prot_dict[p2]),max(prot_dict[p1],prot_dict[p2])
            dist2 = dist2 + blosum62_mat[m1][m2]

    else:
        dist1 = -1
    if dist1 != -1:
        distances.append(dist1)
        distances2.append(dist2)
print(distances2)

