import copy
from ens210_project_lib import prot_dict, fastareader

sequences, names = fastareader("sequences_1000_v2.fas")

seq_len = len(sequences[0])
print(seq_len)
print(len(sequences))

prob_mat = [[]]*int(seq_len)
for i in range(len(prob_mat)):
    prob_mat[i] = [0]*22

diff_seq = []

for i in range(seq_len):
    temp = []
    for j in range(len(sequences)):
        this_prot = prot_dict[sequences[j][i]]
        # if this_prot != 20:
        #     print(i,j)
        #     print(this_prot)
        prob_mat[i][this_prot] = prob_mat[i][this_prot] + 1
        if this_prot != 20:
            temp.append(j)
    if len(temp) <3:
        diff_seq.append(copy.deepcopy(temp))

problem_seq = []
for i in range(len(sequences)):
    problem_seq.append([i,0])

for row in diff_seq:
    for ele in row:
        problem_seq[ele] = problem_seq[ele][0] , problem_seq[ele][1] + 1

problem_seq.sort(reverse = True, key=lambda x: x[1])

print(problem_seq[:10]) #did not found too problematic sequence