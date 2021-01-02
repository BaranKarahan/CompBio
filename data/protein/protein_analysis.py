import copy
from ens210_project_lib import prot_dict, fastareader

sequences, names = fastareader("sequences_1000.fas")

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

print(problem_seq[:5]) #found 3 sequences where this value is bigger than or equal to 80. I concluded that printing 10 of the biggest would show the result.
    
to_be_erased = []

for i in range(3):
    to_be_erased.append(problem_seq[i][0])

new_sequences = []
for ele in sequences:
    new_sequences.append(ele.replace('-',''))

out = open('new_protein_1000.fasta', 'w')
for i in range(len(names)):
    if i in to_be_erased:
        continue
    else:
        out.write('>' + names[i] + '\n')
        for j in range((len(new_sequences[i])//80)+1):
            out.write(new_sequences[i][(j-1)*80 : j*80] + '\n')
        if j*80 < len(new_sequences[i]):
            out.write(new_sequences[i][j*80:] + '\n')