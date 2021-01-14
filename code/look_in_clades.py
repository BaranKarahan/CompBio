from look_in_nwk import create_clades
from ens210_project_lib import prot_dict, fastareader
import copy

def compare_aa_group(aa1,aa2):
    if aa1 == 'R' or aa1 == 'H' or aa1 == 'K' or aa1 == 'D' or aa1 == 'E': match = 0
    elif aa1 == 'S' or aa1 == 'T' or aa1 == 'N' or aa1 == 'Q': match = 1
    elif aa1 == 'C' or aa1 == 'U' or aa1 == 'G' or aa1 == 'P': match = 2
    else: match = 3

    if aa2 == 'R' or aa2 == 'H' or aa2 == 'K' or aa2 == 'D' or aa2 == 'E': match2 = 0
    elif aa2 == 'S' or aa2 == 'T' or aa2 == 'N' or aa2 == 'Q': match2 = 1
    elif aa2 == 'C' or aa2 == 'U' or aa2 == 'G' or aa2 == 'P': match2 = 2
    else: match2 = 3
    if match == match2:
        return 1
    else:
        return 0


def iterate_clades(filename):
    sequences, names = fastareader(filename)
    leaves = create_clades("tree.nwk")
    ind_set = []
    seq_set = []
    aa_occurance = []
    for i in range(10):
        for leaf in leaves[i]:
            this_name = leaf.replace('_',' ').replace("'",'')
            ind_set.append(names.index(this_name))
            if 'Homo_sapiens' in leaf:
                hom_ind = names.index(leaf.replace('_',' ').replace("'",''))
        hom_ind2 = ind_set.index(hom_ind)
        hom_ind3 = ind_set[hom_ind2]
        for ind in ind_set:
            seq_set.append(sequences[ind])
        prot_mat = []
        for j in range(len(seq_set[0])):
            prot_vec = [0]*22
            total = 0
            for seq in seq_set:
                this_prot = prot_dict[seq[j]]
                prot_vec[this_prot] = prot_vec[this_prot]+1
                total = total + 1
            for e in range(len(prot_vec)):
                prot_vec[e] = prot_vec[e] / total
            prot_mat.append(copy.deepcopy(prot_vec))
        hum_seq = sequences[hom_ind3]
        hum_var = []
        for j in range(len(hum_seq)):
            this_prot = prot_dict[hum_seq[j]]
            if this_prot < 20:
                hum_var.append(prot_mat[j][this_prot])
        
        aa_occurance.append(copy.deepcopy(hum_var))
    return ([aa_occurance,hum_seq])




    


    