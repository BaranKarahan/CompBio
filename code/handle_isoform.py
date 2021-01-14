import copy
from ens210_project_lib import prot_dict, fastareader

sequences, names = fastareader("seqdump.fas")

count = 0
count2 = 0
isoforms = []
non_isoformed = []
for i in range(len(names)):
    if 'isoform'  in names[i]:
        count = count + 1
        isoforms.append([i,names[i]])
    else:
        count2 = count2 + 1
        non_isoformed.append(i) 


isoform_list = []
species = []
for i in range(len(isoforms)):
#for i in range(1):
    this_sp = isoforms[i][1].split('isoform')[1]
    temp = this_sp.split('[')
    iso_no = temp[0]
    this_species = temp[1][:-1]
    isoform_list.append([isoforms[i][0],this_species,iso_no])
    species.append(this_species)


iso_occ = ([ (i,species.count(i)) for i in set(species) ])
spe_num = list(set(species))
print(len(spe_num))
preferred = {}
for iso in isoform_list:
    spe = iso[1]
    seq_len = len(sequences[iso[0]])
    if spe in preferred:
        if preferred[spe][1] < seq_len:
            preferred[spe] = [iso[0],seq_len]
    else:
        preferred[spe] = [iso[0],seq_len]
if len(spe_num) != len(preferred):
    print('one preferred isoform for each of the species did not found')


sequences_left = count2 + len(preferred)
print('number of sequences left: ' + str(sequences_left))





# newFasta = open('new_seqdump.fas', 'w')
# for noni in non_isoformed:
#     ind = int(noni)
#     newFasta.write('>'+names[ind]+'\n')
#     newFasta.write(sequences[ind]+'\n')
#     newFasta.flush()

# for pref in preferred:
#     ind = int(preferred[pref][0])
#     newFasta.write('>'+names[ind]+'\n')
#     newFasta.write(sequences[ind]+'\n')
#     newFasta.flush()

