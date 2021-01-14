from look_in_clades import iterate_clades
from parse_populace import get_matrices
from ens210_project_lib import prot_dict

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
mode = 1 ####
db_info = False
ch = 4
harmful_file = 'clinvar_result.txt'
nonharmful_file = 'human_excel.csv'
[aa_occurance,hum_seq] = iterate_clades('sequences.fas')
[harmful_mat,nonharmful] = get_matrices(nonharmful_file,harmful_file,10,db_info)
harmful = []
for i in range(mode,3):
    harmful.extend(harmful_mat[i])
print('harmful set size: ' + str(len(harmful)) + ' nonharmful set size: '+str(len(nonharmful)))
results = []

if ch==0:
    se=range(9,10)
elif ch==1:
    se = range(10)
elif ch==2:
    se = range(10)
elif ch==3:
    se = range(7,9)
elif ch==4:
    se = range(8,10)
t_ = range(11)
if ch == 2:
    t_ = range(7,9)



for t in t_:
    threshold = 0.5 + (t*0.05)
    for i in se:
        res_mat = [[0,0],[0,0]]
        occ = aa_occurance[i]
        for m in harmful:
            if m[2]!= 'fs':
                harm_res = float(occ[int(m[1])])
                if harm_res > threshold:
                    if m[2]== '-' or compare_aa_group(m[0],m[2]) == 0:
                        res_mat[0][0] = res_mat[0][0]+1
                    else:
                        res_mat[0][1] =  res_mat[0][1]+1
                    
                else:
                    res_mat[0][1] = res_mat[0][1]+1
            else: res_mat[0][0] = res_mat[0][0]+1
        for n in nonharmful:
            non_res = float(occ[int(n[1])])
            if non_res > threshold:
                if compare_aa_group(n[0],n[2]) == 1:
                    res_mat[1][1] = res_mat[1][1]+1
                else:
                    res_mat[1][0] = res_mat[1][0]+1 
            else:
                res_mat[1][1] = res_mat[1][1]+1
        results.append([i,threshold,res_mat])
accuracies = []
for res in results:
    this_mat = res[2]
    tot = this_mat[0][0] + this_mat[0][1] + this_mat[1][0] + this_mat[1][1]
    acc = (this_mat[0][0] + this_mat[1][1]) / tot
    sens = (this_mat[0][0] ) / (this_mat[0][0] + this_mat[0][1]) #
    spec = (this_mat[1][1] ) / (this_mat[1][1] + this_mat[1][0])

    th = res[1]
    if th>0.85 and th<0.9:
        th = 0.85
    if th>0.6 and th<0.65:
        th = 0.6
    if th>0.7 and th<0.75:
        th = 0.7
    #str(round(answer, 2))

    accuracies.append([th,round(acc,2),round(sens,2),round(spec,2),res[0]])
for acc in accuracies:
    #print(acc)
    if ch > 0:
        if acc[0]*100 % 10 >0:
            print('threshold = ' + str(acc[0]) + ', clade height = ' + str(acc[4]) +' | Accuracy = '+ str(acc[1])+', Sensitivity = ' + str(acc[2])+', Specificity = ' + str(acc[3]))
        else:
            print('threshold = ' + str(acc[0]) + ',  clade height = ' + str(acc[4]) + ' | Accuracy = '+ str(acc[1])+', Sensitivity = ' + str(acc[2])+', Specificity = ' + str(acc[3]))
    else:    
        if acc[0]*100 % 10 >0:
            print('threshold = ' + str(acc[0]) +' | Accuracy = '+ str(acc[1])+', Sensitivity = ' + str(acc[2])+', Specificity = ' + str(acc[3]))
        else:
            print('threshold = ' + str(acc[0]) + '  | Accuracy = '+ str(acc[1])+', Sensitivity = ' + str(acc[2])+', Specificity = ' + str(acc[3]))




# print(len(aa_occurance[5]))

# print(aa_occurance[5])
# for i in range(len(aa_occurance[5])):
#     if aa_occurance[5][i] < 0.75:
#         print(i)


