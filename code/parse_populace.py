def get_matrices(nonharmful_file,harmful_file,allele_min,db_info):
    three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','TER': '-'}
    my_file = open(nonharmful_file,'r')
    non_harmful_all = []
    non_harmful_count = 0
    for line in my_file:
        if(len(line)<0):
            continue
        if(line[0] != '2'):
            continue
        temp = line.split(';')
        conseq = temp[8]
        if conseq[0] == 'c':
            if '-' in conseq or '+' in conseq or '*' in conseq:
                continue
            else:
                print('there is an error')
        if conseq[0] == 'p':
            prev_aa = three_to_one[conseq[2:5].upper()]
            pos = conseq[5:-3]
            if conseq[-1:] == '?':
                continue
            if 'LeuextTer' in conseq:
                continue
            if '513_Glu514' in conseq:
                continue
            if '3_Thr4insAlaArgAlaSerMetGly' in conseq:
                continue
            if 'fs' in conseq:
                continue
            if (conseq[-3:] == 'del'):
                new_aa = '-'
            else:
                new_aa = three_to_one[conseq[-3:].upper()]
            non_harmful_all.append([prev_aa,pos,new_aa,temp[13]])
            non_harmful_count = non_harmful_count + int(temp[13])
    
    non_harmful = []
    for ele in non_harmful_all:
        if int(ele[3]) > allele_min:
            non_harmful.append(ele)
    if(db_info==True):
        print('number of total unique nonharmful mutations: '+ str(len(non_harmful_all)))
        print('number of unique nonharmful mutations with occurance threshold: '+ str(len(non_harmful)))

    missense = 0
    my_file2 = open(harmful_file)
    harmful0 = []
    harmful1 = []
    harmful2 = []
    for line in my_file2:
        if(len(line)<0):
            continue
        conseq = line.split('\t')[0]
        if '-' in conseq or '+' in conseq or '*' in conseq: 
            continue
        t = -1
        if 'diabetes' in line:
            t = 0
            if 'Likely pathogenic' in line:
                t = 1
            elif 'Pathogenic' in line:
                t = 2
        if t == -1:
            continue
        temp2 = conseq.split('(')[2][2:-1]
        # if 'Uncertain significance' in line or 'Benign' in line or 'Likely benign' in line or 'Conflicting interpretations of pathogenicity' in line:
        #     continue
        if '=' in conseq:
            continue
        if 'fs' in conseq:
            prev_aa = three_to_one[temp2[:3].upper()]
            pos = temp2[3:-2]
            if temp2[-2:] != 'fs':
                print('there is a problem')
            else:
                if t==0:
                    harmful0.append([prev_aa,pos,'fs'])
                elif t==1:
                    harmful1.append([prev_aa,pos,'fs'])
                else:
                    harmful2.append([prev_aa,pos,'fs'])
        else:
            prev_aa = three_to_one[temp2[:3].upper()]
            pos = temp2[3:-3]
            next_aa = three_to_one[temp2[-3:].upper()]
            if t == 0:
                harmful0.append([prev_aa,pos,next_aa])
            elif t==1:
                harmful1.append([prev_aa,pos,next_aa])
            else:
                harmful2.append([prev_aa,pos,next_aa])
            missense = missense + 1
    if(db_info==True):
        print('number of non pathogenic mutations with symptoms: ' + str(len(harmful0)))
        print('number of likely pathogenic mutations with symptoms: ' + str(len(harmful1)))
        print('number of pathogenic mutations with symptoms: ' + str(len(harmful2)))
    harmful = [harmful0,harmful1,harmful2]
    return [harmful,non_harmful]
        