filename = "clinvar_result.txt"
my_file = open(filename,'r')
patient1 = [] #most sure
patient2 = [] #less sure
patient3 = [] #least sure
patient_info = []
pat_no = -1
for line in my_file:
    pat_no = pat_no+1
    dont_know = True
    line = line.replace('\n', '')
    
    tmp = line.find('\t')
    line_ = line[:tmp]
    tmp = line_.find(':c.')
    line_ = line_[tmp+3:]
    if '(p.' in line_:
         line_ = line_[:line_.find('(p.')-1]
    if '*' in line_:
        line_ = line_[1:]
    if '-' in line_:
        line_ = line_[1:]
    mtype = -1
    if '>' in line_:
        mtype = 0
        tmp = line_.find('>')
        temp1 = line_[:tmp-1]
        temp2 = line_[tmp-1:]
    elif 'dup' in line_:
        mtype = 1
        tmp = line_.find('>')
        temp1 = line_[:tmp-2]
        temp2 = ""
    elif 'del' in line_:
        mtype = 2
        tmp = line_.find('del')
        temp1 = line_[:tmp-1]  
        temp2 = "" 
    elif 'ins' in line_: #no example of this
        mtype = 3
        tmp = line_.find('ins')
        temp1 = line_[:tmp-1]
        temp2 = ""
    else:
        continue #I could not understand this HGVS  for example: 3008_*3011TTGT[1] (there is no >,del,dup,ins)
    if mtype > -1:
        patient_info.append([mtype,temp1,temp2])

    patient3.append(pat_no)
    if "Hypogonadism, diabetes mellitus, alopecia, mental retardation and electrocardiographic abnormalities" in line:
        patient2.append(pat_no)

    if  ("diabetes" in line) and ("pathogenic" in line)  and ("conficting not in line"):
            patient1.append(pat_no)















ORF = [[]]*43
ORF[6]= ['+',	353, 1915]
ORF[42]= ['-',	514, 17]
ORF[28]= ['-',	456,	223]
ORF[5]= ['+',	4552,	4779]
ORF[32]= ['-',	3461,	3258]
ORF[38]= ['-',	4021,	3824]
ORF[37]= ['-',	5584,	5396]
ORF[31]= ['-',	3965,	3795]
ORF[35]= ['-',	2003,	1845]
ORF[27]= ['-',	1194,	1057]
ORF[20]= ['+',	5691,	5828]
ORF[23]= ['-',	4002,	3868]
ORF[34]= ['-',	2423,	2292]
ORF[14]= ['+',	2235,	2366]
ORF[22]= ['-',	4257,	4129]
ORF[30]= ['-',	4664,	4539]
ORF[9]= ['+',	4271,	4393]
ORF[11]= ['+',	5627,	5749]
ORF[36]= ['-',	695,	576]
ORF[24]= ['-',	3417,	3298]
ORF[29]= ['-',	5843,	5724]
ORF[12]= ['+',	534,	644]
ORF[17]= ['+',	3957,	4061]
ORF[8]= ['+',	3920,	4024]
ORF[25]= ['-',	3249,	3145]
ORF[2]= ['+',3385,	3486]
ORF[33]= ['-',3032,	2931]
ORF[15]= ['+',	2646,	2741]
ORF[39]= ['-',	2248,	2156]
ORF[26]= ['-',	2667,	2482]
ORF[7]= ['+',	3020,	3112]
ORF[19]= ['+',	5010,	5102]
ORF[13]= ['+',	1230,	1319]
ORF[21]= ['-',	4728,	4639]
ORF[41]= ['-',	1306,	1217]
ORF[16]= ['+',	3861,	3947]
ORF[10]= ['+',	4682,	4765]
ORF[1]= ['+',	622,	705]
ORF[4]= ['+',	3964,	4044]
ORF[3]= ['+',	3880,	3960]
ORF[18]= ['+',	4131,	4211]
ORF[40]= ['-',	1531,	1454]


counter = 0
mem = [0] * len(ORF)
pat_set = patient3 #patient1/patient2/patient3
for i in range(len(pat_set)):
    m_info = patient_info[i]
    test = False
    change_loc = m_info[1]
    if ('+' in change_loc):
        temp = change_loc.find('+')
        change_loc = int(change_loc[:temp]) + int(change_loc[temp+1:])

    elif ('-' in change_loc):
        temp = change_loc.find('-')
        change_loc = int(change_loc[:temp]) - int(change_loc[temp+1:])

    change_loc = int(change_loc)
    for j in range(1,len(ORF)):
        if ORF[j][0] == '+':
            pl = True
        else:
            pl = False

        if pl==True:
            mi = ORF[j][1]
            ma = ORF[j][2]
        else:
            ma = ORF[j][1]
            mi = ORF[j][2]

        if (test==False) and (change_loc < ma) and (change_loc > mi):
            test = True
            mem[j] = mem[j]+1

    if(test == False):
        #print(m_info)
        counter = counter+1
print(counter)
print(len(pat_set))
print(mem) #shows how many mutations are on each ORF




#NM_025000.4 is the referenced genome for all of the patient data it spans from . i wrote it as gen1
#print(len(gen1P),len(gen1N))
#print(patient1)