def create_clades(filename):
    #filename = "tree.nwk"
    my_file = open(filename,'r')
    newick = ""
    for line in my_file:
        newick=line
    ind= newick.find('Homo_sapiens')
    lind = newick[:ind].rfind('(')-1
    rind = newick[ind:].find(')') + ind+1


    l = lind
    r= rind
    rend = False
    lend = False
    len_nwk = len(newick)
    clades = []
    clades.append(newick[l:r])
    #homo sapiens are in 10 clades
    for i in range(9):
        l_dist = 1
        r_dist = 1
        while l_dist!=0:
            l = l - 1
            if newick[l] == '(':
                l_dist = l_dist - 1
            elif newick[l] == ')':
                l_dist = l_dist + 1
            if l == 0:
                lend = True
                break
        
        while r_dist!=0:
            r = r + 1
            if newick[r] == '(':
                r_dist = r_dist + 1
            elif newick[r] == ')':
                r_dist = r_dist - 1
            if r == len_nwk:
                rend = True
                break
        if lend == True:
            if rend == True:
                break
            else: print('OOPS')
        else:
            if(rend==True):
                print('OOPS')
        clades.append(newick[l:r+1])
    leaves = []
    this_leaf = ""
    for i in range(len(clades)):
        temp = []
        for letter in clades[i]:
            if letter == '(' or letter == ')' or letter ==',':
                if len(this_leaf)>0:
                    temp.append(this_leaf)
                    this_leaf = ""
            else:
                this_leaf = this_leaf+letter
        leaves.append(temp)
    return(leaves)