prot_dict = {'G':0 , 'A':1 , 'L':2 , 'M':3 , 'F':4 , 'W':5 , 'K':6 , 'Q':7 , 'E':8 , 'S':9 , 'P':10 , 'V':11 , 'I':12 , 'C':13 , 'Y':14 , 'H':15 , 'R':16 , 'N':17 , 'D':18 , 'T':19 , '-':20, 'X':21 }

def fastareader(filename):
    sequences = []
    my_file = open(filename,'r')
    seq = ""
    start = True
    names = []

    for line in my_file:
        line = line.replace('\n', '')
        if len(line)==0:{}    #I encountered some fasta files generated where there is a blank line between two sequences 
        if(line[0]=='>'):
            if(start):
                start = False
                names.append(line[1:])
            else:
                sequences.append(seq)
                seq = ""
                names.append(line[1:])
        else:
            seq = seq + line
    sequences.append(seq)
    return sequences, names


blosum62_mat = [[]]*22
for i in range(len(blosum62_mat)):
    blosum62_mat[i] = [[]]*22

blosum62_mat[4][5] = 1; blosum62_mat[2][16] = -2; blosum62_mat[9][10] = -1; blosum62_mat[11][19] = 0
blosum62_mat[7][7] = 5; blosum62_mat[1][17] = -2; blosum62_mat[5][16] = -3
blosum62_mat[1][7] = -1; blosum62_mat[9][18] = 0; blosum62_mat[15][15] = 8; blosum62_mat[9][15] = -1
blosum62_mat[15][18] = -1; blosum62_mat[2][17] = -3; blosum62_mat[1][5] = -3; blosum62_mat[3][14] = -1
blosum62_mat[0][16] = -2; blosum62_mat[12][14] = -1; blosum62_mat[8][14] = -2
blosum62_mat[1][14] = -2; blosum62_mat[11][18] = -3; blosum62_mat[14][14] = 7
blosum62_mat[0][17] = 0; blosum62_mat[8][13] = -4; blosum62_mat[8][14] = -1
blosum62_mat[1][11] = 0; blosum62_mat[13][13] = 9; blosum62_mat[3][16] = -1; blosum62_mat[8][11] = -2
blosum62_mat[17][19] = 0; blosum62_mat[10][10] = 7; blosum62_mat[11][12] = 3; blosum62_mat[9][11] = -2
blosum62_mat[3][11] = 1; blosum62_mat[4][19] = -2; blosum62_mat[7][11] = -2
blosum62_mat[6][6] = 5; blosum62_mat[10][18] = -1; blosum62_mat[12][15] = -3; blosum62_mat[12][18] = -3
blosum62_mat[16][19] = -1; blosum62_mat[2][10] = -3; blosum62_mat[0][6] = -2; blosum62_mat[3][17] = -2
blosum62_mat[10][15] = -2; blosum62_mat[4][7] = -3; blosum62_mat[2][21] = -1
blosum62_mat[3][19] = -1; blosum62_mat[15][21] = -1; blosum62_mat[16][18] = -2
blosum62_mat[18][21] = -1; blosum62_mat[1][4] = -2
blosum62_mat[4][8] = -3; blosum62_mat[17][18] = 1
blosum62_mat[21][21] = -1; blosum62_mat[4][12] = 0; blosum62_mat[19][21] = 0
blosum62_mat[3][4] = 0
blosum62_mat[9][9] = 4; blosum62_mat[2][7] = -2; blosum62_mat[5][8] = -3; blosum62_mat[7][16] = 1
blosum62_mat[17][17] = 6; blosum62_mat[3][5] = -1; blosum62_mat[7][13] = -3; blosum62_mat[5][12] = -3
blosum62_mat[9][13] = -1; blosum62_mat[1][2] = -1; blosum62_mat[0][9] = 0; blosum62_mat[2][8] = -3
blosum62_mat[5][7] = -2; blosum62_mat[0][15] = -2; blosum62_mat[6][9] = 0; blosum62_mat[7][17] = 0
blosum62_mat[16][17] = 0; blosum62_mat[13][15] = -3; blosum62_mat[14][17] = -2; blosum62_mat[0][7] = -2
blosum62_mat[4][14] = 3; blosum62_mat[1][13] = 0; blosum62_mat[2][11] = 1; blosum62_mat[0][8] = -2
blosum62_mat[0][1] = 0; blosum62_mat[6][16] = 2; blosum62_mat[8][18] = 2; blosum62_mat[14][16] = -2
blosum62_mat[3][7] = 0; blosum62_mat[12][19] = -1; blosum62_mat[13][18] = -3; blosum62_mat[4][11] = -1
blosum62_mat[1][19] = 0; blosum62_mat[10][19] = -1; blosum62_mat[8][19] = -1
blosum62_mat[11][17] = -3; blosum62_mat[0][10] = -2; blosum62_mat[1][3] = -1; blosum62_mat[6][15] = -1
blosum62_mat[11][16] = -3; blosum62_mat[10][13] = -3; blosum62_mat[3][8] = -2; blosum62_mat[2][6] = -2
blosum62_mat[11][11] = 4; blosum62_mat[3][12] = 1; blosum62_mat[7][19] = -1; blosum62_mat[0][12] = -4
blosum62_mat[6][10] = -1; blosum62_mat[3][3] = 5; blosum62_mat[6][18] = -1; blosum62_mat[12][13] = -1
blosum62_mat[4][16] = -3; blosum62_mat[6][21] = -1; blosum62_mat[7][18] = 0
blosum62_mat[0][21] = -1; blosum62_mat[13][21] = -2
blosum62_mat[4][4] = 6; blosum62_mat[5][21] = -2
blosum62_mat[1][18] = -2; blosum62_mat[2][9] = -2; blosum62_mat[9][21] = 0
blosum62_mat[4][17] = -3; blosum62_mat[9][16] = -1; blosum62_mat[5][18] = -4; blosum62_mat[11][14] = -1
blosum62_mat[2][5] = -2; blosum62_mat[15][16] = 0; blosum62_mat[5][15] = -2; blosum62_mat[15][17] = 1
blosum62_mat[5][19] = -2; blosum62_mat[19][19] = 5; blosum62_mat[4][9] = -2; blosum62_mat[5][10] = -4
blosum62_mat[2][18] = -4; blosum62_mat[2][15] = -3; blosum62_mat[9][17] = 1
blosum62_mat[2][2] = 4; blosum62_mat[6][14] = -2; blosum62_mat[7][8] = 2
blosum62_mat[0][14] = -3; blosum62_mat[13][14] = -2; blosum62_mat[0][18] = -1
blosum62_mat[1][8] = -1; blosum62_mat[5][14] = 2; blosum62_mat[8][8] = 5
blosum62_mat[9][14] = -2; blosum62_mat[13][17] = -3; blosum62_mat[11][13] = -1; blosum62_mat[15][19] = -2
blosum62_mat[10][16] = -2; blosum62_mat[0][11] = -3; blosum62_mat[2][19] = -1; blosum62_mat[6][11] = -2
blosum62_mat[6][7] = 1; blosum62_mat[1][16] = -1; blosum62_mat[12][16] = -3; blosum62_mat[18][19] = -1
blosum62_mat[4][10] = -4; blosum62_mat[12][17] = -3; blosum62_mat[6][12] = -3; blosum62_mat[3][18] = -3
blosum62_mat[5][11] = -3; blosum62_mat[5][5] = 11; blosum62_mat[3][15] = -2; blosum62_mat[10][17] = -2
blosum62_mat[1][6] = -1; blosum62_mat[2][3] = 2; blosum62_mat[6][8] = 1
blosum62_mat[17][21] = -1; blosum62_mat[4][21] = -1
blosum62_mat[6][13] = -3
blosum62_mat[4][13] = -2; blosum62_mat[0][4] = -3
blosum62_mat[11][21] = -1; blosum62_mat[4][6] = -3
blosum62_mat[16][21] = -1; blosum62_mat[18][18] = 6; blosum62_mat[0][5] = -2
blosum62_mat[7][9] = 0; blosum62_mat[5][13] = -2; blosum62_mat[5][6] = -3; blosum62_mat[7][15] = 0
blosum62_mat[2][13] = -1; blosum62_mat[5][17] = -4; blosum62_mat[1][9] = 1; blosum62_mat[0][2] = -4
blosum62_mat[5][9] = -3; blosum62_mat[8][9] = 0; blosum62_mat[8][15] = 0; blosum62_mat[9][12] = -2
blosum62_mat[1][15] = -2; blosum62_mat[3][9] = -1; blosum62_mat[2][14] = -1; blosum62_mat[14][15] = 2
blosum62_mat[14][18] = -3; blosum62_mat[8][16] = 0; blosum62_mat[10][21] = -2; blosum62_mat[0][0] = 6
blosum62_mat[0][13] = -3; blosum62_mat[8][17] = 0; blosum62_mat[14][19] = -2; blosum62_mat[10][14] = -3
blosum62_mat[6][19] = -1; blosum62_mat[1][1] = 4; blosum62_mat[7][10] = -1; blosum62_mat[13][19] = -1
blosum62_mat[11][15] = -3; blosum62_mat[0][19] = -2; blosum62_mat[7][12] = -3
blosum62_mat[13][16] = -3; blosum62_mat[10][11] = -2; blosum62_mat[8][10] = -1; blosum62_mat[3][13] = -1
blosum62_mat[6][17] = 0; blosum62_mat[12][12] = 4; blosum62_mat[1][10] = -1; blosum62_mat[0][3] = -3
blosum62_mat[9][19] = 1; blosum62_mat[8][12] = -3; blosum62_mat[3][10] = -2; blosum62_mat[3][6] = -1
blosum62_mat[1][12] = -1; blosum62_mat[10][12] = -3; blosum62_mat[16][16] = 5; blosum62_mat[3][21] = -1
blosum62_mat[2][12] = 2; blosum62_mat[12][21] = -1; blosum62_mat[8][21] = -1
blosum62_mat[1][21] = 0
blosum62_mat[4][18] = -3; blosum62_mat[14][21] = -1; blosum62_mat[4][15] = -1
blosum62_mat[2][4] = 0; blosum62_mat[7][21] = -1

for i in range(21):
    blosum62_mat[i][20] = 0

blosum62_mat[20][21] = 0

# for i in range(len(blosum62_mat)):
#     print(blosum62_mat[i])