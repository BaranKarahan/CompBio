def fastareader(filename):
    seqDict = {}
    # Code here
    my_file = open(filename,'r')
    head = ""
    seq = ""
    start = True
    debugg = 0
    for line in my_file:
        line = line.replace('\n', '')
        if len(line)==0:    #I encountered some fasta files generated where there is a blank line between two sequences 
            continue
        if(line[0]=='>'):
            if(start):
                start = False
                head = line[1:]
            else:
                seqDict[head] = seq
                head = line[1:]
                seq = ""
        else:
            seq = seq + line
    seqDict[head] = seq
    return seqDict


seq_dict = fastareader("Dcaf17_refseq_transcript.fasta")
print(seq_dict['NM_025000.4 Homo sapiens DDB1 and CUL4 associated factor 17 (DCAF17), transcript variant 1, mRNA'])