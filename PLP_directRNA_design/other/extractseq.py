def extractseq(goi,ref):
    in_f =open (ref, 'r+')
    seqs=[]
    su=0
    for line in in_f:
    #    print(line)
        if su==1:
            if line[0]=='>':
                su=0
                if lout[0]!='>':
    #                print(lout)
                    seqs.append(lout)
            else:
                lout=lout+line[:-1]       
        if line in goi:
            seqs.append(line)
            su=1
            lout=''
    return seqs