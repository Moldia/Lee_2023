def extract_seqs(genes,ref,column='Gene'):
    in_f =open (ref, 'r+')
    headers=list()
    for line in in_f:
    #    print(line)
        if line[0] == '>':
            headers.append(line)
    #        print(line)
    dictionary=dict()
    header=[]
    seq=[]
    genesexp=genes[column].unique()
    listo=[]
    lista=[]
    for ici in genesexp:
        icis='('+ici+')'
        matching = [s for s in headers if icis in s]
        lista.append(len(matching))
        listo.append(matching)
    
    return(genesexp,listo,lista)