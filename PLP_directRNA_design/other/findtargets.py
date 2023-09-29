def findtargets (mrna,refpath,ie,outfiles,plp_length=30,gc_min=50,gc_max=65):
    targets = pd.DataFrame(columns=['Gene', 'Position', 'Sequence'])
    end = len(mrna)-(plp_length-1)
    #print (end)
    for i in range (0, end):
        #print (mrna[i:i+30])
    #The next line checks if position 16 (remember python is 0-indexed) is a C or G
        if mrna.seq[i+round(plp_length/2)] == 'C' or mrna.seq[i+round(plp_length/2)] == 'G' :
            #The next line filters out any probe with GC content <= 50 and >=65
            if GC(mrna.seq[i:i+plp_length]) > gc_min:
                if GC(mrna.seq[i:i+plp_length]) < gc_max:
                    if mrna.seq[i:i+plp_length].count("AAA")==0 and mrna.seq[i:i+plp_length].count("TTT")==0 and mrna.seq[i:i+plp_length].count("GGG")==0 and mrna.seq[i:i+plp_length].count("CCC")==0:
                    #Here I create a dataframe with all the suitable targets, where column 1 is the start position and column 2 is the actual sequence.
                        #print (GC(mrna.seq[i:i+30]))
                        targets = targets.append({'Gene': mrna.id, 'Position': i, 'Sequence':mrna.seq[i:i+plp_length]}, ignore_index=True)  
                        pato=refpath+ '/target_regions_'+mrna.id+'_'+str(ie)+'.csv'
                        outfiles.append(pato)
                        targets.to_csv(pato)
    return [targets,outfiles]   