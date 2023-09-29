from weakref import ref


def extractseq(goi,ref):

    '''
    The extractseq function takes two arguments: goi, a list of specific genes of interest, 
    and ref, the path to a reference file containing biological sequences. 
    The function reads through the reference file line by line and extracts the sequences 
    that correspond to the genes specified in the goi list. The extracted sequences are stored 
    in a list called seqs, which is then returned. 
    The function is designed to work with files where each biological sequence starts 
    with a '>' character followed by metadata, and the actual sequence spans one or more 
    lines until the next '>' character appears. Note that the function expects the 
    reference file to be formatted in such a way that each gene and its corresponding 
    sequence are adjacent lines in the file.
    '''
    
    import pandas as pd
    from pandas import DataFrame
    import Bio
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    
    # Open the reference file for reading and writing ('r+')
    in_f = open(ref, 'r+')
    
    # Initialize an empty list 'seqs' to store extracted sequences
    seqs = []
    
    # Initialize a variable 'su' to keep track of sequence extraction
    su = 0
    
    # Loop through each line of the reference file
    for line in in_f:
        # Check if su is 1, which means the sequence should be captured
        if su == 1:
            # If the line starts with '>', it indicates the end of the sequence
            if line[0] == '>':
                su = 0
                # Append the full sequence to 'seqs' if it does not start with '>'
                if lout[0] != '>':
                    seqs.append(lout)
            # Otherwise, concatenate the current line to the ongoing sequence
            else:
                lout = lout + line[:-1]
                
        # If the current line is in 'goi' (genes of interest), set su to 1 and reset lout
        if line in goi:
            seqs.append(line)
            su = 1
            lout = ''
            
    # Return the list of extracted sequences
    return seqs

def extract_seqs(genes, ref, column='Gene'):


    '''
    The function takes in three arguments:

    genes: A DataFrame containing gene information. The DataFrame must contain a column specified by the column parameter (default to 'Gene').
    ref: The path to a file that contains sequence information in FASTA format.
    column: The column name to look for in the genes DataFrame when extracting sequences (default is 'Gene').

    The function extract all the mRNA sequences corresponding to a specific gene, by doing this:

    Open the reference file specified by ref and read its content.
    Initialize several lists (headers, header, seq, listo, lista) and a dictionary (dictionary)
    as placeholders for various data.
    Loop through each line in the reference file to find the headers (lines starting with >), 
    appending each to the headers list.
    Fetch the unique gene identifiers (genesexp) based on the specified column in the genes DataFrame.
    Loop through each unique gene identifier to find matching headers in the headers list.
    Populate the lists lista and listo with the number of matching headers and the headers 
    themselves, respectively.
    Return the three lists: genesexp (unique gene identifiers), listo (matching headers), 
    and lista (counts of matching headers).


    '''
    # Import required libraries
    import pandas as pd
    from Bio import SeqIO, SeqUtils, Seq, SeqRecord, AlignIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    import random
    import numpy as np

    # Open the reference file for reading
    in_f = open(ref, 'r+')

    # Initialize a list to hold the headers from the reference file
    headers = list()

    # Loop through each line in the reference file
    for line in in_f:
        # Check if the line starts with '>', indicating it's a header
        if line[0] == '>':
            headers.append(line)

    # Initialize variables for later use
    dictionary = dict()
    header = []
    seq = []
    # Extract unique gene names from the DataFrame based on the specified column
    genesexp = genes[column].unique()
    
    listo = []  # To hold the matching headers
    lista = []  # To hold the count of matching headers

    # Loop through each unique gene identifier
    for ici in genesexp:
        icis = '(' + ici + ')'
        # Find headers that contain the unique gene identifier
        matching = [s for s in headers if icis in s]
        # Add the count of matching headers to lista
        lista.append(len(matching))
        # Add the matching headers themselves to listo
        listo.append(matching)

    # Return the unique gene identifiers, matching headers, and counts
    return genesexp, listo, lista

def findtargets (mrna,refpath,ie,outfiles,plp_length=30,gc_min=50,gc_max=65,ligationsite_GC=True):
    '''
    The function extracts, for a given mRNA, all the possible k-mers.
    k (k-mer length) is specified by plp_length (int, defaults to 30).
    The k-mers are then filtered by upper and lower GC content (defaults = 50-65),
    and screened for G or C at the ligation site (beneficial but not crucial, optional but default =True)
    '''
    import pandas as pd 
    from pandas import DataFrame
    import Bio
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    targets = pd.DataFrame(columns=['Gene', 'Position', 'Sequence'])
    end = len(mrna)-(plp_length-1)
    #print (end)
    for i in range(0, end):
        # Conditional check for the central base ('C' or 'G') based on ligationsite_GC flag
        if ligationsite_GC:
            central_condition = mrna.seq[i+round(plp_length/2)] == 'C' or mrna.seq[i+round(plp_length/2)] == 'G'
        else:
            central_condition = True
        
        if central_condition:
            if GC(mrna.seq[i:i+plp_length]) > gc_min and GC(mrna.seq[i:i+plp_length]) < gc_max:
                if mrna.seq[i:i+plp_length].count("AAA")==0 and mrna.seq[i:i+plp_length].count("TTT")==0 and mrna.seq[i:i+plp_length].count("GGG")==0 and mrna.seq[i:i+plp_length].count("CCC")==0:
                    targets = targets.append({'Gene': mrna.id, 'Position': i, 'Sequence':mrna.seq[i:i+plp_length]}, ignore_index=True)
                    pato = refpath + '/target_regions_' + mrna.id + '_' + str(ie) + '.csv'
                    outfiles.append(pato)
                    targets.to_csv(pato)
    return [targets,outfiles] 

def plot_alignment(refpath,alignment,common):
    """
    Plots common regions through sequence variants.

    Parameters:
    refpath (str): The directory path where the plot will be saved.
    alignment (Bio.Align.MultipleSeqAlignment): The multiple sequence alignment object.
    common (list): List or array specifying common regions in alignment.
    """
    import pandas as pd 
    from pandas import DataFrame
    import Bio
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(10,3))
    plt.plot(range(0,alignment.get_alignment_length()),common,c='grey')
    plt.hlines(len(alignment)+0.1,linestyles='--',xmin=0,xmax=alignment.get_alignment_length(),colors='r')
    plt.xlim([0,alignment.get_alignment_length()])
    plt.savefig(refpath+'/common_regions_though_variants.png')
    plt.close(fig)
    
def extract_and_align(genes, ref, path, pathclustal, plp_length=30, gc_min=50, gc_max=65, ligationsite_GC=True):
    '''
    The purpose of this function is to extract and align multiple sequences 
    for a given list of genes from a reference.
    The function requires the user to have ClustalW2 installed and specify its path
    If multiple sequences are available for a gene, the function uses ClustalW to perform 
    sequence alignment. It also calculates the most common characters at each position 
    and uses this information for further analysis or visualization.

    Some arguments as passed to the function findtargets, so please refer to that function for specs.

    notfound: List of genes that were not found.
    '''
    import pandas as pd 
    import os
    from pandas import DataFrame
    import Bio
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    import matplotlib.pyplot as plt
    isExist = os.path.exists(path+'/gene_info')
    if not isExist:
        os.makedirs(path+'/gene_info') 
    outfiles=[]
    genesexp,listoref,lista=extract_seqs(genes,ref)
    listaref=lista
    lista=[]
    listo=[]
    for el in listoref:
        ls=[]
        for ele in el:
            check=ele.find('PREDICTED')==-1
            if check==True:
                ls.append(ele)
        lista.append(len(ls))
        listo.append(ls)
    notfound=[]
    for holi in range(0,len(genesexp)):
        ie=0
        refpath=path+"/gene_info/"+genesexp[holi]
        if lista[holi]<2:
            if lista[holi]==0:
                print("Gene "+genesexp[holi] +" was not found")
                notfound.append(genesexp[holi])
            if lista[holi]==1:
                gen=genesexp[holi]
                goi=listo[holi]
                print('Starting '+gen)
                refpath=path+'/gene_info'+'/'+gen
                import os
                if not os.path.exists(refpath):
                    os.makedirs(refpath)
                seqs=extractseq(goi,ref)
                with open(refpath+'/seqs.fasta', 'w') as f:
                    comseq=1
                    for item in seqs:
                        f.write(">"+ gen+ " Seq"+str(comseq)+ "\n" )
                        f.write("%s\n" % item)
                        comseq=comseq+1
                    records=[]
                for ia in range(0,len(goi)-1):
                    records.append(SeqRecord(Seq(seqs[ia]),id=goi[ia]))       
                ie=0
                for seq_record in SeqIO.parse(refpath+"/seqs.fasta", "fasta"):
                    ie=ie+1
                    seq_record = seq_record.upper()
                    targetsall, outfiles = findtargets(seq_record, refpath, ie, outfiles, plp_length, gc_min, gc_max, ligationsite_GC)       
        else:   
            gen=genesexp[holi]
            goi=listo[holi]
            print('Starting '+gen)
            print()
            refpath=path+'/gene_info'+'/'+gen
            import os
            if not os.path.exists(refpath):
                os.makedirs(refpath)
            seqs=extractseq(goi,ref)
            with open(refpath+'/seqs.fasta', 'w') as f:
                for item in seqs:
                    f.write("%s\n" % item)
            clustalw_exe = pathclustal
            cmd = ClustalwCommandline(clustalw_exe,
            infile=refpath+'/seqs.fasta')
            stdout, stderr = cmd()
            alignment = AlignIO.read(refpath+'/seqs.aln', "clustal")
            st=''
            cseqs=[]
            common=[]
            for esa in range(0,alignment.get_alignment_length()):
                un=alignment[:,esa]
                col=collections.Counter(un).most_common(1)[0]
                common.append(col[1])
                if col[1]==len(alignment):
                    st=st+str(col[0])
                else:
                    if len(st)>35:
                        cseqs.append(st)
            #            cseqs=[]
                        st=''
            if len(st)>35:
                        cseqs.append(st)
                        st=''
            freqseq=np.zeros([len(alignment[:,1]),int(alignment.get_alignment_length())])
            for es2 in range(0,alignment.get_alignment_length()):
                        un=alignment[:,es2]
                        for pos2 in range(0,len(un)):
                            if un[pos2]=='-':
                                freqseq[pos2,es2]==0
                            else:
                                cnt=0
                                for el in un:
                                    cnt=cnt+(un[pos2]==el)*1
                                freqseq[pos2,es2]=cnt

            plot_alignment(refpath,alignment,common)
            plot_alignment_of_variants(refpath,freqseq,alignment)
            with open(refpath+'/aligned_seqs.fasta', 'w') as f:
                comseq=1
                for item in cseqs:
                    f.write(">"+ gen+ " Seq"+str(comseq)+ "\n" )
                    f.write("%s\n" % item)
                    comseq=comseq+1
            records=[]
            for ia in range(0,len(goi)-1):
                records.append(SeqRecord(Seq(seqs[ia]),id=goi[ia]))       
            # Loads fasta
            ie=0
            for seq_record in SeqIO.parse(refpath+"/aligned_seqs.fasta", "fasta"):
                seq_record = seq_record.upper()
                #print (seq_record)
                ie=ie+1
                targetsall, outfiles = findtargets(seq_record, refpath, ie, outfiles, plp_length, gc_min, gc_max, ligationsite_GC)
    #           print(targetsall)
    of=pd.DataFrame(outfiles)
    of.to_csv(path+'/outfiles.csv')
    selected,unigene=retrieve_targets(outfiles,path)
    hits=dict(zip(genesexp,listaref))
    selected['exp_hits']=selected['Gene'].map(hits)
    return selected,unigene,notfound

def retrieve_targets(outfiles,path):
    print("Extracting all possible targets for all genes")
    import numpy as np
    import pandas as pd
#    outfiles=pd.DataFrame(outfiles)
#    outfiles.to_csv(path+'/outfiles.csv')
    fil=np.unique(outfiles)
    pan=pd.read_csv(fil[0])
    for s in fil:
       pe=pd.read_csv(s)
       pan=pd.concat([pan,pe])
    selected=pd.DataFrame(pan,columns=pan.columns)
    unigene=pan['Gene'].unique()
    return selected,unigene



def select_sequences(path,selected,genes_required,number_of_selected,subgroup=1):
    '''
    The function takes in five parameters:
    The function select_sequences is intended to filter and select a specified number of 
    sequences from a given dataset based on given genes and other parameters.

    path: The location where the resulting CSV will be saved.
    selected: A DataFrame containing sequences.
    genes_required: A list of genes to be included.
    number_of_selected: The number of sequences to be selected for each gene.
    subgroup: An optional parameter to label the resulting CSV file; it defaults to 1.
    '''
    import pandas as pd
    import random 
    pan=selected
    selected2=pd.DataFrame(columns=selected.columns)
    unigene=genes_required
    for e in unigene:
        print(e)
        ele=pan[pan['Gene']==e]
        if ele.shape[0]<number_of_selected:
            selec=ele
        else:    
            for num in range(0,number_of_selected):
                if ele.shape[0]>0:
                    randomlist = random.sample(range(0, ele.shape[0]), 1)
                    sele=ele.iloc[randomlist,:]
                    try:
                        seleall=pd.concat([seleall,sele])
                    except:
                        seleall=sele
                    exclude=list(range(int(sele['Position']-20),int(sele['Position']+20)))
                    ele=ele.loc[~ele['Position'].isin(exclude),:]
            selec=seleall
        selected2=pd.concat([selected2,selec])
    selected2.to_csv(path+'/selected_targets_group'+str(subgroup)+'.csv')
    return selected2

def check_plps(bcf_all,final_designed,genes,path,subgroup=1):
    '''
    This function checks the output of the mapping function, to exclude the non-specific
    Padlock probes.
    '''
    import pandas as pd
    import numpy as np
    import random
    # Filter bcf_all dataframe to only include rows where the 'Gene' value is in the genes dataframe
    bcf_all2=bcf_all[bcf_all['Gene'].isin(genes['Gene'])]
    # Add a 'same' column that checks if 'exp_hits' is the same as 'number_of_hits' and encodes True as 1 and False as 0
    bcf_all['same']=(bcf_all['exp_hits']==bcf_all['number_of_hits'])*1
    outp=[]
    for s in bcf_all.index:
        r=str(bcf_all.loc[s,'hits']).count(bcf_all.loc[s,'Gene'])
        outp.append(r)
    bcf_all['number_hits_same']=list(outp)
    bcf_all['same2']=(bcf_all['exp_hits']==bcf_all['number_hits_same'])*1
    bcf=bcf_all.loc[bcf_all['same2'].isin([1])]
    bcf=bcf.loc[bcf['same'].isin([1])]
    bcfexc=bcf_all.loc[~bcf_all['same'].isin([1])]
    bcfexc_todesign=bcfexc[~bcfexc['Gene'].isin(bcf['Gene'])]
    bcfexc_todesign=bcfexc_todesign[bcfexc_todesign['Gene'].isin(genes['Gene'])]
    nondesigned=bcf_all2[~bcf_all2['Gene'].isin(bcf['Gene'])]
    ND=nondesigned[nondesigned['Gene'].isin(genes['Gene'])]
    ND.to_csv(path+'/untargeted_genes.csv')
    selected=pd.DataFrame(columns=bcf.columns)
    unigene=bcf['Gene'].unique()
    for e in unigene:
        ele=bcf[bcf['Gene']==e]
        if ele.shape[0]<final_designed:
            sele=ele
        else:    
            randomlist = random.sample(range(0, ele.shape[0]), final_designed) 
            sele=ele.iloc[randomlist,:]
        selected=pd.concat([selected,sele])
    bcf=selected


    bcf.to_csv(path+'/good_targets'+str(subgroup)+'.csv')
    genes_too_low_PLPs=bcf.groupby(['Gene']).count()[bcf.groupby(['Gene']).count()['hits']<final_designed]
    genes_too_low_PLPs=genes_too_low_PLPs.iloc[:,0:1]
    genes_too_low_PLPs.columns=['number_of_PLPs']
    genes_good_PLPs=bcf.groupby(['Gene']).count()[bcf.groupby(['Gene']).count()['hits']>(final_designed-1)]
    genes_good_PLPs=genes_good_PLPs.iloc[:,0:1]
    genes_good_PLPs.columns=['number_of_PLPs']
    genes_too_low_PLPs=bcf.groupby(['Gene']).count()[bcf.groupby(['Gene']).count()['hits']<final_designed]
    genes_too_low_PLPs= genes_too_low_PLPs.iloc[:,0:1]
    genes_too_low_PLPs.columns=['number_of_PLPs']
    genes_no_PLPs=bcfexc_todesign['Gene'].unique()
    bcf.to_csv(path+'/specific_targets_'+str(subgroup)+'.csv')
    return bcf,genes_good_PLPs, genes_too_low_PLPs, genes_no_PLPs

def build_plps(path,specific_seqs_final,L_probe_library,plp_length,how='start',on=201): #starts, end, customized
    import pandas as pd
    import numpy as np
    from Bio.Seq import Seq
    import Bio
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    bcf=specific_seqs_final
    if how=='customized':
#        print(customized)
        selected_ID_to_gene=pd.read_csv(on,sep=',')
        column_names = ["Larm", "idseq", "anchor", "Rarm","LbarID", "AffyID", "Gene"]
        probesP1 = pd.DataFrame (columns = column_names)
        gene_names_ID_columns = ['gene', "idseq", 'Lbar_ID', 'AffyID']
        gene_names_ID = pd.DataFrame(columns=gene_names_ID_columns)
        sbh=pd.read_csv(L_probe_library)
        gene2ID=dict(zip(selected_ID_to_gene['Gene'],selected_ID_to_gene['Lbar_ID']))
        gname = bcf['Gene']
        gname = gname.unique()
        column_names = ["Larm", "idseq", "anchor", "Rarm", "AffyID", "Gene"]
        n=0
        for g in gname:
            gene_names_ID = gene_names_ID.append({"gene": g, "idseq" : np.array(sbh.loc[sbh['Lbar_ID']==gene2ID[g],'ID_Seq'])[0], "Lbar_ID" : str(np.array(sbh.loc[sbh['Lbar_ID']==gene2ID[g],'Lbar_ID'])[0]), "AffyID" : np.array(sbh.loc[sbh['Lbar_ID']==gene2ID[g],'L_Affy_ID'])[0] }, ignore_index=True)
            n=n+1
        gene_names_ID2=gene_names_ID.set_index("gene", drop = False)
        dictiocodes=dict(zip(sbh['L_Affy_ID'],sbh['Barcode_Combi']))
        gene_names_ID2['code']=list(gene_names_ID2['AffyID'].map(dictiocodes))
        gene_names_ID2.to_csv(path+'/codebook.csv')
    if how=='start':
        #this step assings a unique barcode to a gene (meaning that all the probes for the same gene will have the same barcode)
        #generate an empty dataframe to populate with the probe sequence, and useful counters to access stats at the end.
        #ID=LbarID-200
        gname = bcf['Gene']
        gname = gname.unique()
        column_names = ["Larm", "idseq", "anchor", "Rarm","LbarID", "AffyID", "Gene"]
        probesP1 = pd.DataFrame (columns = column_names)
        gene_names_ID_columns = ['gene', "idseq", 'Lbar_ID', 'AffyID']
        gene_names_ID = pd.DataFrame(columns=gene_names_ID_columns)
        ID=on
        sbh=pd.read_csv(L_probe_library)
        n=0
        for g in gname:
            gene_names_ID = gene_names_ID.append({"gene": g, "idseq" : np.array(sbh.loc[sbh['number']==ID+n,'ID_Seq'])[0], "Lbar_ID" : str(np.array(sbh.loc[sbh['number']==ID+n,'Lbar_ID'])[0]), "AffyID" : np.array(sbh.loc[sbh['number']==ID+n,'L_Affy_ID'])[0] }, ignore_index=True)
            n=n+1
        gene_names_ID2=gene_names_ID.set_index("gene", drop = False)
        dictiocodes=dict(zip(sbh['L_Affy_ID'],sbh['Barcode_Combi']))
        gene_names_ID2['code']=list(gene_names_ID2['AffyID'].map(dictiocodes))
        gene_names_ID2.to_csv(path+'/codebook.csv')
    if how=='end':
        #this step assings a unique barcode to a gene (meaning that all the probes for the same gene will have the same barcode)
        #generate an empty dataframe to populate with the probe sequence, and useful counters to access stats at the end.
        #ID=LbarID-200
        gname = bcf['Gene']
        gname = gname.unique()
        column_names = ["Larm", "idseq", "anchor", "Rarm","LbarID", "AffyID", "Gene"]
        probesP1 = pd.DataFrame (columns = column_names)
        gene_names_ID_columns = ['gene', "idseq", 'Lbar_ID', 'AffyID']
        gene_names_ID = pd.DataFrame(columns=gene_names_ID_columns)
        ID=on-len(gname)
        sbh=pd.read_csv(L_probe_library)
        n=0
        for g in gname:
            gene_names_ID = gene_names_ID.append({"gene": g, "idseq" : np.array(sbh.loc[sbh['number']==ID+n,'ID_Seq'])[0], "Lbar_ID" : str(np.array(sbh.loc[sbh['number']==ID+n,'Lbar_ID'])[0]), "AffyID" : np.array(sbh.loc[sbh['number']==ID+n,'L_Affy_ID'])[0] }, ignore_index=True)
            n=n+1
        gene_names_ID2=gene_names_ID.set_index("gene", drop = False)
        dictiocodes=dict(zip(sbh['L_Affy_ID'],sbh['Barcode_Combi']))
        gene_names_ID2['code']=list(gene_names_ID2['AffyID'].map(dictiocodes))
        gene_names_ID2.to_csv(path+'/codebook.csv')
    for index, row in bcf.iterrows():
        r = (row['Gene'])
        x = Seq(row['Sequence'])
        y = x.reverse_complement()
        y = y.upper()
        probesP1 = probesP1.append({"Rarm": str(y[0:round(plp_length/2)]), "Larm": str(y[round(plp_length/2):plp_length]), "anchor" : "TGCGTCTATTTAGTGGAGCC", "idseq" : gene_names_ID2.loc[r]['idseq'], "Lbar_ID" : gene_names_ID2.loc[r]['Lbar_ID'], "AffyID" : gene_names_ID2.loc[r]['AffyID'], "Gene" : gene_names_ID2.loc[r]['gene'] }, ignore_index=True)
        n=n+1

    print ("I just processed",n,"unique target sequences. I am done")
    #The "probes" dataframe at this stage contains the sequences of the probes, before transforming the last base to RNA for ordering 
    probe_col = ["sequence","Lbar_ID", "AffyID", "Gene"]
    probes = pd.DataFrame (columns = probe_col)
    probes["sequence"] = probesP1 ["Larm"] + probesP1 ["idseq"]+ probesP1 ["anchor"]+ probesP1["Rarm"]
    probes["Lbar_ID"] = probesP1 ["Lbar_ID"]
    probes["AffyID"] = probesP1 ["AffyID"]
    probes["Gene"]= probesP1["Gene"]
    # In[ ]: this bit of the script extracts the 3' terminal base of each probe, converts it to RNA code for IDT and rewrites the sequence in the dataframe. Then it outputs a CSV file.
    rnaprobes = []
    for pad in probes.itertuples():
        capt = pad.sequence
        rnabase = (capt[-1])
        plp = (capt [0:-1]+"r"+rnabase)
        rnaprobes.append (plp)
        #print (plp)
    #print (rnaprobes)
    probes["sequence"] = rnaprobes
    #print (probes)
    probes['code']=list(probes['AffyID'].map(dictiocodes))
    probes.to_csv(path+'/designed_PLPs_final.csv')
    return probes

def extract_align_variants(genes, ref, path, pathclustal, selection, plp_length=30, gc_min=50, gc_max=65, ligationsite_GC=True):
    import pandas as pd 
    import os
    from pandas import DataFrame
    import Bio
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    import matplotlib.pyplot as plt
    genes=genes.loc[genes['Gene']==selection,:]
    isExist = os.path.exists(path+'/gene_info')
    if not isExist:
        os.makedirs(path+'/gene_info') 
    outfiles=[]
    genesexp,listoref,lista=extract_seqs(genes,ref)
    listaref=lista
    lista=[]
    listo=[]
    for el in listoref:
        ls=[]
        for ele in el:
            check=ele.find('PREDICTED')==-1
            if check==True:
                ls.append(ele)
        lista.append(len(ls))
        listo.append(ls)
    
    
    for holi in range(0,len(genesexp)):
        ie=0
        notfound=[]
        refpath=path+"/gene_info/"+genesexp[holi]
        if lista[holi]<2:
            if lista[holi]==0:
                print("Gene "+genesexp[holi] +" was not found")
                notfound.append(genesexp[holi])
            if lista[holi]==1:
                gen=genesexp[holi]
                goi=listo[holi]
                print('Starting '+gen)
                refpath=path+'/gene_info'+'/'+gen
                import os
                if not os.path.exists(refpath):
                    os.makedirs(refpath)
                seqs=extractseq(goi,ref)
                with open(refpath+'/seqs_variants.fasta', 'w') as f:
                    comseq=1
                    for item in seqs:
                        f.write(">"+ gen+ " Seq"+str(comseq)+ "\n" )
                        f.write("%s\n" % item)
                        comseq=comseq+1
                    records=[]
                for ia in range(0,len(goi)-1):
                    records.append(SeqRecord(Seq(seqs[ia]),id=goi[ia]))       
                ie=0
                for seq_record in SeqIO.parse(refpath+"/seqs_variants.fasta", "fasta"):
                    ie=ie+1
                    seq_record = seq_record.upper()
                    targetsall, outfiles = findtargets(seq_record, refpath, ie, outfiles, plp_length, gc_min, gc_max, ligationsite_GC)       
        else:   
            gen=genesexp[holi]
            goi=listo[holi]
            print('Starting '+gen)
            print()
            refpath=path+'/gene_info'+'/'+gen
            import os
            if not os.path.exists(refpath):
                os.makedirs(refpath)
            seqs=extractseq(goi,ref)
            with open(refpath+'/seqs_variants.fasta', 'w') as f:
                for item in seqs:
                    f.write("%s\n" % item)
            clustalw_exe = pathclustal
            cmd = ClustalwCommandline(clustalw_exe,
            infile=refpath+'/seqs_variants.fasta')
            stdout, stderr = cmd()
            alignment = AlignIO.read(refpath+'/seqs_variants.aln', "clustal")
            st=''
            cseqs=[]
            common=[]
            for esa in range(0,alignment.get_alignment_length()):
                un=alignment[:,esa]
                col=collections.Counter(un).most_common(1)[0]
                common.append(col[1])
                if col[1]==len(alignment):
                    st=st+str(col[0])
                else:
                    if len(st)>35:
                        cseqs.append(st)
            #            cseqs=[]
                        st=''
            if len(st)>35:
                        cseqs.append(st)
                        st=''
            freqseq=np.zeros([len(alignment[:,1]),int(alignment.get_alignment_length())])
            for es2 in range(0,alignment.get_alignment_length()):
                        un=alignment[:,es2]
                        for pos2 in range(0,len(un)):
                            if un[pos2]=='-':
                                freqseq[pos2,es2]==0
                            else:
                                cnt=0
                                for el in un:
                                    cnt=cnt+(un[pos2]==el)*1
                                freqseq[pos2,es2]=cnt

            plot_alignment(refpath,alignment,common)
            plot_alignment_of_variants(refpath,freqseq,alignment)
    return genesexp,listo,lista

def plot_alignment_of_variants(refpath,freqseq,alignment):
    import pandas as pd 
    import Bio
    import matplotlib.pyplot as plt
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    colors=['r','b','g','c']#'o','p','y','r','b','g','c','o','r','b','g','c','o','r','b','g','c','o','r','b','g','c','o','p','y','r','b','g','c','o','r','b','g','c','o','r','b','g','c','o']
    if np.size(freqseq,axis=0)<4:
        fig, ax = plt.subplots(figsize=(10,(5*(np.size(freqseq,axis=0)))),nrows=np.size(freqseq,axis=0),ncols=1)  
        for s in range(0,np.size(freqseq,axis=0)):
            ax[s].plot(range(0,alignment.get_alignment_length()),freqseq[s,:],c=colors[s])
            ax[s].set_title(alignment[s].id)
    #        ax[s].set_xlab('Bases of the transcript')
    #        ax[s].set_ylab('Numb. of variants matching the base')
            ax[s].hlines(len(alignment)+0.01,linestyles='--',xmin=0,xmax=alignment.get_alignment_length(),colors=[0,0,0])

        plt.savefig(refpath+'/common_regions_though_variants.png')


def extract_seqs_for_variants(path,genesexp,listo,lista,ref,pathclustal):
    import pandas as pd 
    import Bio
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline
    import collections
    from Bio import AlignIO
    import random
    import numpy as np
    outfiles=[]
    for holi in range(0,len(genesexp)):
        ie=0
        notfound=[]
        refpath=path+"/gene_info/"+genesexp[holi]
        if lista[holi]<2:
            if lista[holi]==0:
                print("Gene "+genesexp[holi] +" was not found")
                notfound.append(genesexp[holi])
            if lista[holi]==1:
                gen=genesexp[holi]
                goi=listo[holi]
                print('Starting '+gen)
                refpath=path+'/gene_info'+'/'+gen
                import os
                if not os.path.exists(refpath):
                    os.makedirs(refpath)
                seqs=extractseq(goi,ref)
                with open(refpath+'/seqs.fasta', 'w') as f:
                    comseq=1
                    for item in seqs:
                        f.write(">"+ gen+ " Seq"+str(comseq)+ "\n" )
                        f.write("%s\n" % item)
                        comseq=comseq+1
                    records=[]
                for ia in range(0,len(goi)-1):
                    records.append(SeqRecord(Seq(seqs[ia]),id=goi[ia]))       
                ie=0
                for seq_record in SeqIO.parse(refpath+"/seqs.fasta", "fasta"):
                    ie=ie+1
                    seq_record = seq_record.upper()
                    targetsall, outfiles = findtargets(seq_record, refpath, ie, outfiles, plp_length, gc_min, gc_max, ligationsite_GC)       
        else:   
            gen=genesexp[holi]
            goi=listo[holi]
            print('Starting '+gen)
            print()
            refpath=path+'/gene_info'+'/'+gen
            import os
            if not os.path.exists(refpath):
                os.makedirs(refpath)
            seqs=extractseq(goi,ref)
            with open(refpath+'/seqs.fasta', 'w') as f:
                for item in seqs:
                    f.write("%s\n" % item)
            clustalw_exe = pathclustal
            cmd = ClustalwCommandline(clustalw_exe,
            infile=refpath+'/seqs.fasta')
            stdout, stderr = cmd()
            alignment = AlignIO.read(refpath+'/seqs.aln', "clustal")
            st=''
            cseqs=[]
            common=[]
            for esa in range(0,alignment.get_alignment_length()):
                un=alignment[:,esa]
                col=collections.Counter(un).most_common(1)[0]
                common.append(col[1])
                if col[1]==len(alignment):
                    st=st+str(col[0])
                else:
                    if len(st)>35:
                        cseqs.append(st)
            #            cseqs=[]
                        st=''
            if len(st)>35:
                        cseqs.append(st)
                        st=''
            plot_alignment(refpath,alignment,common)
            with open(refpath+'/aligned_seqs.fasta', 'w') as f:
                comseq=1
                for item in cseqs:
                    f.write(">"+ gen+ " Seq"+str(comseq)+ "\n" )
                    f.write("%s\n" % item)
                    comseq=comseq+1
            records=[]
            for ia in range(0,len(goi)-1):
                records.append(SeqRecord(Seq(seqs[ia]),id=goi[ia]))       
            # Loads fasta
            ie=0
            for seq_record in SeqIO.parse(refpath+"/aligned_seqs.fasta", "fasta"):
                seq_record = seq_record.upper()
                #print (seq_record)
                ie=ie+1
                targetsall,outfiles = findtargets(seq_record,refpath,ie,outfiles)
    #           print(targetsall)
    of=pd.DataFrame(outfiles)
    of.to_csv(path+'/outfiles.csv')
    selected,unigene=retrieve_targets(outfiles,path)
    hits=dict(zip(genesexp,lista))
    selected['exp_hits']=selected['Gene'].map(hits)
    return selected,unigene,notfound

def map_sequences(selected, subgroup, mis_threshold, path, transcriptome=ref):
    import pandas as pd 
    kmers =list(selected['Sequence'])
    #transcriptome = (ref)
    seqlist = []
    hitlist = []
    lenlist = []
    s=0
    for sequence in kmers:
        s=s+1
        mis_threshold = str(mis_threshold)
        print ('Looking for sequence ('+str(s)+'/'+str(len(kmers))+'): '+sequence + ' allowing ' + (mis_threshold) + ' mismatches')
        ########################################MODIFY CUTADAPT PATH IF NEEDED#############################
        import subprocess

        command = [
            "cutadapt", "-j", "0", "-a", str(sequence), "--overlap", "30", 
            "--untrimmed-output", "/dev/null", str(transcriptome), 
            "--no-indels", "-e", str(mis_threshold), "--action=none"
        ]
        print("Executing:", " ".join(command))
        #output= !cutadapt -j 0 -a $sequence --overlap 30 --untrimmed-output /dev/null $transcriptome --no-indels -e $mismatches --action=none 
        n=0
        output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode('utf-8')        
        c2 = [line for line in output.split('\n') if line.startswith('>')]
        print ('Found '+str(len(c2))+' hits')
        seqlist.append (sequence)
        hitlist.append (c2)
        lenlist.append (len(c2))
    expoutput = pd.DataFrame(list(zip(seqlist, hitlist, lenlist)),
                   columns =['sequence_with_Ns', 'hits', 'number_of_hits'])
    bcf_all=pd.concat([selected.reset_index(),expoutput],axis=1)
    bcf_all.to_csv(path+'mapped_sequences'+str(subgroup)+'.csv')
    return bcf_all