def plot_alignment(refpath,alignment,common):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10,3))
    plt.plot(range(0,alignment.get_alignment_length()),common,c='grey')
    plt.hlines(len(alignment)+0.1,linestyles='--',xmin=0,xmax=alignment.get_alignment_length(),colors='r')
    plt.xlim([0,alignment.get_alignment_length()])
    plt.savefig(refpath+'/alignment_of_bases.png')