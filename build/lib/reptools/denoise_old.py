

#read in csv file(s), line by line
    #store (a) sequence, (b) count, (c) file source, (d) length
#send a list of files which are technical replicates
#loop through all files, assigning a unique ID to each CDR3, and collapsing sequences only within techincal replicates
#TODO - pool "Unknown" Vs and Js into more common defined VJs

def denoiseTCR(filelist, outdir, cluster_file, genes = ['V','J','C'], multiplier = 10):
    from builtins import str
    import os
    import io
    import Levenshtein
    import numpy as np
    import itertools as it
    from reptools import FASTAparser, DerepCount, is_derepFas
    from reptools import makelevmatrix, denoise_poolseqs
    flatten = lambda l: [item for sublist in l for item in sublist]
    #
    #
    def checktitleline(fn,genes):
        for (title,seq) in FASTAparser(io.open(fn)):
            break
        try:
            if not is_derepFas(title):
                raise IOError("%s is not a derep'd fasta file." % fn)
        except:
            print(fn)
            raise
        split_title = title.split(';')
        for gene in genes:
            if gene not in [s.split('=')[0] for s in split_title] :
                raise IOError(
                    '%s has an incorrect title line format:\n' 
                    'Gene %s was specified in the function call, '
                    'but is not present in the title line.' % (fn,gene)
                    )
    #
    #
    def getgeneorder(fn,genes):
        gene_indices = []
        for (title,seq) in FASTAparser(io.open(fn)):
            break
        split_title = [t.split('=') for t in title.split(';')]
        for gene in genes:
            gene_found = [True if (title[0]==gene and len(title)==2) else False for title in split_title]
            if sum(gene_found)>1:
                raise IOError('File %s has ambiguous title lines: Gene %s appears more than once.' % (fn,gene))
            else:
                gene_indices.append(gene_found.index(True))
        return(gene_indices)
    
    
    repsCDR3counts = {}
    repsCDR3set = set()
    for fn in filelist:
        fileCDR3set = set()
        if os.path.getsize(fn) == 0 :
            repsCDR3counts[fn]={}
            print(('Warning:\n%s is an empty file.\n' % fn)) #ignore empty files
        else:
            checktitleline(fn,genes) #check the title line structure
            geneorder = getgeneorder(fn,genes) #get the gene order in the title line
            repsCDR3counts[fn]={}
            oldfileCDR3set_size = 0
            #
            for (title,seq) in FASTAparser(io.open(fn)):
                splittitle = title.split(';')
                genes_present = [splittitle[i] for i in geneorder]
                clone = tuple([seq]+[g for g in it.chain(genes_present)])
                fileCDR3set.add(clone)
                newCDR3set_size = len(fileCDR3set)
                if newCDR3set_size == oldfileCDR3set_size:
                    raise ValueError(
                            'The same clone is present twice in %s.  '
                            'Has this file been properly dereplicated?'
                            % fn
                            )
                else:
                    oldfileCDR3set_size = newCDR3set_size
                repsCDR3counts[fn][clone] = DerepCount(title)
            repsCDR3set.update(fileCDR3set)
    #
    #build a numpy array
    repsCDR3list = sorted(list(repsCDR3set)) #to ensure stable ordering
    repsCDR3lookup = {CDR3:n for n,CDR3 in enumerate(repsCDR3list)}
    counts = np.zeros(shape=(len(filelist),len(repsCDR3list)),dtype='int')
    #
    for n,fn in enumerate(filelist): #fill a 2D numpy array with the counts.  Rows=files, columns=CDR3
        for p,CDR3 in enumerate(repsCDR3list):
            try: #try-except rather than key search, to attempt speed-up 
                counts[n,p] = repsCDR3counts[fn][CDR3]
            except KeyError:
                pass
    #
    #to build clonotypes, combine the replicates [this can now be done just by summing columns in the numpy array]
    #
    repsCDR3lengths = np.array([len(CDR3[0]) for CDR3 in repsCDR3list])
    roundedCDR3lengths = np.around(repsCDR3lengths/3.)*3
    #
    #find all gene/length combos
    combos = {}
    for CDR3 in repsCDR3list:
        try:
            combos[(CDR3[1:],int(roundedCDR3lengths[repsCDR3lookup[CDR3]]))].append(CDR3) 
        except KeyError:
            combos[(CDR3[1:],int(roundedCDR3lengths[repsCDR3lookup[CDR3]]))] = [CDR3] 
    #
    #TODO - use weighted Levenshtein instead, if it isn't too slow
    #
    counts_out = np.copy(counts)
    #
    changesmatrixCDR3s = {}
    for combo in combos:
        seqs_selected = combos[combo]#
        seqindices = [repsCDR3lookup[k] for k in seqs_selected]
        levmatrix = makelevmatrix(seqs_selected)
        counts_np_selected = counts[:,seqindices] 
        counts_out_selected = counts_out[:,seqindices]  #creates a copy, not a view
        #start with the rarest sequences, and count upwards - roll-up the output into input after each iteration
        ascendingcountorder = np.copy(np.argsort(np.sum(counts_np_selected,axis=0) )) #the copy seems to be needed to make the next line work
        changematrix = {n:set([n]) for n in range(len(seqindices))} #initial change matrix where each cluster contains only one CDR3
        #seqsexamined = [] 
        for i in ascendingcountorder:
            changematrix = denoise_poolseqs(1,multiplier,levmatrix,i,counts_np_selected,counts_out_selected,changematrix)
            counts_np_selected = np.copy(counts_out_selected) 
            #seqsexamined.append(i) 
            #ascendingcountorder = np.copy(np.argsort(np.sum(counts,axis=0) ))
            #ascendingcountorder = ascendingcountorder[(np.in1d(ascendingcountorder,seqsexamined,invert=True))] 
        
        counts_out[:,seqindices] = counts_out_selected #copy the results back to counts_out
        #convert changematrix indexing
        reindexedchanges = {seqindices[n]:list([seqindices[p] for p in changematrix[n]]) for n in list(changematrix.keys())}
        #
        changesmatrixCDR3s.update(reindexedchanges)  #output the change matrix
    #
    #
    #write output fasta
    for n,fn in enumerate(filelist): #loop through the files=rows of numpy array
        outfile = os.path.join(outdir,os.path.split(fn)[1])
        with io.open(outfile,'w') as fasta_handle:
            for n,(CDR3,count) in enumerate(zip(repsCDR3list,counts_out[n,:])):
                if count>0:
                    fasta_handle.write('>%s;%s;size=%s;\n%s\n' % (str(n), ';'.join(CDR3[1:]), str(count), CDR3[0]))
    #
    #
    #save the cluster information from changesmatrixCDR3s
    with io.open(cluster_file,'w') as fn:
        for n in changesmatrixCDR3s:
            clustered = [repsCDR3list[ID] for ID in sorted(changesmatrixCDR3s[n])]
            if len(clustered)>0:
                fn.write('#%s\n%s\n' % (str(n),'\n'.join(['%s\t%s\t' % ('\t'.join(clone[1:]),clone[0] ) for clone in clustered])))


def makelevmatrix(CDR3s): #CDR3s to be tuples, with item [0] a sequence string
    import itertools as it#v12 change (to speed up)
    import Levenshtein
    import numpy as np
    levmatrix = np.empty(shape=(len(CDR3s),len(CDR3s)),dtype=int)
    for n in range(len(CDR3s)):
        levmatrix[n,n] = 0 #v12 change (n/p equality not produced in loop, as I'm now using itertools)
    for n,p in it.combinations(list(range(len(CDR3s))),2):  #v12 change (to speed up)
        levmatrix[n,p] = Levenshtein.distance(CDR3s[n][0],CDR3s[p][0])
        levmatrix[p,n] = levmatrix[n,p]
    return(levmatrix)


def denoise_poolseqs(distance,mult,levmatrix,i,counts_np_selected,counts_out_selected,changematrix):
    import numpy as np
    if np.sum(counts_np_selected[:,i])==0:     #no counts here
        #print('returning unchanged')
        return(changematrix)
    dist_np = np.where(levmatrix[i,:]==distance)[0]
    if dist_np.shape[0]==0:
        return(changematrix) #no sequences at distance
    counts_at_dist = counts_np_selected[:,dist_np] #creates a copy
    counts_out_at_dist = counts_out_selected[:,dist_np] #as this is selecting with an index array, it creates a copy, not a view
    dist_comparator_array = np.sum(counts_at_dist,axis=0,dtype='float')/np.sum(counts_np_selected[:,i])
    if np.amax(dist_comparator_array) >= mult:
        target = np.argmax(dist_comparator_array) #note that ties are broken blindly (take the first instance), rather than using an intelligent algorithm
        counts_out_at_dist[:,target] = counts_at_dist[:,target] + counts_np_selected[:,i] #add the donor reads to the recipient
        counts_out_selected[:,i] = 0 #delete the donor reads
        counts_out_selected[:,dist_np] = counts_out_at_dist #now assign back to counts_out_selected
        changematrix[dist_np[target]].update(changematrix[i]) #modify the changematrix
        changematrix[dist_np[target]].add(i) 
        changematrix[i] = set([])
        return(changematrix)
    else:
        return(changematrix)


