from numba import vectorize
from numba import njit
import numpy as np
import collections
import reptools

def derep_FASTQ(fn,clust_file):
    seqs = collections.defaultdict(dict)
    gene_ids = collections.defaultdict(list) # stored in a list to avoid repetition - slower, but saves memory
    changes = collections.defaultdict(list)
    with open(fn) as in_handle:
        for title,seq,qual in reptools.FASTQparser(in_handle):
            if any([nt not in ['A','T','G','C','a','t','g','c'] for nt in seq.strip()]):
                continue #omit ambiguous sequences
            title_list = title.strip().strip(';').split(';')
            seqlen = len(seq)
            #seqprobs = [-float(Q)/10 for Q in [ord(c)-33 for c in qual]] #for logs
            seqprobs = [10**(-float(Q)/10) for Q in [ord(c)-33 for c in qual]]
            if title_list[1:] not in gene_ids[seqlen]: 
                #if the geneid is new, the sequence is new (by definition)
                gene_ids[seqlen].append(title_list[1:])
                seqs[seqlen][(seq,len(gene_ids[seqlen])-1)] = [
                                                1, #1, because this is the first time this sequence has been found
                                                title_list[0], #the first title found is stored for output
                                                seqprobs
                                               ] 
            elif (seq,gene_ids[seqlen].index(title_list[1:])) not in seqs[seqlen]: 
                #if the sequence string is new, but not geneid, the geneid code can be reused
                seqs[seqlen][(seq,gene_ids[seqlen].index(title_list[1:]))] = [
                                                              1,
                                                              title_list[0],
                                                              seqprobs
                                                             ]
            else: 
                #not a unique sequence
                #add to changes
                changes[seqs[seqlen][(seq,gene_ids[seqlen].index(title_list[1:]))][1]].append(title_list[0])
                #calculate new probs list
                newprobs = [
                            old * new for old,new in zip( 
                                                          seqs[seqlen][(seq,gene_ids[seqlen].index(title_list[1:]))][2],
                                                          seqprobs
                                                         )
                            ]
                seqs[seqlen][(seq,gene_ids[seqlen].index(title_list[1:]))][0] += 1
                seqs[seqlen][(seq,gene_ids[seqlen].index(title_list[1:]))][2] = newprobs
                
    #very high counts can result in probs of zero (float underrun), which is invalid; so change them to Phred=120
    #TODO = change the probability handling to working with the log probabilities, which will also save
    #memory, as I can then use float16 - N.B. I can't, because numba doesn't support float16 (yet)
    for seqlen in seqs:
        for k in seqs[seqlen]:
            if min(seqs[seqlen][k][2]) < 0.000000000001:
                seqs[seqlen][k][2] = [prob if prob >=0.000000000001 else 0.000000000001 for prob in seqs[seqlen][k][2]]
    
    if clust_file:
        with open(clust_file,'w') as clust_handle:
            for recipient in changes:
                clust_handle.write('{}\t{}\n'.format(recipient,'\t'.join(changes[recipient])))
    
    return(dict(seqs),dict(gene_ids))


def process_geneids(gene_ids):
    """
    #geneids is a dictionary keyed by seqlen, with each dictionary containing a list of ids.  The ids themselves are
    # a list of strings, one for each gene
    #returns a dictionary of lists of sets of allelles (one set per gene, one list per entry in the input list).
    """
    new_gene_ids = {}
    for seqlen in gene_ids:
        new_gene_ids[seqlen] = {}
        for n,gene_list in enumerate(gene_ids[seqlen]):
            new_gene_ids[seqlen][n] = [set(gene.split('=')[1].split('+')) for gene in gene_list]
    gene_labels = [gene.split('=')[0] for gene in gene_list]
    return(new_gene_ids,gene_labels)


def build_arrays(subseqs):
    """
    Input: a dictionary containing sequence info from sequences of the same length, but different genes
    Returns:
        a list of numpy arrays:
            [seq_array, gene_id_array, prob_array, counts_array, seq_names
    """
    #build numpy arrays
    seq_array = reptools.build_seq_array([k[0] for k in subseqs])
    gene_ids = [k[1] for k in subseqs]
    if max(gene_ids) > 65535: raise ValueError('Too many unique gene id combinations (>65535)')
    gene_id_array = np.array([k[1] for k in subseqs],dtype=np.uint16)
    prob_array = np.array([subseqs[k][2] for k in subseqs],dtype=np.float32)
    counts_array = np.array([subseqs[k][0] for k in subseqs])
    seq_names = np.array([subseqs[k][1] for k in subseqs])
    
    #find out how big an integer COULD be needed for counts (if all counts are assigned to just one read)
    counts_dtype = np.sum(counts_array).dtype
    counts_array = counts_array.astype(counts_dtype)
    
    return([seq_array, gene_id_array, prob_array, counts_array, seq_names])


def build_seq_array(seqlist,dtype=np.int8):
    """
    all sequences must be of equal length
    """
    char_array = np.empty([len(seqlist),len(seqlist[0])], dtype=dtype)
    
    for n,seq in enumerate(seqlist):
        char_array[n,:] = [ord(char) for char in seq]
    
    return(char_array)


def update_counts_dtype(seqs_np,counts_idx=3):
    """
    check that the counts array dtype is always big enough to cope with pooling from adjacent seqlen
    and rewrite with changed dtype if it isn't
    """
    previous_total = 0
    total = None
    for n,seqlen in enumerate(sorted(seqs_np)):
        if total is None:
            total = np.sum(seqs_np[seqlen][counts_idx])
        else:
            total = next_total
        try:
            next_total = np.sum(seqs_np[sorted(seqs_np)[n+1]][counts_idx])
        except IndexError: #when there is no next value (i.e., final item)
            next_total = 0
        desired_dtype = np.sum(previous_total+total+next_total).dtype
        seqs_np[seqlen][counts_idx] = seqs_np[seqlen][counts_idx].astype(desired_dtype)
        previous_total = total
        #return(seqs_np)


#now goes from most to least common
def denoise_substitutions(seqs_np,gene_ids_sets,threshold,clust_file,weight_by_qual=True,seqlens=False):
    """
    denoise reads of same length:
        where a sequence is (hamming distance)*threshold less numerous than another, add its reads to the more common
        sequence IF the two sequences have at least one gene segment ID in common for each gene in the id line.
        The daughter record will keep the gene segment IDs of the more numerous parent.
        The rarer parent will have its read count set to zero.
        If weight_by_qual=True (the default), adjusts hamming distance according to qual scores.
        TODO: daughter qual score is the psoterior probability of the input qual scores (weighted by read count)
        Input:
        seqs_np = a dictionary of lists of numpy arrays, one entry per sequence length, with the value being a list
        of numpy arrays for: sequences, gene id codes, probabilities, rea counts, sequence names
        gene_ids_sets = a dictionary of dictionaries, one entry per sequence length, with the value being a dictionary
        where keys match the integers in the gene id codes numpy array, and the values are a list of sets, giving the
        gene segments identified for each gene (one set per gene).
        Output:
            None.  Numpy arrays are modified in place
    """
    if not seqlens: #if False, process all
        seqlens = list(seqs_np.keys())
    seqlens = [ln for ln in seqlens if ln >1] #single nt sequences break the array creation code, and are not useful
    with open(clust_file,'w') if clust_file else reptools.dummy_context_mgr() as clust_handle:
        for seqlen in seqlens:
            changes = {}
            seq_array, gene_id_array, prob_array, counts_array, seq_names  = seqs_np[seqlen]
            changes_made=True
            #loopcounter=0
            while changes_made: #iterate until no more improvements
                changes_made=False
                #loopcounter+=1
                #print('Seqlen={}, iteration={}'.format(str(seqlen),str(loopcounter)))
                reptools.sort_by_freq(
                               seq_array,
                               gene_id_array,
                               prob_array,
                               counts_array,
                               seq_names
                             )
                if np.any(counts_array==0):                                                                  
                    firstzero = np.where(counts_array==0)[0][0]
                else:
                    firstzero = seq_array.shape[0]
                np.seterr(over='ignore') #to not report overflow errors - see below
                for n in range(0,firstzero): #iterate from most to least common (but non-zero) 
                    matching = np.equal(seq_array[n+1:,:],seq_array[n,:]) #compare this row with all subsequent rows
                    if weight_by_qual:
                        #weight by probabilities
                        weights = generic_chance_of_miss(prob_array[n,:],prob_array[n+1:,:],matching[:,:]) #TODO - change generic_chance_of_miss() for logs
                        dists = np.sum(weights,axis=1) #gives hamming distances weighted by qual scores
                    else:
                        dists = np.sum(np.invert(matching),axis=1)
                    #find reads which are "close enough"
                    #the next line may cause overflows, but if the the value is too high, it will be insanely large,  
                    #and set to Inf, so the comparison should still work
                    hits = np.where( counts_array[n+1:] >= np.power(threshold,dists,dtype=np.float32) * counts_array[n] ) 
                    #loop through the hits, checking that the genes match, and continue until a hit is found where they do.
                    for hit in hits[0]:
                        #check gene ids
                        query_genes = gene_ids_sets[seqlen][gene_id_array[n]]
                        target_genes = gene_ids_sets[seqlen][gene_id_array[n+1+hit]]
                        if sum([len(_q.union(_t))==len(_q)+len(_t) for _q,_t in zip(query_genes,target_genes)]) == 0:
                            #all genes have at least one allele in common
                            #v0.14.1: modify qual score of target
                            # calculate prob of base in source ACTUALLY being base in target
                            #   (1-prob_array[n])/3 because a 1/3 chance of each uncalled base
                            source_prob_misread = ((1-prob_array[n])/3)**counts_array[n]  #TODO - change for logs
                            #multiply target prob by the calculated source misread prob
                            prob_array[n+1+hit] = prob_array[n+1+hit]*source_prob_misread #TODO - change for logs
                            #add counts from query to target
                            counts_array[n+1+hit] += counts_array[n]  
                            #set query counts to zero
                            counts_array[n] = 0 
                            #record changes
                            try:
                                changes[seq_names[n+1+hit]].append(seq_names[n])
                            except KeyError:
                                changes[seq_names[n+1+hit]]=[seq_names[n]]
                            changes_made=True
                            break
                np.seterr(over='warn')
            #print('Seqlen {} denoised'.format(seqlen,loopcounter))
            reptools.save_changes(clust_handle,changes,seq_names,counts_array) 
    return #numpy arrays will have been modified in-place, and gene_ids,gene_ids_sets have not been modified


#@profile
def denoise_indelonly(seqs_np,gene_ids_sets,threshold,clust_file,seqlens=False):
    """
    denoise reads differing length, removing indels only (and no more than one indel):
        where a sequence is threshold less numerous than another and has only a single indel difference,
        add its reads to the more common sequence IF the two sequences have at least one gene segment ID in common
        for each gene in the id line.
        The daughter record will keep the gene segment IDs of the more numerous parent.
        The rarer parent will have its read count set to zero.
        TODO: Is there a way to modify qual socres with indels?  Nothing obvious (perhaps with insertions?)
        Input:
        seqs_np = a dictionary of lists of numpy arrays, one entry per sequence length, with the value being a list
        of numpy arrays for: sequences, gene id codes, probabilities, rea counts, sequence names
        gene_ids_sets = a dictionary of dictionaries, one entry per sequence length, with the value being a dictionary
        where keys match the integers in the gene id codes numpy array, and the values are a list of sets, giving the
        gene segments identified for each gene (one set per gene).
        Output:
            None.  Numpy arrays are modified in place
    """
    if not seqlens:
        seqlens = list(seqs_np.keys())
    seqlens = [ln for ln in seqlens if ln >1] #single nt sequences break the array creation code, and are not useful
    seqlens = sorted(seqlens)
    if len(seqlens)<2:
        print('indel denoising cannot be performed when all sequences are the same length')
        return
    previous_seq_array = None
    with open(clust_file,'w') if clust_file else reptools.dummy_context_mgr() as clust_handle:
        for seqlen in seqlens:
            seq_array, gene_id_array, prob_array, counts_array, seq_names  = seqs_np[seqlen]
            changes = {}
            if previous_seq_array is not None:
                if seqlen - previous_seqlen==1: #only look for indels when seqlength delta==1
                    #the previous seq_array is one base shorter than the present one
                    #so want to mask each position in the current one in turn, and get the hamming distance
                    changes_made=True
                    #loopcounter=0
                    while changes_made: #iterate until no more improvements
                        changes_made=False
                        #loopcounter+=1
                        #print('Seqlen={}, iteration={}'.format(str(seqlen),str(loopcounter)))
                        reptools.sort_by_freq(
                                       seq_array,
                                       gene_id_array,
                                       prob_array,
                                       counts_array,
                                       seq_names
                                     )
                        if np.any(counts_array==0):                                                                  
                            firstzero = np.where(counts_array==0)[0][0]
                        else:
                            firstzero = seq_array.shape[0]
                        for n in range(0,firstzero): #iterate from most to least common (but non-zero) 
                            indels = np.full(previous_seq_array.shape[0],False)
                            b = np.full(seq_array.shape[1],True)
                            for x in range(seq_array.shape[1]): #loop through the columns, sliding the missing column across
                                c = np.copy(b)
                                c[x] = False #set missing column
                                #look for perfect matches (with the missing column excluded)
                                indels = find_indels(indels,previous_seq_array,seq_array[n,c])
                            indels = np.logical_and(indels,previous_counts_array>0) #exclude read zero seqs
                            if np.any(indels):
                                indels_where = np.where(indels)[0]
                                cur_into_prev_ratios = previous_counts_array[indels_where]/(counts_array[n]*threshold) #rewrite to avoid division?
                                prev_into_cur_ratios = counts_array[n]/(previous_counts_array[indels_where]*threshold)
                                #these two arrays will index indels_where
                                #a value of >1 meets the threshold
                                #check genes for all where >1
                                gene_matches = [False]*len(indels_where)
                                for p1,p2 in enumerate(indels_where):
                                    qry_genes = gene_ids_sets[seqlen][gene_id_array[n]]
                                    targ_genes = gene_ids_sets[previous_seqlen][previous_gene_id_array[p2]]
                                    if sum([len(_q.union(_t))==len(_q)+len(_t) for _q,_t in zip(qry_genes,targ_genes)]) == 0:
                                        gene_matches[p1] = True #gene_matches will reference indels_where, and also the ratios
                                if np.any(gene_matches):
                                    cur_into_prev_ratios = cur_into_prev_ratios * gene_matches #set ratios to 0 where no match
                                    prev_into_cur_ratios = prev_into_cur_ratios * gene_matches
                                    cur_into_prev_max_idx = np.argmax(cur_into_prev_ratios) #get best ratio
                                    prev_into_cur_max_idx =  np.argmax(prev_into_cur_ratios)
                                    #which direction do we prefer to move the reads?
                                    if (
                                            cur_into_prev_ratios[cur_into_prev_max_idx]
                                            >
                                            prev_into_cur_ratios[prev_into_cur_max_idx]
                                        ):
                                        if cur_into_prev_ratios[cur_into_prev_max_idx] >= 1:
                                            #print((cur_into_prev_ratios[cur_into_prev_max_idx]))
                                            targetrow = indels_where[cur_into_prev_max_idx]
                                            previous_counts_array[targetrow] += counts_array[n]
                                            counts_array[n] = 0
                                            try:
                                                changes[previous_seq_names[targetrow]].append(seq_names[n])   
                                            except KeyError:
                                                changes[previous_seq_names[targetrow]] = [seq_names[n]]                            
                                            changes_made = True
                                    else:
                                        if prev_into_cur_ratios[prev_into_cur_max_idx] >= 1:
                                            #print((prev_into_cur_ratios[prev_into_cur_max_idx]))
                                            targetrow = indels_where[prev_into_cur_max_idx]
                                            counts_array[n] += previous_counts_array[targetrow]
                                            previous_counts_array[targetrow] = 0
                                            try:
                                                changes[seq_names[n]].append(previous_seq_names[targetrow])
                                            except KeyError:
                                                changes[seq_names[n]] = [previous_seq_names[targetrow]]
                                            changes_made = True
                reptools.save_changes(clust_handle,previous_changes,previous_seq_names,previous_counts_array) 
                
                reptools.sort_by_freq( #resort, to allow omission of zeros
                           seq_array,
                           gene_id_array,
                           prob_array,
                           counts_array,
                           seq_names
                           )
            
            if np.any(counts_array==0):                                                                  
                firstzero = np.where(counts_array==0)[0][0]
            else:
                firstzero = seq_array.shape[0]
            
            previous_seqlen = seqlen
            #the next four lines previously used np.copy.  I don't think this is necessary, and is not desirable for modify
            #in place
            previous_seq_array = seq_array[0:firstzero,:]
            previous_gene_id_array = gene_id_array[0:firstzero]
            previous_counts_array = counts_array[0:firstzero]
            previous_seq_names = seq_names[0:firstzero]
            previous_changes = changes
        #output: save current changes (because they will never be saved as the previous changes)
        reptools.save_changes(clust_handle,changes,seq_names,counts_array) 
        
    return #modification should have occurred in place


#@profile
def simplify_genes(seqs_np,gene_ids_sets,clust_file,seqlens=False):
    """
    remove unnecessary ambiguities in gene segment IDs:
        where a sequence is identical to another, the rarer sequence is combined with the more common, IF the two
        sequences have at least one gene segment ID in common for each gene in the id line.
        The daughter record will keep the gene segment IDs of the more numerous parent.
        The rarer parent will have its read count set to zero.
        If the two records have an qual read count, the simpler (i.e. shorter) set of possible gene segments is taken
        for each gene (so that the output selection may take some genes from one sequence, and some from the other).
        TODO: daughter qual score is the posterior probability of the input qual scores (weighted by read count)
        Input:
        seqs_np = a dictionary of lists of numpy arrays, one entry per sequence length, with the value being a list
        of numpy arrays for: sequences, gene id codes, probabilities, rea counts, sequence names
        gene_ids_sets = a dictionary of dictionaries, one entry per sequence length, with the value being a dictionary
        where keys match the integers in the gene id codes numpy array, and the values are a list of sets, giving the
        gene segments identified for each gene (one set per gene).
        Output:
            gene_ids_sets.  Modified gene set dictionary.  Numpy arrays are modified in place
    """
    if not seqlens:
        seqlens = list(seqs_np.keys())
    seqlens = [ln for ln in seqlens if ln >1] #single nt sequences break the array creation code, and are not useful
    with open(clust_file,'w') if clust_file else reptools.dummy_context_mgr() as clust_handle:
        for seqlen in seqlens:
            #print('Seqlen={}'.format(str(seqlen)))
            changes = {}
            seq_array, gene_id_array, prob_array, counts_array, seq_names  = seqs_np[seqlen]
            reptools.sort_by_freq(
                           seq_array,
                           gene_id_array,
                           prob_array,
                           counts_array,
                           seq_names
                        )
            #iterate from least to most common (missing the very most common, as there will be nothing to compare it with)
            for n in range(len(counts_array)-1,0,-1):
                #compare this row with all at least as common (i.e. earlier) rows
                identical = np.where(identical_rows(seq_array[:n,],seq_array[n,:]))[0]
                #work through the identical rows, from the first (most common) on, checking for gene identity
                if len(identical)>0:
                    for p in identical:
                        qry_genes = gene_ids_sets[seqlen][gene_id_array[n]]
                        targ_genes = gene_ids_sets[seqlen][gene_id_array[p]]
                        #if there is a match
                        if sum([len(_q.union(_t))==len(_q)+len(_t) for _q,_t in zip(qry_genes,targ_genes)]) == 0:
                            #check if the read count is equal (should never be less)
                            if counts_array[p] == counts_array[n]:
                                #if so, take the simplest gene set for each
                                newgeneset = []
                                for seg in range(len(qry_genes)):
                                    if len(qry_genes[seg]) < len(targ_genes[seg]):
                                        newgeneset.append(qry_genes[seg])
                                    else:
                                        newgeneset.append(targ_genes[seg])
                                #is this combination of gene ids already in the dictionary?
                                if newgeneset in list(gene_ids_sets[seqlen].values()): 
                                    gene_id_array[p] = [
                                                        k for k in gene_ids_sets[seqlen] if gene_ids_sets[seqlen][k]==newgeneset
                                                        ][0] #if so, set gene_id code
                                else: #if not, create an entry
                                    newentry = max(gene_ids_sets[seqlen])+1
                                    if newentry > 65535:
                                        raise ValueError('Too many unique gene id combinations (>65535)')
                                    gene_ids_sets[seqlen][newentry] = newgeneset
                                    gene_id_array[p] = newentry
                            
                            #combine read counts
                            counts_array[p] += counts_array[n]
                            counts_array[n] = 0
                            #record changes
                            try:
                                changes[seq_names[p]].append(seq_names[n])
                            except KeyError:
                                changes[seq_names[p]]=[seq_names[n]]
                            changes_made=True
                            break #if a match was found, continue to the next rarest row
            reptools.save_changes(clust_handle,changes,seq_names,counts_array) 
    return(gene_ids_sets)


#sort by descending frequency
def sort_by_freq(seq_array, gene_id_array, prob_array, counts_array, seq_names):
    sort_order = np.argsort(counts_array)[::-1]
    seq_array[:,:] = seq_array[sort_order]
    gene_id_array[:] = gene_id_array[sort_order]
    prob_array[:,:] = prob_array[sort_order]
    counts_array[:] = counts_array[sort_order]
    seq_names[:] = seq_names[sort_order]
    return #modification in-place


@vectorize
def generic_chance_of_miss(query,target,matching):
    return((9 - 3*query - 3*target + 4*query*target + matching * (-9 + 12*query + 12*target - 16*query*target))/ 9) ##TODO - change for logs


@njit
def identical_rows(previous_seq_array,seq_array_line):
    not_equal_array =  np.not_equal(previous_seq_array,seq_array_line)
    a = np.full(not_equal_array.shape[0],True,dtype=np.bool_)
    # short-circuiting replacement for np.all(,axis=1)
    for row in range(not_equal_array.shape[0]):
        for x in not_equal_array[row,:]:
            if x:
                a[row] = False
                break 
    return(a)


@njit
def find_indels(indels,previous_seq_array,seq_array_line):
    return(
           np.logical_or(
                           indels,
                           reptools.identical_rows(previous_seq_array,seq_array_line)
                        )
    )


def save_changes(listhandle,changes,seq_names,counts_array):
    """
    Saves cluster information in a simple tab-delimited format, where the first column in the name of the retained
    sequence, and subsequent columns contain the names of sequences merged with it.
    """
    non_zero_seqnames = set(seq_names[counts_array>0].flat)
    clustered_seqs = set([k for k in changes])
    unclustered_seqs = non_zero_seqnames - clustered_seqs
    for seqname in clustered_seqs:
        listhandle.write('{}\t{}\n'.format(seqname,'\t'.join(changes[seqname])))
    for seqname in unclustered_seqs:
        pass
        #listhandle.write('{}\n'.format(seqname))


def saveFASTX(seqs_np,gene_ids_sets,FASTAout,FASTQout,gene_labels):
    #make gene_strings
    gene_ids = make_gene_strings(gene_ids_sets,gene_labels)
    #find maximum count
    maxcount = np.amax([np.amax(seqs_np[seqlen][3]) for seqlen in seqs_np]) 
    #descend through counts
    with open(FASTAout,'w') if FASTAout else reptools.dummy_context_mgr() as fasta_handle:
        with open(FASTQout,'w') if FASTQout else reptools.dummy_context_mgr() as fastq_handle:
            for count in range(maxcount,0,-1): #don't want to write where count=0
                for seqlen in seqs_np:
                    towrite = np.where(seqs_np[seqlen][3]==count)[0]
                    for x in towrite:
                        #try:
                        title = '{};{};size={}'.format(
                                                       seqs_np[seqlen][4][x],
                                                       gene_ids[seqlen][seqs_np[seqlen][1][x]],
                                                       seqs_np[seqlen][3][x]
                                                       )
                        #except TypeError: #if there are no seqs with this count, the iterator returns an empty numpy
                        #    #array, which breaks the indexing
                        #   continue
                        seq = ''.join( [ chr(c) for c in seqs_np[seqlen][0][x] ] )
                        try:
                            qual = ''.join(
                                              [
                                                 chr(int(c+33))
                                                 if c <= 93
                                                 else chr(126)
                                                 for c in #prob_toqual(seqs_np[seqlen][2][x])
                                                         [
                                                           np.around(np.multiply(np.log10(prob),-10),decimals=0) ##TODO - change for logs
                                                              for prob in seqs_np[seqlen][2][x] #TODO - change for logs
                                                          ]
                                               ]
                                           )
                        except:
                            print(seqlen)
                            print(x)
                            #print(seqs_np[seqlen][2])
                            print((seqs_np[seqlen][2][x]))
                            raise
                        fasta_handle.write('>{}\n{}\n'.format(title,seq))
                        fastq_handle.write('@{}\n{}\n+\n{}\n'.format(title,seq,qual))


#@vectorize(["int16(float32)"])
#def prob_toqual(prob):
#    np.around(np.multiply(np.log10(prob),-10),decimals=0)


def make_gene_strings(gene_ids_sets,gene_labels):
    gene_ids = {}
    for seqlen in gene_ids_sets:
        gene_ids[seqlen] = {}
        for id in gene_ids_sets[seqlen]:
            gene_ids[seqlen][id] = ';'.join([
                                             '{}={}'.format(
                                                            label,
                                                            '+'.join(
                                                                     [seg for seg in gene]
                                                                    )
                                                            )
                                              for gene,label in zip(gene_ids_sets[seqlen][id],gene_labels)
                                             ]
                                             )
    return(gene_ids)


