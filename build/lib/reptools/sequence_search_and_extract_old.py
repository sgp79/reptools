import tempfile
import os
import reptools

def retrieve_tophit(hitsfile,filetype,mincols=False,evalue=False,title_split_char='',useScore=False):
    """
    Returns the hit(s) with the lowest evalue, optionally filtering by minimum alignment span and maximum evalue.
    If evalue is not provided (or if useScore==True), takes hit(s) with highest score
    Returns a single hsp object, breaking ties.
    """
    tophits = {}
    #first, check that evalues are present, if required
    if not useScore:
        #first, find a result with a hit
        for result in reptools.hitsParser(hitsfile,filetype,mincols=mincols,evalue=evalue):
            if len(result.hits)>0:
                if len(result.hits[0].hsps)>0:
                    if result.hits[0].hsps[0].evalue is None:
                        useScore = True
                    break
    
    for result in reptools.hitsParser(hitsfile,filetype,mincols=mincols,evalue=evalue):
        #loop through the hits, recording the hit names
        minevalue=False
        maxscore=0
        for hit in result.hits:
            for hsp in hit.hsps:
                if useScore: #if we're looking at bitscore
                    if hsp.bitscore>maxscore:
                        tophits[result.id.split(title_split_char)[0]] = hsp
                        maxscore = hsp.bitscore
                    else:
                        pass
                else: #if we're looking at evalue
                    if hsp.evalue<minevalue or minevalue is False:
                        tophits[result.id.split(title_split_char)[0]] = hsp
                        minevalue = hsp.evalue
                    else:
                        pass
    return(tophits)


def retrieve_tophit_names(hitsfile,filetype,mincols=False,evalue=False,title_split_char='',useScore=False):
    """
    Returns the hit(s) with the lowest evalue, optionally filtering by minimum alignment span and maximum evalue.
    If evalue is not provided (or if useScore==True), takes hit(s) with highest score
    """
    hitnames = {}
    #first, check that evalues are present, if required
    if not useScore:
        #first, find a result with a hit
        for result in reptools.hitsParser(hitsfile,filetype,mincols=mincols,evalue=evalue):
            if len(result.hits)>0:
                if len(result.hits[0].hsps)>0:
                    if result.hits[0].hsps[0].evalue is None:
                        useScore = True
                    break
    
    for result in reptools.hitsParser(hitsfile,filetype,mincols=mincols,evalue=evalue):
        hitnames[result.id.split(title_split_char)[0]]=set()
        #loop through the hits, recording the hit names
        minevalue=False
        maxscore=0
        for hit in result.hits:
            for hsp in hit.hsps:
                if useScore: #if we're looking at bitscore
                    if hsp.bitscore>maxscore:
                        hitnames[result.id.split(title_split_char)[0]] = set([hit.id])
                        maxscore = hsp.bitscore
                    elif hsp.bitscore==maxscore:
                        hitnames[result.id.split(title_split_char)[0]].add(hit.id)
                    else:
                        pass
                else: #if we're looking at evalue
                    if hsp.evalue<minevalue or minevalue is False:
                        hitnames[result.id.split(title_split_char)[0]] = set([hit.id])
                        minevalue = hsp.evalue
                    elif hsp.evalue==minevalue:
                        hitnames[result.id.split(title_split_char)[0]].add(hit.id)
                    else:
                        pass
    return(hitnames)


def derep_seq_cli():    
    def parse_args():
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('infile')
        parser.add_argument('outfile')
        parser.add_argument('--approach',default='memory')
        args = parser.parse_args()
        return(args)
    
    args = parse_args()
    
    derep_seq(infile = args.infile, outfile = args.outfile, approach = args.approach)


def derep_seq(infile,outfile,approach='memory',**kwargs):
    """
    Takes a fastq or fasta file, and outputs a deplicated file, with usearch-format annotation
    Where the title is ;-seperated, it treats all text after the first ; as a gene ID string, and therefore
    an identifier string for dereplication (identical sequences with different gene ID strings will not be
    combined during dereplication.
    """
    
    filetype = checkFiletype(infile)
    
    if filetype=='fastq':
        from reptools import FASTQparser as Parser
    else:
        from reptools import FASTAparser as Parser
    
    #read data
    seqs = {}
    
    if approach == 'memory': #store geneid strings in a list, to save memory
        gene_ids = [] # stored in a list to avoid repetition - slower, but saves quite a bit of memory
        with open(infile) as in_handle:
            for seq_tuple in Parser(in_handle):
                title_list = seq_tuple[0].split(';')
                #gene_string = ';'.join(title_list[1:]) #rejoined to save memory - python3 has large per-string overhead
                if title_list[1:] not in gene_ids: #if the geneid is new, the sequence is new (by definition)
                    gene_ids.append(title_list[1:])
                    seqs[(seq_tuple[1],len(gene_ids)-1)] = [1,title_list[0]] #the first title found is stored for output
                elif (seq_tuple[1],gene_ids.index(title_list[1:])) not in seqs: 
                    #if the sequence string is new, but not geneid, the geneid code can be reused
                    seqs[(seq_tuple[1],gene_ids.index(title_list[1:]))] = [1,title_list[0]]
                else: #not a unique sequence
                    seqs[(seq_tuple[1],gene_ids.index(title_list[1:]))][0] +=1
        
        with open(outfile,'w') as out_handle:
            for seq_key in seqs:
                gene_id = gene_ids[seq_key[1]]
                if len(gene_id)>0:
                    title = ';'.join([seqs[seq_key][1],';'.join(gene_id)])
                else:
                    title = seqs[seq_key][1]
                out_handle.write(
                                '>{};size={}\n{}\n'.format(title, str(seqs[seq_key][0]) , seq_key[0]) 
                                )
    
    elif approach == 'speed':
        with open(infile) as in_handle:
            for seq_tuple in Parser(in_handle):
                title_list = seq_tuple[0].split(';')
                gene_string = ';'.join(title_list[1:]) #rejoined to save memory - python3 has large per-string overhead
                if (seq_tuple[1],gene_string) not in seqs: #a new sequence
                    seqs[(seq_tuple[1],gene_string)] = [1,title_list[0]] #the first title found is stored for output
                else: #not a unique sequence
                    seqs[(seq_tuple[1],gene_string)][0] +=1
    
        with open(outfile,'w') as out_handle:
            for seq_key in seqs:
                gene_string = seq_key[1]
                if len(gene_string)>0:
                    title = ';'.join([seqs[seq_key][1],gene_string])
                else:
                    title = seqs[seq_key][1]
                out_handle.write(
                                '>{};size={}\n{}\n'.format(title, str(seqs[seq_key][0]), seq_key[0])
                                )
    else:
        raise ValueError('unknown approach')



