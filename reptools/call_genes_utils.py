import os
import tempfile
import shutil
import reptools

def assign_filepairs(fns, pairsuffixes=['_1','_2']):
    """
    Returns a dictionary of filepairs, keyed by file stem
    Throws an error if any files are unpaired
    """
    if len(pairsuffixes[0])!=len(pairsuffixes[1]):
        raise ValueError('The read 1 and read 2 suffixes supplied are not the same length.')
    
    #check that all files have a suffix
    for fn in fns:
        if os.path.splitext(fn)[0][-len(pairsuffixes[0]):] not in pairsuffixes:
            raise ValueError('Input file {} does not have one of the specified filepair suffixes.'.format(fn))
    
    #get a list of filename stems present
    basenames = [os.path.splitext(fn)[0] for fn in fns]
    unique_stems = set([fn[:-len(pairsuffixes[0])] for fn in basenames])
    
    #pair the files
    pairfiles = {}
    for stem in unique_stems:
        pair = [fn for fn in fns if os.path.splitext(fn)[0][:-len(pairsuffixes[0])] == stem]
        if len(pair)==1:
            raise ValueError('No pair file was found for {}'.format(pair[0]))
        elif len(pair)>2:
            raise ValueError('More than two files were found with the same namestem: {}\n'
                             'This may occur if more than one directory is being used as input, but is not\n'
                             'permitted because it would not be possible to differentiate the output.\n'
                             'The recommended workaround is to rename the input files, perhaps by adding\n'
                             'a suffix which identifies their source'.format(', '.join(pair)))
        else:
            pairfiles[stem] = pair
    
    return(pairfiles)



def tiebreak(result1,result2,tiebreaker):
    tiebreaker_values = [
                            getattr(result1[0].hsps[0],tiebreaker[0]),
                            getattr(result2[0].hsps[0],tiebreaker[0])
                        ]
    if tiebreaker_values[0] == tiebreaker_values[1]:
        tiebreaker_choice = 0
    elif tiebreaker[1].lower() == 'lower':
        tiebreaker_choice = (tiebreaker_values.index(min(tiebreaker_values)))+1
    elif tiebreaker[1].lower() == 'higher':
        tiebreaker_choice = (tiebreaker_values.index(max(tiebreaker_values)))+1
    else:
        raise ValueError(
    'tiebreaker must be a list or tuple where the 1st element is the property to use, '
    'and the 2nd element is either "higher" or "lower".'
    )
    return(tiebreaker_choice)


def VsegSlicer(hit,genedict):
    direction = hit.strand
    slicepoint = reptools.get_position(
                           hit,
                           pos_H=genedict['VregionC104start'][hit.hit_id],
                           gene='V'
                           )
    ###make tuple
    if direction == 1:
        start = 0
        end = slicepoint + (-2*(direction==-1))  + (2*direction)
    else:
        start = slicepoint + (-2*(direction==-1))  + (2*direction)
        end = -1
    return(start,end,direction)


def get_position(hsp,pos_H,hitAnchor=False,gene=False):
    """
    hitAnchor should be set to start or end, to use hit_start or hit_end as the reference point
    if hitAnchor is False, but gene is supplied, the default is to use hit_end for 'V' and hit_start for 'J'
    """
    if gene and hitAnchor:
        raise ValueError('Supply only ONE of hitAnchor and gene.')
    elif gene:
        if gene.upper()=='V':
            hitAnchor = 'end'
        elif gene.upper()=='J':
            hitAnchor = 'start'
        else:
            raise ValueError('If gene is not False, it must be "V" or "J".')
    else:
        pass
    if hitAnchor == 'start':
        if hsp.strand == 1:
            align_Q,align_H = hsp.query_start,hsp.hit_start
        elif hsp.strand == -1:
            align_Q,align_H = hsp.query_end,hsp.hit_end
        else:
            raise ValueError('hsp object must have strand set to 1 or -1')
    elif hitAnchor == 'end':
        if hsp.strand == 1:
            align_Q,align_H = hsp.query_end,hsp.hit_end
        elif hsp.strand == -1:
            align_Q,align_H = hsp.query_start,hsp.hit_start
        else:
            raise ValueError('hsp object must have strand set to 1 or -1')
    else:
        raise ValueError('hitAnchor must be "start" or "end"')
    return(align_Q - hsp.strand*(align_H - pos_H))


def filterHits(hitList,mincols=False,id=False,evalue=False):
    if mincols:
        hitList = [hit for hit in hitList if max([int(hsp.aln_span) for hsp in hit.hsps])>=mincols] #filter by mincols
    if id:
        hitList = [hit for hit in hitList if max([float(hsp.ident_fract) for hsp in hit.hsps])>=id] #filter by id
    if evalue:
        hitList = [hit for hit in hitList if min([float(hsp.evalue) for hsp in hit.hsps])<=evalue] #filter by evalue
    return(hitList)


def parseHitfile(fn,aligner,title_split=' ',clusterID_position=0,mincols=False,id=False,evalue=False):
    results = {}
    for result in reptools.hitsParser(fn,aligner):
        results[result.id.split(title_split)[clusterID_position]] = reptools.filterHits(
                                                                                        result.hits,
                                                                                        mincols=mincols,
                                                                                        id=id,
                                                                                        evalue=evalue
                                                                                        )
    results = {k:results[k] for k in results if len(results[k])>0} #remove queries left empty by filter
    return(results)


def get_top_hits(
                 hitsfile,filetype='blast',
                 title_split=None,clusterID_position=0,criteria=['bitscore','high'],
                 mincols=False,id=False,evalue=False
                 ):
    """
    Filetype should be "stellar","usearch", "ublast", "local", global", "blast", "swipe", or "dict"
    Returns a dictionary keyed by query id, containing the best hit or hits for each query
    """
    besthits = {}
    if criteria[1].lower().startswith('high'):
        selectionoperator = max
    elif criteria[1].lower().startswith('low'):
        selectionoperator = min
    else:
        raise ValueError(
            '"criteria" must be a list of len=2, with the first item giving the name of the property to select on, and '
            'the second being either "low" or "high"'
        )
    if filetype == 'dict':
        for result in hitsfile:
            besthits[result.split(title_split)[clusterID_position]] = [
                                #this nested list comp. selects hits if their max (or min) hsp value is equal to the 
                                #max/min value of any hit
                                hit for hit in hitsfile[result]
                                  if
                                    selectionoperator([getattr(hsp,criteria[0]) for hsp in hit.hsps])
                                    ==
                                    selectionoperator([
                                                      selectionoperator([getattr(hsp,criteria[0]) for hsp in hit2.hsps])
                                                       for hit2 in hitsfile[result]
                                                      ])
                                    ]
    else:
        for result in reptools.hitsParser(hitsfile,filetype=filetype):
            title = result.id.split(title_split)[clusterID_position]
            besthits[title] = [
                                #this nested list comp. selects hits if their max (or min) hsp value is equal to the 
                                #max/min value of any hit
                                hit for hit in result.hits 
                                  if 
                                    selectionoperator([getattr(hsp,criteria[0]) for hsp in hit.hsps])
                                    ==
                                    selectionoperator([
                                                      selectionoperator([getattr(hsp,criteria[0]) for hsp in hit2.hsps])
                                                       for hit2 in result.hits
                                                      ])
                               ]
    
    besthits = {k:reptools.filterHits(besthits[k],mincols=mincols,id=id,evalue=evalue) for k in besthits}
    besthits = {k:besthits[k] for k in besthits if len(besthits[k])>0} #remove queries left empty by filter
    return(besthits)


def find_nearest_to_overlap_hits(
                                 hitsfile,gene,genedictionary,dict_entry,
                                 filetype='blast',evalue=False,title_split=None,clusterID_position=0,
                                 criteria=('bitscore','high'),
                                 top_hits_only=False,returnhitstrings=False,
                                 mincols=False, id=False
                                 ):
    """
    Filetype should be "stellar","usearch", "ublast", "local", global", "blast", or "swipe"
    """
    if returnhitstrings and not top_hits_only:
        raise ValueError('returnhitstrings=True must be used with top_hits_only=True')
    if criteria[1].lower().startswith('high'):
        #selectionoperator = operator.gt
        selectionoperator = max
    elif criteria[1].lower().startswith('low'):
        #selectionoperator = operator.lt
        selectionoperator = min
    else:
        raise ValueError(
            '"criteria" must be a list of len=2, with the first item giving the name of the property to select on, and '
            'the second being either "low" or "high"'
        )
    besthits = {}
    hitstrings = {}
    for result in reptools.hitsParser(hitsfile,filetype=filetype,evalue=evalue):
        closesthsp={}
        #loop through the hits, and from each, select the hsp closest to the junction
        if len(result.hits) == 0: continue
        hits = result.hits
        if top_hits_only: #first, select only the top hits, if instructed
            hits = [
                    h for h in hits if 
                     selectionoperator([getattr(hsp,criteria[0]) for hsp in h.hsps])
                     ==
                     selectionoperator([
                                        selectionoperator([getattr(hsp,criteria[0]) for hsp in h2.hsps])
                                         for h2 in hits
                                       ])
                   ]
        hits = reptools.filterHits(hits, mincols=mincols, id=id) #filter hits by mincols and id
        if len(hits) == 0 : continue
        if returnhitstrings:
            hitstrings[result.id.split(title_split)[clusterID_position]] = '+'.join([h.id for h in hits])
        for hit in hits:
            pos_H = genedictionary[gene][dict_entry][hit.id]
            distancefromjunction=[]
            for hsp in hit.hsps:
                #does hsp surround hit?
                if (hsp.hit_start<pos_H and hsp.hit_end>pos_H) or (hsp.hit_start>pos_H and hsp.hit_end<pos_H): 
                    distancefromjunction.append(0)
                else: #hsp doesn't surround hit
                    distancefromjunction.append(min(abs(hsp.hit_start-pos_H),abs(hsp.hit_end-pos_H)))
            #select the hsps with the minimal distance from the junction
            closest_dist = min(distancefromjunction)
            candidatehsps = [(hit.hsps[n],d) for n,d in enumerate(distancefromjunction) if d==closest_dist]
            scores = [getattr(hsptuple[0],criteria[0]) for hsptuple in candidatehsps]
            closesthsp[hit.id] = candidatehsps[scores.index(selectionoperator(scores))]
        #now select the best overall hit: ties are broken
        #need to pick closest hit by distance, then the best of the closest by selection criteria
        dists = [closesthsp[hit][1] for hit in closesthsp]
        closesthits = [closesthsp[hit][0] for hit in closesthsp if closesthsp[hit][1]==min(dists)]
        scores = [getattr(hsp,criteria[0]) for hsp in closesthits]
        besthits[result.id.split(title_split)[clusterID_position]]=closesthits[scores.index(selectionoperator(scores))]
    
    if returnhitstrings:
        return(tuple([besthits,hitstrings]))
    else:
        return(besthits)


def getForwardreadpos(targetPos,queryStart,hitStart,hitStrand):
    """
    This takes a nucleotide position which relates to the database hit
    (i.e the V region gene from IMGT) and identifies its position in the query sequence.
    It requires a target position (i.e. the first base of start codon, where the first base=1) and
    queryStart and hitStart data from the SearchIO object.
    It returns the Pythonesque position in the query.
    """
    if hitStrand == 1:
        offset =  hitStart - queryStart
    elif hitStrand == -1:
        offset =  - (hitStart + queryStart)
    else:
        raise ValueError('hitStrand must be -1 or 1')
    queryPos = targetPos + offset
    return(queryPos)


def genedictionaryreader(fn):
    """
    Need to put csv layout info here
    """
    import csv
    #
    VregionID = {"None":0,"Multiple":-1} #an ID number for each sequence in the BLAST DB
    Vgenes = dict()
    #a dictionary containing each V gene, with each BLAST DB sequence associated with it (to allow for splice variants, etc.)
    VgeneID = {"None":0,"Multiple":-1}
    Vfamilies = dict()
    VfamilyID ={"None":0,"Multiple":-1}
    Vseqs = set()
    pseudoVgenes = set()
    VregionC104start = dict()
    startcodon = dict()
    Vstart=dict()
    #
    JregionID = {"None":0,"Multiple":-1} #an ID number for each sequence in the BLAST DB
    JgeneID = {"None":0,"Multiple":-1} #an ID number for each sequence in the BLAST DB
    Jgenes = dict()
    #a dictionary containing each V gene, with each BLAST DB sequence associated with it (to allow for splice variants, etc.)
    Jseqs = set()
    pseudoJgenes = set()
    JregionF118start = dict()
    #
    CregionID = {"None":0,"Multiple":-1}
    Cgenes = dict()
    #a dictionary containing each V gene, with each BLAST DB sequence associated with it (to allow for splice variants, etc.)
    Cseqs = set()
    pseudoCgenes = set()
    #
    with open(fn) as dictfile: 
        for ln,row in enumerate(csv.reader(dictfile)):
            if ln==0:
                if 'Sequence' in row:
                    seqcol = row.index('Sequence')
                else:
                    raise IOError('Specified dictionary file does not have a "Sequence" column')
                if 'Gene' in row:
                    genecol = row.index('Gene')
                else:
                    raise IOError('Specified dictionary file does not have a "Gene" column')
                if 'Family' in row:
                    familycol = row.index('Family')
                else:
                    raise IOError('Specified dictionary file does not have a "Family" column')
                if 'Type' in row:
                    typecol = row.index('Type')
                else:
                    raise IOError('Specified dictionary file does not have a "Type" column')
                if 'Pseudo' in row:
                    pseudocol = row.index('Pseudo')
                else:
                    pseudocol = False
                    pseudoVgenes = False
                if 'C104' in row:
                    C104col = row.index('C104')
                else:
                    C104col = False
                    VregionC104start = False
                if 'F118' in row:
                    F118col = row.index('F118')
                else:
                    F118col = False
                    JregionF118start = False
                if 'ATG' in row:
                    ATGcol = row.index('ATG')
                else:
                    ATGcol = False
                    startcodon = False
                if 'Vstart' in row:
                    Vstartcol = row.index('ATG')
                else:
                    Vstartcol = False
                    Vstart = False
            else: #all rows apart from first
                seqname = row[seqcol].strip() #strip any leading or trailing whitespace, and use this as the dictionary key in subsequent lines
                genename = row[genecol].strip() #strip any leading or trailing whitespace, and use this as the dictionary key in subsequent lines
                familyname = row[familycol].strip() #strip any leading or trailing whitespace, and use this as the dictionary key in subsequent lines
                if row[typecol]=='V':
                    VregionID[seqname] = len(VregionID)-1 #assign an ID number
                    Vseqs.add(seqname)
                    if genename in list(Vgenes.keys()):
                        Vgenes[genename].add(seqname)
                    else:
                        VgeneID[genename] = len(VgeneID)-1 #assign an ID number
                        Vgenes[genename] = set()
                        Vgenes[genename].add(seqname)
                    if familyname in list(Vfamilies.keys()):
                        Vfamilies[familyname].add(genename)
                    else:
                        VfamilyID[familyname] =len(VfamilyID)-1
                        Vfamilies[familyname] = set()
                        Vfamilies[familyname].add(genename)
                    if pseudocol:
                        if int(row[pseudocol])==1:
                            pseudoVgenes.add(genename)                
                    if C104col: VregionC104start[seqname] = int(row[C104col])
                    if ATGcol:
                        try:
                            startcodon[seqname] = int(row[ATGcol])
                        except ValueError:
                            startcodon[seqname] = int(0)
                    if Vstartcol:
                        try:
                            Vstart[seqname] = int(Vstartcol)
                        except ValueError:
                            Vstart[seqname] = int(0)
                #
                elif row[typecol]=='J':
                    JregionID[seqname] = len(JregionID)-1
                    Jseqs.add(seqname)
                    if genename in list(Jgenes.keys()):
                        Jgenes[genename].add(seqname)
                    else:
                        JgeneID[genename] = len(JgeneID)-1
                        Jgenes[genename] = set()
                        Jgenes[genename].add(seqname)
                    if F118col: JregionF118start[seqname] = int(row[F118col])
                    if pseudocol:
                        if int(row[pseudocol])==1:
                            pseudoJgenes.add(genename)
                #
                elif row[typecol]=='C':
                    CregionID[seqname] = len(CregionID)-1
                    Cseqs.add(seqname)
                    if genename in list(Cgenes.keys()):
                        Cgenes[genename].add(seqname)
                    else:
                        Cgenes[genename] = set()
                        Cgenes[genename].add(seqname)
                    if pseudocol:
                        if int(row[pseudocol])==1:
                            pseudoCgenes.add(genename)
                else:
                    pass
    #        
    lookup_VregionID = dict([[v,k] for k,v in list(VregionID.items())])
    lookup_Vgenes = {seq:g for g in list(Vgenes.keys()) for seq in Vgenes[g]}
    lookup_VgeneID = dict([[v,k] for k,v in list(VgeneID.items())])
    lookup_Vfamilies = {gene:f for f in list(Vfamilies.keys()) for gene in Vfamilies[f]}
    lookup_VfamilyID = dict([[v,k] for k,v in list(VfamilyID.items())])
    #
    lookup_JregionID = dict([[v,k] for k,v in list(JregionID.items())])
    lookup_Jgenes = {seq:g for g in list(Jgenes.keys()) for seq in Jgenes[g]}
    lookup_JgeneID = dict([[v,k] for k,v in list(JgeneID.items())])
    #
    lookup_CregionID = dict([[v,k] for k,v in list(CregionID.items())])
    #lookup_Cgenes = {seq:g for g in Cgenes.keys() for seq in Cgenes[g]}
    #
    return({
            'C':{'CregionDict':CregionID,'Cseqs':Cseqs,'Cgenes':Cgenes,'pseudoCgenes':pseudoCgenes,'lookup_CregionDict':lookup_CregionID},
            'V':{'VregionID':VregionID,'Vseqs':Vseqs,'Vgenes':Vgenes,'VgeneID':VgeneID,'Vfamilies':Vfamilies,'VfamilyID':VfamilyID,'pseudoVgenes':pseudoVgenes,'lookup_VregionID':lookup_VregionID,'lookup_Vgenes':lookup_Vgenes,'lookup_VgeneID':lookup_VgeneID,'lookup_VfamilyID':lookup_VfamilyID,'lookup_Vfamilies':lookup_Vfamilies,'VregionC104start':VregionC104start,'startcodon':startcodon,'Vstart':Vstart},
            'J':{'JregionID':JregionID,'Jseqs':Jseqs,'Jgenes':Jgenes,'JgeneID':JgeneID,'pseudoJgenes':pseudoJgenes,'lookup_JregionID':lookup_JregionID,'lookup_Jgenes':lookup_Jgenes,'lookup_JgeneID':lookup_JgeneID,'JregionF118start':JregionF118start}
            })


def make_shorter_db(Vdb,genedict,length):
    """
    Function to trim a V gene fasta file to leave only the 3' end, and to produce a gene dictionary which has been
    modified accordingly.
    For use by reptools.CDR3slice_dir()
    Takes filenames for the original V gene fasta dictionary, the sequences in which should terminate at the conserved
    cysteine residue codon; the gene dictionary csv filename; the number of bases to keep in the output.
    Returns the filenames for the modified V gene fasta dictionary and gene dictionary.
    """
    Vlengths = {}
    newVdb = reptools.createTempfile()
    with open(Vdb) as old_db, open(newVdb,'w') as new_db:
        for title,seq in reptools.FASTAparser(old_db):
            new_db.write('>%s\n%s\n' % (title,seq[-length:]))
            Vlengths[title] = len(seq)
    
    genedictionary = genedictionaryreader(genedict)
    temp_genedict = reptools.createTempfile()
    with open(temp_genedict,'w') as new_dict:
        #minimal gene dictionary structure required here
        new_dict.write('Sequence,Gene,Family,Type,C104,F118\n')
        for Vseq in genedictionary['V']['Vseqs']:
            new_dict.write(
                '%s,,,V,%s,\n' % (
                    Vseq,
                    #calculate position in new fasta
                    (genedictionary['V']['VregionC104start'][Vseq]-Vlengths[Vseq]+length)
                                  )
                          )
        for Jseq in genedictionary['J']['Jseqs']:
            new_dict.write(
                '%s,,,J,,%s\n' % (Jseq,genedictionary['J']['JregionF118start'][Jseq])
            )
    
    return(newVdb,temp_genedict)


def setCaller(aligner):
    a = aligner.lower()
    if a in ['usearch','ublast']:
        from reptools import uSearchcall as _call
        _ext = '.u14'
    elif a == 'vsearch':
        from reptools import vsearchCall as _call
        _ext = '.u14'
    elif a == 'stellar':
        from reptools import stellarCall as _call
        _ext = '.gff'
    elif a in ['blastn','blast']:
        from reptools import blastCall as _call
        _ext = '.b6'
    elif a == 'swipe':
        _ext = '.swipe'
        from reptools import swipeCall as _call
    elif a == 'ssw':
        _ext = '.ssw'
        from reptools import sswCall as _call
    else:
        raise ValueError('unknown aligner')
    return(_call,_ext)


def getReadIDs(in1,in2,title_split=' ',clusterID_position=0):
    readIDs = []
    with open(in1) as inhandle1, open(in2) as inhandle2:
        for (title1,seq1,qual1),(title2,seq2,qual2) in zip(
                                                           reptools.FASTQparser(inhandle1),
                                                           reptools.FASTQparser(inhandle2)
                                                           ): #loop through the fastq lines
            if title1.split(title_split)[clusterID_position] != title2.split(title_split)[clusterID_position]:
                raise IOError('Sequence titles do not match between files:\n{}\n{}\n'.format(title1,title2))
            else:
                readIDs.append(title1.split(title_split)[clusterID_position])
    return(readIDs)


def initialchecks(infiles,aligner=None):
    import os
    
    if aligner is not None: #None skips this check, allowing initialchecks() to be used to get filetype info only
        if aligner.lower() not in ['usearch','vsearch','ublast','blastn','swipe','stellar','ssw']:
            raise ValueError('Unknown aligner')
    
    exts = [os.path.splitext(fn)[1].lower() for fn in infiles]
    
    if len(set(exts))>1:
        raise IOError(
         'infiles do not appear to be of same type (different filename extensions)'
        )
    
    if exts[0] in ['.fastq','.fq']:
        from reptools import FASTQparser as Parser
        def Writer(outhandle,seq_tuple):
            outhandle.write('@{}\n{}\n+\n{}\n'.format(seq_tuple[0],seq_tuple[1],seq_tuple[2]))
    else:
        raise IOError(
         "fastq handling; fasta to be implemented later (but won't allow for EE filtering)"
        )
    
    return(Writer,Parser,exts[0])


