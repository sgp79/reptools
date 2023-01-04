def v081_CDR3slice_dir(
                    indir,
                    outdir,
                    genedict,
                    db_files,
                    db_dir=False,
                    hitsDir=False,
                    store_search_out=False,
                    mincols=(20,15),
                    id=(0.93,0.93),
                    strand='both',
                    evalue=(0.001,0.001),
                    genes=('V','J'),
                    locations = ('VregionC104start','JregionF118start'),
                    filetype='fastq',
                    overwrite=False,
                    usearchpath='usearch',
                    stellarPath='stellar',
                    swipePath='swipe',
                    blastPath = 'blastn',
                    makeblastdbPath = 'makeblastdb',
                    alnout=False, #for debugging use only
                    algorithm='swipe',
                    Vdb_length=30,
                    threads = 10,
                    verbose=True,
                    **kwargs
                    ):
    #modification: the databases and gene dict are now modified to look at the last Vdb_length bases of V only
    import os
    import tempfile
    import reptools
    import reptools.test
    filetypes = reptools.select_filetypes(filetype)
    typefiles = [fn for fn in os.listdir(indir) if os.path.splitext(fn)[1] in filetypes]
    if len(typefiles)==0:
        print('No files of specified type found.\n')
        reptools.ensure_dir(outdir)
        return
    
    if overwrite:
        reptools.remove_dir(outdir)
        if alnout: reptools.remove_dir(alnout)
    reptools.cautious_mkdir(outdir)
    if alnout: reptools.cautious_mkdir(alnout)
    
    if store_search_out is True: 
        #note that store_searchout can be True, False, or contain a path to write uSearch files to
        store_search_out = os.path.join(outdir,'search_reports')
    
    if store_search_out:
        reptools.cautious_mkdir(store_search_out)
    
    if db_dir:
        db_files = [os.path.join(db_dir,fn) for fn in db_files]
    
    counts= {}
    
    if algorithm.lower()=='swipe':
        hitExt='.tsv'
    elif algorithm.lower()=='stellar':
        hitExt='.gff'
    elif algorithm.lower() in ['local','ublast']:
        hitExt='.u14'
    elif algorithm.lower() == 'blast':
        hitExt='.b6'
    else:
        raise ValueError(
            'Unknown algorithm: must be "swipe","stellar","blast","local" or "ublast". '
            'local uses usearch local.'
                        )
    
    #make modified V database and gene dictionary, to work with only last 30 bases (by default - value in Vdb_length)
    (temp_Vdb, temp_genedict) = reptools.make_shorter_db(db_files[0],genedict,Vdb_length)
    db_files = [temp_Vdb,db_files[1]] #new db_files list
    
    #call reptools.CDR3slice()
    for fn in typefiles:
        if hitsDir:
            hitFiles = []
            hitFiles = [fn for fn in os.listdir(hitsDir) if os.path.splitext(fn)[0]==hitExt]
            for gene in genes:
                hitFiles = hitFiles + [
                    os.path.join(hitsDir,fn) for fn in hitFiles if 
                        os.path.splitext(fn)[0][-1]==gene
                    ]
            if len(hitFiles)!= len(genes):
                raise ValueError('Missing hits file for %s' % fn)
        else:
            hitFiles=False
        infile = os.path.join(indir,fn)
        outfile = os.path.join(outdir,fn)
        alnoutfile = False
        if alnout: alnoutfile = os.path.join(alnout,os.path.splitext(fn)[0]+'.aln')
        reptools.test.v081_CDR3slice(
            infile,outfile,db_files=db_files,genedict=temp_genedict,hitFiles=hitFiles,
            store_search_out=store_search_out,mincols=mincols,id=id,strand=strand,
            evalue=evalue,genes=genes,locations=locations,
            usearchpath=usearchpath,stellarPath=stellarPath,swipePath=swipePath,
            blastPath=blastPath,makeblastdbPath=makeblastdbPath,
            alnout=alnoutfile,algorithm=algorithm,threads=threads,verbose=verbose,
            **kwargs
        )
    
    os.remove(temp_Vdb)
    os.remove(temp_genedict)


def v081_CDR3slice(
                infile,
                outfile,
                db_files,
                genedict,
                hitFiles = False,
                store_search_out = False,
                mincols = (20,15),
                id = (0.93,0.93),
                strand = 'both',
                evalue = (0.001,0.001),
                genes = ('V','J'),
                locations = ('VregionC104start','JregionF118start'),
                usearchpath = 'usearch',
                stellarPath = 'stellar',
                swipePath = 'swipe',
                blastPath = 'blastn',
                makeblastdbPath = 'makeblastdb',
                alnout=False,
                algorithm='swipe',
                threads = 10,
                verbose=True,
                **kwargs
                ):
    
    import tempfile
    import os
    import io
    import reptools
    
    genedictionary = reptools.genedictionaryreader(genedict)
    
    ext = os.path.splitext(infile)[1].lower()
    if ext in ['.fastq','.fq']:
        from reptools import FASTQparser as Parser
        def Writer(outhandle,seq_tuple):
            outhandle.write('@%s\n%s\n+\n%s\n' % (seq_tuple[0],seq_tuple[1],seq_tuple[2]))
    else:
        from reptools import FASTAparser as Parser
        def Writer(outhandle,seq_tuple):
            outhandle.write('>%s\n%s\n' % (seq_tuple[0],seq_tuple[1]))
    
    if isinstance(algorithm,str): 
        algorithm = [algorithm for n in range(len(db_files))]
    
    Calls = []
    Exts = []
    for a in algorithm:
        if a.lower() in ['local','ublast']:
            from reptools import uSearchcall as _call
            Exts.append('.u14')
        elif a.lower() == 'stellar':
            from reptools import stellarCall as _call
            Exts.append('.gff')
        elif a.lower() == 'swipe':
            from reptools import swipeCall as _call
            Exts.append('.tsv')
        elif a.lower() == 'blast':
            from reptools import blastCall as _call
            Exts.append('.b6')
        else:
            raise ValueError(
                'Please set algorithm to "swipe" (default), "stellar", "blast" (for blast+ blastn), '
                '"ublast" or "local".\n'
                '"local" uses usearch local search.'
                        )
        Calls.append(_call)
    
    #call usearch/stellar/swipe/blast
    if not hitFiles: #if we need to run a search
        hitFiles=[]
        
        for gene,this_db,this_id,this_mincols,this_alg,this_ext,this_Call in zip(
                                                                      genes,db_files,id,mincols,algorithm,Exts,Calls
                                                                                ):        
            if store_search_out:
                reptools.ensure_dir(store_search_out)
                hits_out=os.path.join(
                    store_search_out,os.path.splitext(os.path.split(infile)[1])[0]+'_'+gene+this_ext
                )
            else:
                temp = tempfile.NamedTemporaryFile(delete=False,suffix=this_ext)
                temp.close()
                hits_out=temp.name
            hitFiles.append(hits_out)
            if alnout:
                alnoutfile = os.path.splitext(alnout)[0]+'_'+gene+os.path.splitext(alnout)[1]
            else:
                alnoutfile = False
            if this_mincols < 16 :
                minhsp = this_mincols
            else:
                minhsp = False
            this_Call(
                        infile=infile,
                        db=this_db,
                        type=this_alg,
                        id=this_id,
                        hitsOut=hits_out,
                        alnout=alnoutfile,
                        strand=strand,
                        top_hits_only=False,
                        mincols=this_mincols,
                        usearchpath = usearchpath,
                        stellarPath = stellarPath,
                        swipePath = swipePath,
                        blastPath = blastPath,
                        makeblastdbPath = makeblastdbPath,
                        minhsp = minhsp,
                        maxrejects = '0',
                        maxaccepts = '0',
                        threads = threads,
                        verbose=verbose,
                        **kwargs
                      )
    
    
    #parse output
    hits=[]
    for fn,gene,location,this_evalue,this_alg in zip(hitFiles,genes,locations,evalue,algorithm):
        hits.append(
                    reptools.find_nearest_to_overlap_hits(
                                                              fn,
                                                              gene,
                                                              genedictionary,
                                                              location,
                                                              filetype=this_alg,
                                                              evalue=this_evalue
                                                          )
                    )
        #this is failing to parse final read in a usearch file #CHECK
    
    if not store_search_out:
        for fn in hitFiles:
            os.remove(fn)
    
    #take slices
    slices = {}
    for title in hits[0]:
        try:
            if hits[0][title].strand==hits[1][title].strand: #if the strand directions don't match, skip this sequence
                this_strand = hits[0][title].strand
                #Jadjust = 4*this_strand #move the end of the slice forward or backward two bases for the full F codon
                if title.split(' ')[0] in slices:
                    raise ValueError(
                                     'Duplicate sequence IDs.  Note that a space acts as an ID terminator during '
                                     'processing, due to SWIPE and BLAST constraints.\n'
                                     )
                slices[title.split(' ')[0]] = tuple([
                                        reptools.get_position(
                                                   hits[0][title],
                                                   pos_H=genedictionary[genes[0]][locations[0]][hits[0][title].hit_id],
                                                   gene='V'
                                                     ),
                                        reptools.get_position(
                                                   hits[1][title],
                                                   pos_H=genedictionary[genes[1]][locations[1]][hits[1][title].hit_id],
                                                   gene='J'
                                                     ),
                                        this_strand
                                       ])
        except KeyError: #if there's no hit for either V or J, we'll get a key error here
            pass
    if len(slices)>0:
        processed = set()
        with io.open(infile) as inhandle,io.open(outfile,'wb') as outhandle:
            for seq_tuple in Parser(inhandle):
                seq_list = [seq_tuple[0]]
                shorttitle = seq_tuple[0].split(' ')[0]
                if shorttitle in processed:
                    raise ValueError(
                                     'Duplicate sequence IDs.  Note that a space acts as an ID terminator during '
                                     'processing, due to SWIPE and BLAST constraints.\n'
                                     )
                else:
                    processed.add(shorttitle)
                if shorttitle in slices:
                    ###make tuple
                    start = slices[shorttitle][0]  -1
                    end = slices[shorttitle][1] + (-2*(slices[shorttitle][2]==-1))  + (2*slices[shorttitle][2])
                    if slices[shorttitle][2] == -1 and end <0: #fix for taking the first character in a reverse slice
                        end = -len(seq_tuple[1])-1
                    #add 2 to slice for forward reads, subtract 4 for reverse reads
                    CDR3 = seq_tuple[1][start:end:slices[shorttitle][2]]  #extract CDR3 seq
                    if slices[shorttitle][2]==-1:
                        CDR3 = reptools.complement(CDR3)
                    if len(CDR3)>0:
                        seq_list.append(CDR3)
                        if len(seq_tuple)==3:
                            seq_list.append(
                                seq_tuple[2][start:end:slices[shorttitle][2]]
                            ) #extract CDR3 qual
                    else:
                        seq_list.append('N')
                        seq_list.append('!')
                else:
                    seq_list.append('N')
                    seq_list.append('!')
                Writer(outhandle,seq_list)                

