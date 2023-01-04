def descr_parser(infile):
    """
    Takes a description file (.descr) produced by miXCR, and parses it, outputting a dictionary of form
    sequenceID:['V'],['J'],['CDR3']
    """
    gene_dict = {}
    import csv
    for row in csv.reader(open(infile),delimiter='\t'):
        if row[0] == 'Index of read':
            pass
        else:
            gene_dict[row[0]]={}
            gene_dict[row[0]]['V']=row[1].split('*')[0]
            gene_dict[row[0]]['J']=row[2].split('*')[0]
            gene_dict[row[0]]['CDR3']=row[3]    
    return(gene_dict)


def countgenehits_hitsfile(
                    fastqFile,
                    hitsFile,
                    type,
                    evalue = False,
                    mincols = False,
                    title_split=' ',
                    verbose = True
                    ):
    import reptools
    import io
    results = {'hit':0, 'ambiguous':0, 'fail':0}
    reptoolsdict = reptools.retrieve_tophits(
                                              hitsFile,
                                              type,
                                              evalue=evalue,
                                              mincols=mincols,
                                              title_split_char=title_split
                                            )
    
    with io.open(fastqFile) as fastq_handle:
        for title,seq,qual in reptools.FASTQparser(fastq_handle):
            trimmed_id = title.split(title_split)[0]
            if trimmed_id in reptoolsdict:
                if len(reptoolsdict[trimmed_id]) == 1:
                    results['hit']+=1
                elif len(reptoolsdict[trimmed_id]) > 1:
                    results['ambiguous']+=1
                else:
                    results['fail']+=1
            else:
                results['fail']+=1
                
    if verbose:
        print(('hit = %s' % results['hit']))
        print(('fail = %s' % results['fail']))
        print(('ambiguous = %s' % results['ambiguous']))
    
    return(results)


def checkgeneID_hitsfile(
                    hitsFile,
                    mock_dict,
                    type,
                    gene,
                    evalue = False,
                    mincols = False,
                    title_split=' ',
                    verbose = True
                    ):
    
    import reptools

    reptoolsdict = reptools.retrieve_tophits(
                                              hitsFile,
                                              type,
                                              evalue=evalue,
                                              mincols=mincols,
                                              title_split_char=title_split
                                            )
    
    reptoolsdict = {k:[s.split('gamma')[0].split('alpha')[0].split('_')[0] for s in reptoolsdict[k]] for k in reptoolsdict }
    
    results = {'hit':0, 'miss':0, 'ambiguous':0, 'fail':0}
    
    for id in mock_dict:
        trimmed_id = id.split(title_split)[0]
        try:
            if len(reptoolsdict[trimmed_id]) == 1:
                if [h for h in reptoolsdict[trimmed_id]][0] == mock_dict[id][gene].split('*')[0]:
                    results['hit']+=1
                elif reptoolsdict[trimmed_id]=='none':
                    results['fail']+=1
                else:
                    results['miss']+=1
            elif len(reptoolsdict[trimmed_id]) == 0:
                results['fail']+=1
            elif len(reptoolsdict[trimmed_id]) > 1:
                if mock_dict[id][gene].split('*')[0] in reptoolsdict[trimmed_id]:
                    results['ambiguous']+=1
                else:
                    results['miss']+=1
            else:
                raise ValueError('How did we get here?')
        except KeyError:
            results['fail']+=1
    
    if verbose:
        print(('hit = %s' % results['hit']))
        print(('miss = %s' % results['miss']))
        print(('ambiguous (including hit) = %s' % results['ambiguous']))
        print(('failed = %s' % results['fail']))
        print(('pct hit = %s' % (results['hit']/float(sum(results.values()))*100)))
        print(('pct hit (including ambiguous hit) = %s' % ((results['hit']+results['ambiguous'])/float(sum(results.values()))*100)))
    
    return(results)
                    
                    

    

def checkgeneID_fastq(
                    mock_dict,gene,fastqfile,title_split=' ',
                    hits_out=False,misses_out=False,failures_out=False,ambiguous_out=False,
                    verbose = True
                    ):
    #TODO: add transtable option to replace the [s.split('gamma')[0].split('alpha')[0].split('_')[0] for s in gene_strings] line
    #or, process the mock_dict first, to match
    import reptools
    reptoolsdict = {}
    with open(fastqfile) as infile:
        for title,seq,qual in reptools.FASTQparser(infile):
            id = title.split(';')[0].split(title_split)[0]
            gene_strings = [s.split('=')[1] for s in title.split(';') if s.split('=')[0]==gene]
            gene_strings = gene_strings[0].split('+')
            gene_strings = [s.split('gamma')[0].split('alpha')[0].split('_')[0] for s in gene_strings]
            reptoolsdict[id] = gene_strings
    
    results = {'hit':0, 'miss':0, 'ambiguous':0, 'fail':0}
    
    with open(hits_out,'wb') if hits_out else reptools.dummy_context_mgr() as hits_handle:
        with open(misses_out,'wb') if misses_out else reptools.dummy_context_mgr() as misses_handle:
            with open(failures_out,'wb') if failures_out else reptools.dummy_context_mgr() as failures_handle:
                with open(ambiguous_out,'wb') if ambiguous_out else reptools.dummy_context_mgr() as ambiguous_handle:
                    for id in mock_dict:
                        try:
                            if len(reptoolsdict[id]) == 1:
                                if reptoolsdict[id][0] == mock_dict[id][gene].split('*')[0]:
                                    results['hit']+=1
                                    if hits_out:
                                        hits_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                                elif reptoolsdict[id]=='none':
                                    results['fail']+=1
                                    failures_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                                else:
                                    results['miss']+=1
                                    misses_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                            elif len(reptoolsdict[id]) == 0:
                                results['fail']+=1
                                failures_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                            elif len(reptoolsdict[id]) > 1:
                                if mock_dict[id][gene].split('*')[0] in reptoolsdict[id]:
                                    results['ambiguous']+=1
                                    ambiguous_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                                else:
                                    results['miss']+=1
                                    misses_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                            else:
                                raise ValueError('How did we get here?')
                        except KeyError:
                            results['fail']+=1
                            failures_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
    if verbose:
        print(('hit = %s' % results['hit']))
        print(('miss = %s' % results['miss']))
        print(('ambiguous (including hit) = %s' % results['ambiguous']))
        print(('failed = %s' % results['fail']))
        print(('pct hit = %s' % (results['hit']/float(sum(results.values()))*100)))
        print(('pct hit (including ambiguous hit) = %s' % ((results['hit']+results['ambiguous'])/float(sum(results.values()))*100)))
    
    return(results)


def checkCDR3_fastq(
                    mock_dict,fastqfile,
                    hits_out=False,misses_out=False,failures_out=False,
                    title_split=' ',verbose = True
                    ):
    """
    Compares the sliced CDR3 in a fastq file with those simulated by MiXCR (and given in a .descr file)
    """
    import reptools
    reptoolsdict = {}
    with open(fastqfile) as infile:
        for title,seq,qual in reptools.FASTQparser(infile):
            id = title.split(';')[0].split(title_split)[0]
            reptoolsdict[id] = seq
    
    results = {'hit':0, 'miss':0, 'fail':0}
    
    with open(hits_out,'wb') if hits_out else reptools.dummy_context_mgr() as hits_handle:
        with open(misses_out,'wb') if misses_out else reptools.dummy_context_mgr() as misses_handle:
            with open(failures_out,'wb') if failures_out else reptools.dummy_context_mgr() as failures_handle:
                for id in mock_dict:
                    try:
                        if reptoolsdict[id].lower() == mock_dict[id]['CDR3'].lower():
                            results['hit']+=1
                            hits_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                        else:
                            if reptoolsdict[id].lower()=='n':
                                results['fail']+=1
                                failures_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                            else:
                                results['miss']+=1
                                misses_handle.write('>%s\n%s\n' % (id,reptoolsdict[id]))
                    except KeyError:
                        results['fail']+=1
                        print(id)
                        failures_handle.write('>%s\n%s\n' % (id,''))
                        
    if verbose:
        print(('hit = %s' % results['hit']))
        print(('miss = %s' % results['miss']))
        print(('failed = %s' % results['fail']))
        print(('pct hit = %s' % (results['hit']/float(sum(results.values()))*100)))
        print(('pct miss = %s' % (results['miss']/float(sum(results.values()))*100)))
    return(results)

    
def checkCDR3_prod(
                    fastqfile,
                    minlen = 3*5, maxlen = 3*30,
                    startchars='C',endchars='FWH',
                    hits_out=False,failures_out=False,frameshift_out=False,long_out=False,short_out=False,
                    stop_out=False,bad_out=False,
                    title_split=' ',verbose = True
                    ):
    """
    This for use where no reference file is available.
    Reports % of CDR3 which are productive or start with C and end with F/W/H, and are within a sensible length range
    Over- and under-length CDR3 are eliminated first
    Then those with a bad start or end residue (not C and F/W/H)
    Then those with a stop
    Then those with a frameshift
    """
    import reptools
    results = {'good':0, 'frameshift':0, 'stop':0, 'bad':0, 'long':0, 'short':0, 'fail':0}

    with open(fastqfile) as infile:
        with open(hits_out,'wb') if hits_out else reptools.dummy_context_mgr() as hits_handle:
            with open(failures_out,'wb') if failures_out else reptools.dummy_context_mgr() as failures_handle:
                with open(frameshift_out,'wb') if frameshift_out else reptools.dummy_context_mgr() as shift_handle:
                    with open(long_out,'wb') if long_out else reptools.dummy_context_mgr() as long_handle:
                        with open(short_out,'wb') if short_out else reptools.dummy_context_mgr() as short_handle:
                            with open(stop_out,'wb') if stop_out else reptools.dummy_context_mgr() as stop_handle:
                                with open(bad_out,'wb') if bad_out else reptools.dummy_context_mgr() as bad_handle:
                                    for title,seq,qual in reptools.FASTQparser(infile):
                                        id = title.split(';')[0].split(title_split)[0]
                                        seq = seq.strip()
                                        if seq.lower()=='n':
                                            results['fail']+=1
                                            failures_handle.write('>%s\n%s\n' % (id,seq))
                                        elif len(seq)>maxlen:
                                            results['long']+=1
                                            long_handle.write('>%s\n%s\n' % (id,seq))
                                        elif len(seq)<minlen:
                                            results['short']+=1
                                            short_handle.write('>%s\n%s\n' % (id,seq))
                                        elif (
                                                reptools.trans(seq[0:3]).lower() not in startchars.lower()
                                                or 
                                                reptools.trans(seq[-3:]).lower() not in endchars.lower()
                                        ):
                                            results['bad']+=1
                                            bad_handle.write('>%s\n%s\n' % (id,seq))
                                        elif '*' in reptools.trans(seq):
                                            results['stop']+=1
                                            stop_handle.write('>%s\n%s\n' % (id,seq))
                                        elif len(seq)%3!=0:
                                            results['frameshift']+=1
                                            shift_handle.write('>%s\n%s\n' % (id,seq))
                                        else:
                                            results['good']+=1
                                            hits_handle.write('>%s\n%s\n' % (id,seq))
    
    if verbose:
        totalreads = float(sum(results.values()))
        print(('over length = %s (%s pct)' % (results['long'], results['long']/totalreads*100)))
        print(('under length = %s (%s pct)' % (results['short'], results['short']/totalreads*100)))
        print(('bad start/end = %s (%s pct)' % (results['bad'], results['bad']/totalreads*100)))
        print(('stop codon = %s (%s pct)' % (results['stop'], results['stop']/totalreads*100)))
        print(('frameshift = %s (%s pct)' % (results['frameshift'], results['frameshift']/totalreads*100)))
        print(('no CDR3 = %s (%s pct)' % (results['fail'], results['fail']/totalreads*100)))
        print(('good CDR3 = %s (%s pct)' % (results['good'], results['good']/totalreads*100)))
    return(results)
