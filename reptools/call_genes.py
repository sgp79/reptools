import os
import tempfile
import shutil
import reptools
import contextlib
import collections
from reptools import fill_defaults, select_filetypes, assign_filepairs, checkFiletype
from reptools import dummy_context_mgr as dummy
import errno

def call_pairslist(
            filepairs,
            dbs,
            genedictfile,
            outdir,
            db_dir=False,
            noCDR3=False,
            notrim=False,
            adaptertrimmed_dir=None,
            adapters=False,
            hits_dir=None,
            CDR3_dir=None,
            Cmiss_dir=False,
            ChitVmiss_dir=False,
            ChitVhitJmiss_dir=False,
            Vsegmentout_dir=False,
            filetype='fastq',
            pairsuffixes=('_1','_2'),
            title_split=' ',
            overwrite=False,
            labels=False,
            mincols=False,
            id=False,
            evalue=False,
            wordlength=False,
            gapopen=False,
            gapextend=False,
            aligners=False,
            aligner_paths=False,
            threads=4,
            Vdb_length=30,
            tiebreaker = ['evalue','lower'],
            clusterID_position = 0,
            blastdb_version=5
            #reward=1, #commented out options are not currently implemented, but could be sent to aligner
            #penalty=-3,
            #targetcov=None,
            #minhsp=False,
            #xDrop=3,
            #alphabet='dna5',
            #numMatches=10,
            #verbose=True,
    ):
    tempdirs = []
    
    #########
    #haven't I done the block below in args?
    #adaptertrimmed_dir,td = reptools.build_path(not notrim, adaptertrimmed_dir, 'adaptertrimmed', outdir)
    #tempdirs.append(td)
    #hits_dir,td = reptools.build_path(True, hits_dir, 'fastq_hits', outdir)
    #tempdirs.append(td)
    #CDR3_dir,td = reptools.build_path(not noCDR3, CDR3_dir, 'rawCDR3', outdir)
    #tempdirs.append(td)
    #
    #if not notrim and not adapters:
    #    raise ValueError('Unless --notrim is set, adapters must be supplied.')
    ##########
    
    #make output directories, deleting pre-existing data if overwrite is set
    for pth in [
                adaptertrimmed_dir,hits_dir,CDR3_dir,
                Cmiss_dir,ChitVmiss_dir,ChitVhitJmiss_dir,Vsegmentout_dir
                ]:
        if pth: reptools.cautious_mkdir(pth,overwrite=overwrite)
    
    if db_dir:
        dbs = {gene:os.path.join(db_dir,dbs[gene]) for gene in dbs}            
    
    #call
    for stem,filepair in filepairs.items():
        in1=filepair[0]
        in2=filepair[1]
        adaptertrimmed1,adaptertrimmed2 = reptools.make_paired_filepaths(adaptertrimmed_dir, stem, pairsuffixes)
        outF,outR = reptools.make_paired_filepaths(hits_dir, stem, ['_F','_R'])
        Cmiss1,Cmiss2 = reptools.make_paired_filepaths(Cmiss_dir, stem, pairsuffixes)
        VmissF, VmissR = reptools.make_paired_filepaths(ChitVmiss_dir, stem, ['_F','_R'])
        JmissF, JmissR = reptools.make_paired_filepaths(ChitVhitJmiss_dir, stem, ['_F','_R'])
        outVsegF, outVsegR = reptools.make_paired_filepaths(Vsegmentout_dir, stem, ['_F','_R'])
        
        outCDR3 = reptools.make_unpaired_filepaths(CDR3_dir, stem)
        
        _ = reptools.call_filepair(
                    in1=in1, in2=in2,
                    dbs=dbs,
                    genedictfile=genedictfile,
                    noCDR3=noCDR3,
                    notrim=notrim,
                    adaptertrimmed1=adaptertrimmed1, adaptertrimmed2=adaptertrimmed2,
                    outF=outF, outR=outR,
                    outCDR3=outCDR3,
                    Cmiss1=Cmiss1, Cmiss2=Cmiss2,
                    JmissF=JmissF, JmissR=JmissR,
                    VmissF=VmissF, VmissR=VmissR,
                    outVsegF=outVsegF, outVsegR=outVsegR,
                    adapters=adapters,
                    labels=labels,
                    aligners=aligners,
                    aligner_paths=aligner_paths,
                    wordlength=wordlength,
                    mincols=mincols,
                    id=id,
                    evalue=evalue,
                    gapopen=gapopen,
                    gapextend=gapextend,
                    title_split=title_split,
                    tiebreaker = tiebreaker,
                    Vdb_length=Vdb_length,
                    threads=threads,
                    overwrite=True,
                    clusterID_position=clusterID_position,
                    blastdb_version=blastdb_version
                    )
    reptools.clean_up_dirs(tempdirs)
    
    return(outCDR3)


def call_filepair(
         in1,
         in2,
         dbs,
         genedictfile,
         noCDR3=False,
         notrim=False,
         adaptertrimmed1=False, adaptertrimmed2=False,
         outF=False, outR=False,
         outCDR3=False,
         Cmiss1=False, Cmiss2=False,
         JmissF=False, JmissR=False, 
         VmissF=False, VmissR=False,
         outVsegF=False, outVsegR=False,
         adapters=False,
         labels=False,
         aligners=False,
         aligner_paths=False,
         wordlength=False,
         mincols=False,
         id=False,
         evalue=False,
         gapopen=False,
         gapextend=False,
         title_split=' ',
         tiebreaker = ['evalue','lower'],
         Vdb_length=30,
         threads=4,
         overwrite=False,
         clusterID_position=0,
         blastdb_version=5
         ):
    """
    """
    hitstrings_dict = collections.OrderedDict()
    adaptertrimmed1_tmp = None #calls to _cleanup fail if these haven't been set
    adaptertrimmed2_tmp = None
    FRmissfiles = {}
    genes_to_call = list(dbs.keys())
    if 'V' in genes_to_call and 'C' in genes_to_call:
        FRmissfiles['V']=(VmissF,VmissR)
    elif 'V' in genes_to_call:
        FRmissfiles['V']=(False,False) #V misses are only FR-sorted if the V search can take FR-sorted files as input
    else:
        pass
    if 'J' in genes_to_call: FRmissfiles['J']=(JmissF,JmissR)
    if 'C' in genes_to_call: FRmissfiles['C']=(False,False) #C misses are not FR-sorted
    
    def _cleanup(notrim,adaptertrimmed1_tmp,adaptertrimmed2_tmp,adaptertrimmed1,adaptertrimmed2,to_remove):
        if not notrim:
            if adaptertrimmed1: reptools.fileSystemSafeMove(adaptertrimmed1_tmp, adaptertrimmed1)
            if adaptertrimmed2: reptools.fileSystemSafeMove(adaptertrimmed2_tmp, adaptertrimmed2)
        reptools.clean_up(to_remove)
    
    #check for existance of target files, to prevent overwrite
    to_check = [outF,outR,Cmiss1,Cmiss2,JmissF,JmissR,VmissF,VmissR,outVsegF,outVsegR]
    if not notrim: to_check.extend([adaptertrimmed1,adaptertrimmed2])
    if not overwrite:
        for fn in to_check:
            if fn:
                if os.path.exists(fn):
                    raise IOError('overwrite=False, but target file {} already exists'.format(fn))
    
    gene_dict = reptools.genedictionaryreader(genedictfile)
    to_remove = []
    if not noCDR3 and not outCDR3:
        raise ValueError('A target outCDR3 file must be supplied unless --noCDR3 is set.  If raw CDR3 '
                         'output is not required, but denoised CDR3 output is, set outCDR3 to a temporary file.')
    
    if noCDR3: outCDR3=None
    
    #fill missing defaults from dictionaries
    
    aligners = fill_defaults(aligners, ['C','J','V','C104'], ['blastn','swipe','blastn','swipe'])
    aligner_paths = fill_defaults(
                                  aligner_paths,
                                  ['blastn','swipe','vsearch','makeblastdb','bbduk'],
                                  ['blastn','swipe','vsearch','makeblastdb','bbduk.sh']
                                  )
    mincols = fill_defaults(mincols, ['C','J','V','V_Fread','C104'], [15,10,50,80,20])
    id = fill_defaults(id, ['C','J','V','C104'], [0.9,0.9,0.9,0.93])
    evalue = fill_defaults(evalue, ['C','J','V','C104'], [0.001,0.001,0.001,0.001])
    wordlength = fill_defaults(wordlength, ['C','J','V','C104'], [9,9,9,9])
    labels = fill_defaults(labels, ['V','J','C'], ['V','J','C'], ordered=True)
    
    #Check genes_to_call and CDR3 switches for valid combination
    if 'V' not in genes_to_call and 'C' not in genes_to_call:
        raise ValueError("At least one of 'V' and 'C' must be present in 'dbs' in call to call_filepair()")
    
    if not noCDR3 and ('V' not in genes_to_call or 'J' not in genes_to_call):
        raise ValueError(
                         "CDR3 cannot be output unless both 'V' and 'J' are present in 'dbs' in call to "
                         "call_filepair()"
                         )
    
    if outVsegF or outVsegR and 'V' not in genes_to_call:
        raise ValueError(
                         "If either 'outVsegF' or 'outVsegR' are set, 'V' must be present in 'dbs' "
                         "in call to call_filepair()"
                         )
    
    #check input values
    if (
        not notrim
        and
        ((adaptertrimmed1 and not adaptertrimmed2) or (adaptertrimmed2 and not adaptertrimmed1))
        ):
        raise ValueError('If either one of adaptertrimmed1 and adaptertrimmed2 are supplied, both must be.')
    
    if not notrim and not adapters:
        raise ValueError('If trimAdapters is a specified action (the default), adapters must be supplied.')
    
    #trim adapters
    if not notrim:
        adaptertrimmed1_tmp = reptools.createTempfile(suffix='.fastq')
        adaptertrimmed2_tmp = reptools.createTempfile(suffix='.fastq')
        to_remove.extend([adaptertrimmed1_tmp,adaptertrimmed2_tmp])
        reptools.remove_adapter_bbduk(
                                infile=in1, pairfile=in2,
                                outfile=adaptertrimmed1_tmp, outpair=adaptertrimmed2_tmp,
                                adapters=adapters, adapter_rc=False,
                                discard_trimmed=False,
                                bbduk_path=aligner_paths['bbduk']
                               )
        in1 = adaptertrimmed1_tmp
        in2 = adaptertrimmed2_tmp
    
    readIDs = reptools.getReadIDs(in1,in2,title_split=title_split,clusterID_position=clusterID_position)
    
    #Csearch to identify C genes and output temporary _F and _R reads (and misses)
    #Cgene is a list of strings of the same length as the input file
    if Cmiss1 or Cmiss2:
        write_Cmisses = True
    else:
        write_Cmisses = False
    
    if min(os.path.getsize(in1),os.path.getsize(in2)) == 0: #empty input file (after trimming), so return immediately
        _cleanup(notrim,adaptertrimmed1_tmp,adaptertrimmed2_tmp,adaptertrimmed1,adaptertrimmed2,to_remove)
        return(None)
    
    if 'C' in genes_to_call:
        hitstrings_dict['C'], tempF, tempR, tempNA_C1, tempNA_C2  = reptools.Csearch(
                                                                      in1,in2,
                                                                      db=dbs['C'],
                                                                      readIDs=readIDs,
                                                                      aligner=aligners['C'],
                                                                      aligner_path=aligner_paths[aligners['C']],
                                                                      mincols=mincols['C'],
                                                                      id=id['C'],
                                                                      evalue=evalue['C'],
                                                                      wordlength=wordlength['C'],
                                                                      title_split = title_split,
                                                                      tiebreaker = tiebreaker,
                                                                      write_Cmisses = write_Cmisses,
                                                                      makeblastdbPath = aligner_paths['makeblastdb'],
                                                                      threads=threads,
                                                                      blastdb_version = blastdb_version
                                                                      )
    else:
        tempF,tempR = in1,in2 #these being the input files for Vsearch
        tempNA_C1,tempNA_C2 = False,False
    
    if 'V' in genes_to_call:
        if 'C' not in genes_to_call or len(hitstrings_dict['C'])>0:
            #hitF_v is only required either:
            #(a) for V segment output, so is False unless at least one of outVsegF and outVsegR are True
            #(b) in the absence of C gene searching, when it will be the full F read out 
            if outVsegF or outVsegR or 'C' not in genes_to_call:
                returnF=True
            else:
                returnF=False
            if 'C' not in genes_to_call:
                strand_searches = ['both','both']
                FRsorted=False
            else:
                strand_searches = ['plus','minus']
                FRsorted=True
            hitstrings_dict['V'],hitF_v,hitR_v = reptools.Vsearch(
                                            tempF, tempR,
                                            db=dbs['V'],
                                            readIDs=readIDs,
                                            aligner = aligners['V'], aligner_path = aligner_paths[aligners['V']],
                                            mincols = mincols['V'], mincols_F = mincols['V_Fread'], id = id['V'],
                                            evalue = evalue['V'], wordlength=wordlength['V'],
                                            title_split = title_split, 
                                            clusterID_position  = clusterID_position,
                                            preferred_read='R',
                                            tiebreaker = tiebreaker,
                                            makeblastdbPath = aligner_paths['makeblastdb'],
                                            returnF = returnF,
                                            threads = threads,
                                            strand_searches = strand_searches,
                                            FRsorted = FRsorted,
                                            blastdb_version = blastdb_version
                                            )
        else:
            hitF_v,hitR_v = False,False
            hitstrings_dict['V']={}
        if 'C' not in genes_to_call:
            tempF,tempR = hitF_v,hitR_v #these being the F/R sorted files used by emitFR
    else:
        hitF_v,hitR_v = False,False
    #because CDR3 selection involves very slow searches, want to restrict following searches to VC hits
    if 'J' in genes_to_call:
        if 'V' in genes_to_call:
            infile = hitR_v
        else:
            infile = tempR
        #Jhits will be a dictionary of reptools.hits objects containing the info of the hits closest to the junction
        if 'V' not in genes_to_call or len(hitstrings_dict['V'])>0 : #run Jsearch unless V was run but returned nothing
            Jhits,hitstrings_dict['J']  = reptools.Jsearch(
                                                            infile,
                                                            gene_dict,
                                                            dbs['J'],
                                                            aligners['J'],aligner_paths[aligners['J']],
                                                            mincols['J'],id['J'],evalue['J'],
                                                            title_split = title_split,
                                                            clusterID_position = clusterID_position,
                                                            selectby = tiebreaker,
                                                            makeblastdbPath = aligner_paths['makeblastdb'],
                                                            threads=threads,
                                                            blastdb_version=blastdb_version
                                                           )
        else:
            Jhits = {}
            hitstrings_dict['J'] = {}
    if outVsegF or outVsegR or ('V' in genes_to_call and 'J' in genes_to_call and not noCDR3):
        if len(hitstrings_dict['V'])>0 and (outVsegF or outVsegR or not noCDR3):
            #make modified V database and gene dictionary, to work with only last 30 bases (default - set by Vdb_length)
            (temp_Vdb, temp_genedictfile) = reptools.make_shorter_db(dbs['V'],genedictfile,Vdb_length)
            to_remove.extend([temp_Vdb,temp_genedictfile])
            temp_genedict = reptools.genedictionaryreader(temp_genedictfile)
            
            #C104hits will be a dictionary of reptools.hits objects containing info of the hits closest to the junction
            C104hits = reptools.V_C104search(
                                     hitR_v,
                                     temp_Vdb,
                                     temp_genedict,
                                     aligners['C104'],aligner_paths[aligners['C104']],
                                     mincols['C104'],id['C104'],evalue['C104'],
                                     title_split = title_split,
                                     clusterID_position = clusterID_position,
                                     makeblastdbPath = aligner_paths['makeblastdb'],
                                     threads = threads,
                                     blastdb_version = blastdb_version
                                     )
        else:
            C104hits={}
    
    #if V segment output is required, search the F read with SWIPE as well
    if len(hitstrings_dict['V'])>0 and (outVsegF or outVsegR):
        #(temp_Vdb, temp_genedictfile) = reptools.reptools.make_shorter_db(dbs['V'],genedictfile,Vdb_length)
        #to_remove.extend([temp_Vdb,temp_genedictfile])
        #temp_genedict = reptools.reptools.genedictionaryreader(temp_genedictfile)
        for_C104hits = reptools.V_C104search(
                                         hitF_v,
                                         temp_Vdb,
                                         temp_genedict,
                                         aligners['C104'],aligner_paths[aligners['C104']],
                                         mincols['C104'],id['C104'],evalue['C104'],
                                         title_split = title_split,
                                         clusterID_position = clusterID_position,
                                         makeblastdbPath = aligner_paths['makeblastdb'],
                                         strand = 'plus',
                                         threads = threads,
                                         blastdb_version = blastdb_version
                                         )
        
        for_Vhits = reptools.V_C104search(
                                 hitF_v,
                                 dbs['V'],
                                 gene_dict,
                                 aligners['V'],aligner_paths[aligners['V']],
                                 mincols['V'],id['V'],evalue['V'],
                                 title_split = title_split,
                                 clusterID_position = clusterID_position,
                                 makeblastdbPath = aligner_paths['makeblastdb'],
                                 strand='plus',
                                 threads=threads,
                                 blastdb_version = blastdb_version
                                 )
        
        rev_Vhits = reptools.V_C104search(
                                 hitF_v,
                                 dbs['V'],
                                 gene_dict,
                                 aligners['V'],aligner_paths[aligners['V']],
                                 mincols['V'],id['V'],evalue['V'],
                                 title_split = title_split,
                                 clusterID_position = clusterID_position,
                                 makeblastdbPath = aligner_paths['makeblastdb'],
                                 strand = 'minus',
                                 threads = threads,
                                 blastdb_version = blastdb_version
                                 )
    
    #output fastq and raw CDR3 files
    #raise ValueError(hitstrings_dict)
    reptools.emitFR(
                    inF = tempF,
                    inR = tempR,
                    outF = outF,
                    outR = outR,
                    genes = collections.OrderedDict(
                                    [
                                     (
                                       gene,
                                       {
                                         'Fmiss':FRmissfiles[gene][0],
                                         'Rmiss':FRmissfiles[gene][1],
                                         'hits':hitstrings_dict[gene]
                                        }
                                      )
                                      for gene in hitstrings_dict 
                                     ]
                                     ),
                    labels = labels,
                    title_split = title_split,
                    clusterID_position = clusterID_position
                    )
    
    if len(hitstrings_dict['V'])>0: #if no V hits were found, no CDR3 or V segment to write
        if not noCDR3 and len(hitstrings_dict['J'])>0:
            temp_outCDR3 = reptools.createTempfile(suffix='.fastq')
            to_remove.append(temp_outCDR3)
            temp_outCDR3 = reptools.emitCDR3(
                                              temp_outCDR3,
                                              hitR_v,
                                              hitstrings_dict,
                                              Jhits,
                                              C104hits,
                                              gene_dict,
                                              temp_genedict,
                                              title_split=title_split,
                                              clusterID_position=clusterID_position,
                                              labels=labels
                                              )
        else:
            temp_outCDR3 = None
        if outVsegF or outVsegR:
            reptools.emitVsegs(
                      outVsegF,
                      outVsegR,
                      hitF_v,
                      hitR_v,
                      hitstrings_dict['C'],
                      hitstrings_dict['J'],
                      hitstrings_dict['V'],
                      for_Vhits,
                      rev_Vhits,
                      C104hits,
                      for_C104hits,
                      temp_genedict,
                      gene_dict,
                      title_split=title_split,
                      clusterID_position=clusterID_position,
                      labels=labels
                      )
    else:
        temp_outCDR3 = None
    
    if temp_outCDR3 is None:
        outCDR3 = None
    
    if Cmiss1 and os.path.isfile(tempNA_C1): reptools.fileSystemSafeMove(tempNA_C1, Cmiss1)
    if Cmiss2 and os.path.isfile(tempNA_C2): reptools.fileSystemSafeMove(tempNA_C2, Cmiss2)
    if not noCDR3 and outCDR3 is not None: reptools.fileSystemSafeMove(temp_outCDR3,outCDR3)
    
    to_remove.extend([tempF,tempR,hitR_v])
    _cleanup(notrim,adaptertrimmed1_tmp,adaptertrimmed2_tmp,adaptertrimmed1,adaptertrimmed2,to_remove)
    
    return(outCDR3)


def Csearch(in1,
            in2,
            readIDs,
            db,
            aligner,
            aligner_path,
            mincols,
            id,
            evalue,
            wordlength,
            title_split=' ',
            clusterID_position=0,
            tiebreaker = ['evalue','lower'],
            write_Cmisses=False,
            makeblastdbPath='makeblastdb',
            threads=False,
            blastdb_version=5
            ):
    
    toremove = []
    tempNA_1 = None
    tempNA_2 = None
    
    Writer,Parser,out_ext = reptools.initialchecks([in1,in2],aligner)
    Caller,hits_ext = reptools.setCaller(aligner)
    
    tempF = reptools.createTempfile(suffix=out_ext)
    tempR = reptools.createTempfile(suffix=out_ext)
    if write_Cmisses:
        tempNA_1 = reptools.createTempfile(suffix=out_ext)
        tempNA_2 = reptools.createTempfile(suffix=out_ext)
    
    hitfiles = [reptools.createTempfile(suffix=hits_ext) for fn in [in1,in2]]
    toremove.extend(hitfiles)
    
    for fn,hitfile in zip([in1,in2],hitfiles):
        Caller(
                    infile = fn,
                    db = db,
                    type = aligner,
                    evalue = evalue,
                    hitsOut=hitfile,
                    strand='both',
                    #id = id, #passing id was missing perfectly good hits - so filter top hits by id
                    #mincols = mincols,
                    output_no_hits = False,
                    top_hits_only = True,
                    top_hit_only = False,
                    maxaccepts = '0',
                    maxrejects = '0',
                    swipePath = aligner_path,
                    blastPath = aligner_path,
                    makeblastdbPath = makeblastdbPath,
                    wordlength = wordlength,
                    threads = threads,
                    #gapopen = gapopen,
                    #gapextend = gapextend,
                    #reward = reward,
                    #penalty = penalty,
                    #targetcov = targetcov,
                    #minhsp = minhsp,
                    #fulldp = fulldp,
                    #xDrop = xDrop,
                    #alphabet = alphabet,
                    #numMatches = numMatches,
                    verbose = True,
                    blastdb_version=blastdb_version
                )
    
    #parse the files
    forstrand,tophits = reptools.parseCV(
                                infiles={'1':hitfiles[0],'2':hitfiles[1]},
                                readIDs=readIDs,
                                mincols={'1':mincols,'2':mincols},
                                id=id,
                                aligner=aligner,
                                title_split=title_split,
                                clusterID_position=clusterID_position,
                                tiebreaker=tiebreaker
                                )
    
    #write the temporary output files
    reptools.writeFRsorted_pairfiles(
                                    in1=in1,
                                    in2=in2,
                                    hitF=tempF,
                                    hitR=tempR,
                                    miss1=tempNA_1,
                                    miss2=tempNA_2,
                                    forstrand=forstrand,
                                    Parser=Parser,
                                    Writer=Writer,
                                    writemisses=write_Cmisses
                                    )
    
    reptools.clean_up(toremove)
    
    return(tophits, tempF, tempR, tempNA_1, tempNA_2)


def writeFRsorted_pairfiles(in1,in2,hitF,hitR,miss1,miss2,forstrand,Parser,Writer,writemisses=False):
    with (open(miss1, 'w') if writemisses else dummy()) as misshandle_1:
        with (open(miss2, 'w') if writemisses else dummy()) as misshandle_2:
            with open(hitF, 'w') as hithandle_F, open(hitR, 'w') as hithandle_R:
                with open(in1) as inhandle1, open(in2) as inhandle2:
                    for seq_tuple1,seq_tuple2,F in zip(Parser(inhandle1),Parser(inhandle2),forstrand):
                        if F==1:
                            Writer(hithandle_F,seq_tuple1)
                            Writer(hithandle_R,seq_tuple2)
                        elif F==2:
                            Writer(hithandle_F,seq_tuple2)
                            Writer(hithandle_R,seq_tuple1)
                        else:
                            Writer(misshandle_1,seq_tuple1)
                            Writer(misshandle_2,seq_tuple2)


def write_unsorted_file(infile,outfile,hits,title_split,clusterID_position,Parser,Writer):
    with open(outfile, 'w') as out_handle:
        with open(infile) as in_handle:
            for seq_tuple in Parser(in_handle):
                if seq_tuple[0].split(title_split)[clusterID_position] in hits:
                    Writer(out_handle,seq_tuple)


def Vsearch(
            inF,
            inR,
            db,
            readIDs,
            aligner,
            aligner_path,
            mincols,
            mincols_F,
            id,
            evalue,
            wordlength,
            title_split=' ',
            clusterID_position=0,
            preferred_read='R',
            top_hits_only=True,
            makeblastdbPath='makeblastdb',
            returnF = False,
            threads=False,
            FRsorted=True,
            tiebreaker = ['evalue','lower'],
            strand_searches = ['plus','minus'],
            blastdb_version=5
            ):
    
    toremove = []
    
    Writer,Parser,out_ext = reptools.initialchecks([inF,inR],aligner)
    Caller,hits_ext = reptools.setCaller(aligner)
    
    hitfiles = [reptools.createTempfile(suffix=hits_ext) for fn in [inF,inR]]
    toremove.extend(hitfiles)
    
    for fn,strand,hitfile in zip([inF,inR],strand_searches,hitfiles):
        Caller( #BLAST+ does not have a top_hits_only option, so will need to filter posthoc
                    infile = fn,
                    db = db,
                    type = aligner,
                    evalue = evalue,
                    hitsOut=hitfile,
                    strand=strand,
                    output_no_hits = False,
                    top_hits_only = top_hits_only,
                    top_hit_only = False,
                    maxaccepts = '0',
                    maxrejects = '0',
                    swipePath = aligner_path,
                    blastPath = aligner_path,
                    makeblastdbPath=makeblastdbPath,
                    wordlength = wordlength,
                    threads = threads,
                    #gapopen = gapopen,
                    #gapextend = gapextend,
                    #reward = reward,
                    #penalty = penalty,
                    #targetcov = targetcov,
                    #minhsp = minhsp,
                    #fulldp = fulldp,
                    #xDrop = xDrop,
                    #alphabet = alphabet,
                    #numMatches = numMatches,
                    verbose = True,
                    blastdb_version=blastdb_version
                )
    #parse the files 
    if top_hits_only and not (aligner.lower().startswith('vsearch') or aligner.lower().startswith('usearch')):
        parse_top_hits_only = True
    else:
        parse_top_hits_only = False # usearch/vsearch can return joint top hits, so further parsing is not required
    
    if FRsorted:
        infiles={'F':hitfiles[0],'R':hitfiles[1]}
        mincolsdict={'F':mincols_F,'R':mincols}
    else:
        infiles={'1':hitfiles[0],'2':hitfiles[1]}
        #if the input files aren't already seperated into F and R, we can't apply a different mincols threshold to F
        mincolsdict={'1':mincols,'2':mincols}
    
    forstrand, Vhits = reptools.parseCV(
                                        infiles=infiles,
                                        readIDs=readIDs,
                                        mincols=mincolsdict,
                                        id=id,
                                        aligner=aligner,
                                        title_split=title_split,
                                        clusterID_position=clusterID_position,
                                        preferred_read=preferred_read,
                                        top_hits_only=parse_top_hits_only,
                                        tiebreaker=tiebreaker
                                        )
    
    
    #write the output files
    #if the input files are already FRsorted, simply output those with hits:
    if FRsorted:
        outR = reptools.createTempfile(suffix=out_ext)
        reptools.write_unsorted_file(
                                     infile=inR,
                                     outfile=outR,
                                     hits=Vhits,
                                     title_split=title_split,
                                     clusterID_position=clusterID_position,
                                     Parser=Parser,
                                     Writer=Writer
                                     )
        #write a paired output file only if returnF=True
        if returnF:
            outF = reptools.createTempfile(suffix=out_ext)
            reptools.write_unsorted_file(
                                     infile=inF,
                                     outfile=outF,
                                     hits=Vhits,
                                     title_split=title_split,
                                     clusterID_position=clusterID_position,
                                     Parser=Parser,
                                     Writer=Writer
                                     )
        else:
            outF = False
    
    else: # if not FR sorted
        outF = reptools.createTempfile(suffix=out_ext)
        outR = reptools.createTempfile(suffix=out_ext)
        reptools.writeFRsorted_pairfiles(
                                in1=inF,
                                in2=inR,
                                hitF=outF,
                                hitR=outR,
                                miss1=False,
                                miss2=False,
                                forstrand=forstrand,
                                Parser=Parser,
                                Writer=Writer,
                                writemisses=False
                                )
    reptools.clean_up(toremove)
    
    return(Vhits, outF, outR)


def V_C104search(
                 tempR_v,
                 db,
                 gene_dict,
                 aligner,aligner_path,
                 mincols,id,evalue,
                 title_split,
                 clusterID_position,
                 makeblastdbPath='makeblastdb',
                 strand='minus',
                 threads=False,
                 blastdb_version=5
                 ):
    toremove = []
    
    Writer,Parser,out_ext = reptools.initialchecks([tempR_v],aligner)
    Caller,hits_ext = reptools.setCaller(aligner)
    
    hitfile = reptools.createTempfile(suffix=hits_ext)
    toremove.append(hitfile)
    
    Caller( #need to sort out defining the aligner path
                infile = tempR_v,
                db = db,
                type = aligner,
                #id = id,
                evalue = evalue,
                hitsOut=hitfile,
                strand=strand,
                #mincols=mincols,
                output_no_hits = False,
                top_hits_only = False, ###NOTE THIS LINE###
                top_hit_only = False,
                maxaccepts = '0',
                maxrejects = '0',
                swipePath = aligner_path,
                blastPath = aligner_path,
                makeblastdbPath=makeblastdbPath,
                threads = threads,
                #wordlength = wordlength,
                #gapopen = gapopen,
                #gapextend = gapextend,
                #reward = reward,
                #penalty = penalty,
                #targetcov = targetcov,
                #minhsp = minhsp,
                #fulldp = fulldp,
                #xDrop = xDrop,
                #alphabet = alphabet,
                #numMatches = numMatches,
                verbose = True,
                blastdb_version=blastdb_version
            )
    
    C104hits = reptools.parseC104(hitfile,gene_dict,aligner,title_split=title_split,
                                  clusterID_position=clusterID_position,mincols=mincols,id=id)
    
    reptools.clean_up(toremove)
    
    return(C104hits)


def Jsearch(
             inR,gene_dict,db,
             aligner,aligner_path,
             mincols,id,evalue,
             makeblastdbPath='makeblastdb',
             title_split=' ',
             clusterID_position=0,
             top_hits_only=True,
             selectby=['evalue','low'],
             threads=False,
             blastdb_version=5
            ):
    
    toremove = []
    
    Writer,Parser,out_ext = reptools.initialchecks([inR],aligner)
    Caller,hits_ext = reptools.setCaller(aligner)
    
    hitfile = reptools.createTempfile(suffix=hits_ext)
    toremove.append(hitfile)
    
    Caller( #need to sort out defining the aligner path
                infile = inR,
                db = db,
                type = aligner,
                #id = id,
                evalue = evalue,
                hitsOut=hitfile,
                strand='minus',
                #mincols=mincols,
                output_no_hits = False,
                top_hits_only = top_hits_only,
                top_hit_only = False,
                maxaccepts = '0',
                maxrejects = '0',
                threads = threads,
                #wordlength = wordlength,
                #gapopen = gapopen,
                #gapextend = gapextend,
                #reward = reward,
                #penalty = penalty,
                #targetcov = targetcov,
                #minhsp = minhsp,
                #fulldp = fulldp,
                #xDrop = xDrop,
                #alphabet = alphabet,
                #numMatches = numMatches,
                verbose = True,
                swipePath = aligner_path,
                blastPath = aligner_path,
                blastdb_version=blastdb_version
            )
    
    #parse the files 
    Jhits,Jhitstrings = reptools.parseJ(
                                hitfile,
                                gene_dict,
                                aligner,
                                title_split=title_split,
                                clusterID_position=clusterID_position,
                                top_hits_only=top_hits_only,
                                selectby=selectby,
                                mincols=mincols, id=id
                                )
    
    reptools.clean_up(toremove)
    
    return(Jhits,Jhitstrings)


def parseCV(
            infiles,
            readIDs,
            mincols,
            id,
            aligner,
            title_split=' ',
            clusterID_position=0,
            preferred_read='R',
            top_hits_only=True,
            tiebreaker = ['evalue','lower']
            ):
    """
    produces a (optionally) a list and a dictionary
        If the infiles are '1' and '2', 'forstrand' is a list which shows read is the forward strand (0 for no hit or
        conflicting hits).  It has the same order as the input files.  If the infiles are 'F' and 'R', None is returned
        for 'forstrand'.
        'tophits' is a dictionary giving the hit name, or when there are multiple top hits all of them,
        seperated by "+" (empty string for no hit or conflicting hits)
        'infiles' should be a dictionary of 2 file names (keyed '1'&'2' or 'F'&'R')
        'mincols' should be a dictionary with the same ids as the file dictionary
        'preferred_read' should be F or R.  Gene segment IDs will be taken from this read, with the other
        read used to reduce ambiguity (default R).
    """
    if set(list(infiles.keys()))==set(['1','2']):
        FRsorted=False
    elif set(list(infiles.keys()))==set(['F','R']):
        FRsorted=True
    else:
        raise KeyError('The infiles dictionary supplied to parseCV() must have the keys "1"/"2" or "F"/"R"')
    
    if preferred_read=='R':
        paired_read='F'
    elif preferred_read=='F':
        paired_read='R'
    else:
        raise KeyError(
       'preferred_read must be set to "F" or "R"'
          )
    
    forstrand = None
    tophits = {}
    hitinfo = {}
    for readID in infiles: #parse the infiles
        if top_hits_only:
            hitinfo[readID] = reptools.get_top_hits(
                                        infiles[readID],
                                        filetype=aligner,
                                        title_split=title_split,
                                        clusterID_position=clusterID_position,
                                        criteria=tiebreaker,
                                        id=id,
                                        mincols=mincols[readID]
                                        )
        else:
            hitinfo[readID] = reptools.parseHitfile(
                                        infiles[readID],
                                        aligner,
                                        title_split=title_split,
                                        clusterID_position=clusterID_position,
                                        id=id,
                                        mincols=mincols[readID]
                                        )
    
    #if the input files weren't pre-sorted into F and R, do that here, and create a list to reference which read, while
    #also allocating read 1 and read 2 hits to F/R
    if not FRsorted:
        forstrand = []
        hitinfo['F'] = {}
        hitinfo['R'] = {}
        for seq in readIDs: #for each readpair
            if not (sum([1 for readID in infiles if seq in hitinfo[readID]])): #no hits
                forstrand.append(0)
            else: #there is a hit on at least one strand
                if seq in hitinfo['1'] and seq not in hitinfo['2']: #hit only in read_1
                    result = hitinfo['1'][seq]
                    if result[0].hsps[0].strand == 1:
                        forstrand.append(1)
                        hitinfo['F'][seq]=result
                    else:
                        forstrand.append(2)
                        hitinfo['R'][seq]=result
                elif seq in hitinfo['2'] and seq not in hitinfo['1']: #hit only in read_2
                    result = hitinfo['2'][seq]
                    if result[0].hsps[0].strand == 1:
                        forstrand.append(2)
                        hitinfo['F'][seq]=result
                    else:
                        forstrand.append(1)
                        hitinfo['R'][seq]=result
                else: #hit in both reads
                    tiebreaker_choice = None #initialise this, so that it doesn't ever have to be calculated twice
                    result1 = hitinfo['1'][seq]
                    result2 = hitinfo['2'][seq]
                    #do strands agree?
                    if result1[0].hsps[0].strand * result2[0].hsps[0].strand == 1: #clash, so consult tiebreaker
                        tiebreaker_choice = reptools.tiebreak(result1,result2,tiebreaker)
                        if tiebreaker_choice == 0:
                            forstrand.append(0)
                        elif tiebreaker_choice == 1:
                            if result1[0].hsps[0].strand == 1:
                                forstrand.append(1)
                                hitinfo['F'][seq]=result1
                                hitinfo['R'][seq]=result2
                            else:
                                forstrand.append(2)
                                hitinfo['F'][seq]=result2
                                hitinfo['R'][seq]=result1
                        elif tiebreaker_choice == 2:
                            if result2[0].hsps[0].strand == 1:
                                forstrand.append(2)
                                hitinfo['F'][seq]=result2
                                hitinfo['R'][seq]=result1
                            else:
                                forstrand.append(1)
                                hitinfo['F'][seq]=result1
                                hitinfo['R'][seq]=result2
                        else: 
                            raise ValueError("Likely bug: tiebreaker_choice wasn't set correctly.")
                    else: #the two strands are different, so no contradiction
                        if result1[0].hsps[0].strand == 1:
                            forstrand.append(1)
                            hitinfo['F'][seq]=result1
                            hitinfo['R'][seq]=result2
                        else:
                            forstrand.append(2)
                            hitinfo['F'][seq]=result2
                            hitinfo['R'][seq]=result1
                            
    #get ID strings
    #TODO - consult gene dictionary here for segment/alleles name, in order to allow multiple transcripts for one allele
    for seq in hitinfo[preferred_read]:
        result = hitinfo[preferred_read][seq]
        if len(result)==0: continue
        hits_preferred = set([hit.id for hit in result])
        if len(hits_preferred) == 1: #there is only one top hit on the R read
            tophits[seq] = ''.join(hits_preferred)
        else: #as there is ambiguity, use the other read to narrow down ids
            if seq in hitinfo[paired_read]: #but only if there were hits on the other read, of course
                paired_result = hitinfo[paired_read][seq]
                hits_paired = set([hit.id for hit in paired_result])
                hits_in_common = hits_preferred.intersection(hits_paired)
                if len(hits_in_common)>0: #use the hits the two reads have in common
                    tophits[seq] = '+'.join(sorted(hits_in_common))
                else: #unless they have none in common, in which case take only the hits from the selected read
                   tophits[seq] = '+'.join(sorted(hits_preferred))
            else: #no hits on the other read, so take the hits from the selected read
                tophits[seq] = '+'.join(sorted(hits_preferred))
    
    return(forstrand,tophits)


def parseJ(
               hitfile,gene_dict,aligner,
               mincols=False,id=False,
               title_split=' ',
               clusterID_position=0,
               top_hits_only=True,
               selectby=['evalue','low']
           ):
    """
    """
    closest_hits,hitstrings = reptools.find_nearest_to_overlap_hits(
                                                          hitfile,
                                                          'J',
                                                          gene_dict,
                                                          'JregionF118start',
                                                          filetype=aligner,
                                                          title_split=title_split,
                                                          clusterID_position=clusterID_position,
                                                          top_hits_only=top_hits_only,
                                                          returnhitstrings=True,
                                                          mincols=mincols,id=id
                                                  )
    return(closest_hits,hitstrings)


def parseC104(hitfile,gene_dict,aligner,mincols=False,id=False,title_split=' ',clusterID_position=0):
    """
    """
    #hitinfo = reptools.reptools.parseHitfile(hitfile,aligner,title_split=title_split)
    #this currently relies on bitscore to break ties - need to check that it is working correctly with blastn and swipe
    #?this is failing to parse final read in a usearch file #CHECK
    #todo = OPTION TO SELECT TOP HITS ONLY?  (SWIPE doesn't have a top hits only setting)
    closest_hits = reptools.find_nearest_to_overlap_hits( 
                                                  hitfile,
                                                  'V',
                                                  gene_dict,
                                                  'VregionC104start',
                                                  filetype=aligner,
                                                  title_split=title_split,
                                                  clusterID_position=clusterID_position,
                                                  mincols=mincols,id=id
                                              )
    return(closest_hits)


def emitFR(
            inF,
            inR,
            outF,
            outR,
            genes,
            labels,
            title_split = ' ',
            clusterID_position = 0
    ):
    """
    Input:
        'inF','inR' = forward and reverse input 
        'outF', 'outR' = destination for pairs with hits to all supplied genes
        'genes' = an ordered dictionary containing a key for each gene to be reported, with each entry being a
                  dictionary containing 'hits' (a dictionary of hit strings keyed by SeqID), 'Fmiss' (destination for F 
                  misses, or None), 'Rmiss' (destination for R misses, or None).  The order informs assumed filtering 
                  order (in other words, if a sequence is saved as a miss for the first gene, it won't be saved as a 
                  miss for subsequent genes'
        'title_split' = a string
        'clusterID_position' = an integer
        'labels' = an ordered dictionary containing the labels for each gene.  The order informs label display order.
    """
    if os.path.getsize(inF)==0 or os.path.getsize(inR)==0:
        return()
    
    Writer,Parser,out_ext = reptools.initialchecks([inF,inR])
    
    #remove any entries from 'genes' where the value for 'hits' is None, and remove values from 'labels' which aren't
    #in 'genes'
    genes = collections.OrderedDict({k:genes[k] for k in genes if genes[k]['hits'] is not None})
    labels = collections.OrderedDict({k:labels[k] for k in labels if k in genes})
    
    with open(outF,'w') if outF else dummy() as outF_h, open(outR,'w') if outR else dummy() as outR_h:
        with contextlib.ExitStack() as outstack:
            Fmiss_handles = {}
            Rmiss_handles = {}
            for gene in genes:
                if genes[gene]['Fmiss']:
                    Fmiss_handles[gene] = outstack.enter_context(open(genes[gene]['Fmiss'],'w'))
                else:
                    Fmiss_handles[gene] = outstack.enter_context(dummy())
                if genes[gene]['Rmiss']:
                    Rmiss_handles[gene] = outstack.enter_context(open(genes[gene]['Rmiss'],'w'))
                else:
                    Rmiss_handles[gene] = outstack.enter_context(dummy())
            
            with open(inF) as inhandleF, open(inR) as inhandleR:
                for seq_tupF,seq_tupR in zip(Parser(inhandleF),Parser(inhandleR)):
                    miss=False
                    titlestem = seq_tupF[0].split(title_split)[clusterID_position]
                    
                    #build title string
                    segments = {}
                    for gene in genes:
                        if titlestem in genes[gene]['hits']:
                            segments[gene] = genes[gene]['hits'][titlestem]
                        else:
                            segments[gene] = 'None'
                            break #subsequent genes aren't recorded at all (even as None), as they haven't been searched
                    
                    gene_string = '{};'.format(
                                               ';'.join(
                                                        [
                                                            '{}={}'.format(labels[gene],segments[gene]) 
                                                            for gene in labels
                                                            if gene in segments
                                                        ]
                                                        )
                                               )
                    
                    for gene in genes: #loop through genes in filtering order, outputting as soon as a miss is found
                        if titlestem not in genes[gene]['hits']:
                            miss=True
                            Writer(Fmiss_handles[gene], tuple(['{};{}'.format(seq_tupF[0],gene_string)]) + seq_tupF[1:])
                            Writer(Rmiss_handles[gene], tuple(['{};{}'.format(seq_tupR[0],gene_string)]) + seq_tupR[1:])
                            break
                    if not miss:
                        Writer(outF_h, tuple(['{};{}'.format(seq_tupF[0],gene_string)]) + seq_tupF[1:])
                        Writer(outR_h, tuple(['{};{}'.format(seq_tupR[0],gene_string)]) + seq_tupR[1:])
    
    _ = reptools.removeemptyfilepair(outF,outR)
    for gene in genes:
        if genes[gene]['Fmiss'] or genes[gene]['Rmiss']:
            _ = reptools.removeemptyfilepair(genes[gene]['Fmiss'], genes[gene]['Rmiss'])


def emitCDR3(
             outCDR3,
             tempR_v,
             hitstrings_dict,
             Jhits,
             C104hits,
             genedict,
             temp_genedict,
             title_split,
             clusterID_position,
             labels
             ):
    hitstrings_dict = collections.OrderedDict({k:v for k,v in hitstrings_dict.items() if v is not None})
    labels = collections.OrderedDict({k:v for k,v in labels.items() if k in hitstrings_dict })
    
    Writer,Parser,out_ext = reptools.initialchecks([tempR_v])
    with open(outCDR3,'w') as outhandle:
        with open(tempR_v) as inhandleR:
            for seq_tuple in Parser(inhandleR):
                titlestem = seq_tuple[0].split(title_split)[clusterID_position]
                try: #is there a result for each, with matching strand directions?
                    #if the strand directions don't match, skip this sequence
                    if Jhits[titlestem].strand==C104hits[titlestem].strand: 
                        direction = Jhits[titlestem].strand
                        sliced = tuple([
                                    reptools.get_position(
                                               C104hits[titlestem],
                                               pos_H=temp_genedict['V']['VregionC104start'][C104hits[titlestem].hit_id],
                                               gene='V'
                                                 ),
                                    reptools.get_position(
                                               Jhits[titlestem],
                                               pos_H=genedict['J']['JregionF118start'][Jhits[titlestem].hit_id],
                                               gene='J'
                                                 )
                                   ])
                        
                        #make full gene string
                        segments = {}
                        for gene in hitstrings_dict:
                            if titlestem in hitstrings_dict[gene]:
                                segments[gene] = hitstrings_dict[gene][titlestem]
                            else:
                                segments[gene] = 'None'
                                break #subsequent genes aren't recorded at all (even as None): no search performed
                        
                        gene_string = '{};'.format(
                           ';'.join(
                                    [
                                        '{}={}'.format(labels[gene],segments[gene]) 
                                        for gene in labels
                                        if gene in segments
                                    ]
                                    )
                           )
                           
                        #make tuple
                        seq_list = ['{};{}'.format(seq_tuple[0],gene_string)]
                        start = sliced[0] -1
                        end = sliced[1] + (-2*(direction==-1))  + (2*direction)
                        if direction == -1 and end <0: #fix for taking the first character in a reverse slice
                            end = -len(seq_tuple[1])-1
                        #add 2 to slice for forward reads, subtract 4 for reverse reads
                        CDR3 = seq_tuple[1][start:end:direction]  #extract CDR3 seq
                        if direction==-1:
                            CDR3 = reptools.complement(CDR3).upper()
                        if len(CDR3)>0:
                            seq_list.append(CDR3)
                            if len(seq_tuple)==3:
                                seq_list.append(seq_tuple[2][start:end:direction]) #extract CDR3 qual
                        else: #if the slice was zero or negative, do not emit this read
                            continue
                        Writer(outhandle,seq_list)
                except KeyError: #if there's no hit for either V or J
                    pass    #may want to add the option to write fails to disk
    return(reptools.removeemptyfile(outCDR3)) #returns None if the file was empty (and has been removed), else the fn


def emitVsegs(
              outVsegF,outVsegR,tempF_v,tempR_v,
              Chits,Jhitstrings,Vhits,for_Vhits,rev_Vhits,C104hits,for_C104hits,
              temp_genedict,genedict,title_split,clusterID_position,labels
              ):
    Writer,Parser,out_ext = reptools.initialchecks([tempF_v])
    
    with open(outVsegF,'w') if outVsegF else dummy() as outhandF, open(outVsegR,'w') if outVsegR else dummy() as outhandR:
        with open(tempF_v) as inhandleF, open(tempR_v) as inhandleR:
            for seq_tupleF,seq_tupleR in zip(Parser(inhandleF),Parser(inhandleR)):
                titlestemF = seq_tupleF[0].split(title_split)[clusterID_position]
                titlestemR = seq_tupleR[0].split(title_split)[clusterID_position]
                if titlestemF!=titlestemR:
                    raise ValueError(
                            'Title disagreement between paired reads:\n{}\n{}\n'.format(seq_tupleF[0],seq_tupleF[1])
                            )
                else:
                    titlestem = titlestemF
                
                #if not any(
                #           titlestem in for_C104hits, titlestem in C104hits,
                #           titlestem in for_Vhits, titlestem in Vhits):
                #    raise ValueError('No hits found for read:{}'.format(titlestem))
                #    continue #no hits at all for this read (shouldn't be the case!)
                
                Vseq_pair = []
                Vqual_pair = []
                for SWIPE,BLAST,seqtuple in zip(
                                                  [for_C104hits,C104hits],
                                                  [for_Vhits,rev_Vhits],
                                                  [seq_tupleF,seq_tupleR]
                                                  ):
                    try:
                        hit = SWIPE[titlestem]
                        gd = temp_genedict['V']
                    except KeyError:
                        try:
                            hit = BLAST[titlestem]
                            gd = genedict['V']
                        except KeyError:
                            hit = False
                    if hit:
                        start,end,direction = VsegSlicer(hit,gd)
                        Vseg_seq = seqtuple[1][start:end]  #extract Vseg seq
                        if len(Vseg_seq)>0:
                            Vseq_pair.append(Vseg_seq[::direction])
                            Vseg_qual = seqtuple[2][start:end]
                            Vqual_pair.append(Vseg_qual[::direction])
                        else:
                            Vseq_pair.append('N')
                            Vqual_pair.append('!')
                    else:
                        Vseq_pair.append('N')
                        Vqual_pair.append('!')
                
                if Vseq_pair[0]=='N' and Vseq_pair[1]=='N': #if both Vsegs are length zero, do not emit this read
                    continue
                
                #make gene string
                try:
                    Jstring = Jhitstrings[titlestem]
                except KeyError:
                    Jstring = 'None'
                gene_string = '{}={};{}={};{}={}'.format(
                                                          labels['V'],
                                                          Vhits[titlestem],
                                                          labels['J'],
                                                          Jstring,
                                                          labels['C'],
                                                          Chits[titlestem]
                                                         )
                #make tuple
                seq_list = ['{};{}'.format(seq_tupleF[0],gene_string)]
                seq_listF = seq_list + [Vseq_pair[0]]
                seq_listR = seq_list + [Vseq_pair[1]]
                if len(seq_tupleF)==3: #if we're working with a FASTQ file
                    seq_listF.append(Vqual_pair[0])
                    seq_listR.append(Vqual_pair[0])
                Writer(outhandF,seq_listF)
                Writer(outhandR,seq_listR)


