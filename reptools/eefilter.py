import reptools
import os

def EEfilter_dir(
                indir,
                outdir=False,
                FASTQout_dir=None,
                FASTAout_dir=False,
                maxee=1,
                overwrite=False,
                filetype='fastq'
                ):
    if not outdir:
        outdir = indir
    
    FASTQout_dir,_ = reptools.build_path(True, FASTQout_dir, 'EEfilteredCDR3', outdir)
    FASTAout_dir = False
    
    if not FASTQout_dir and not FASTAout_dir:
        raise ValueError('Please supply one or both of FASTQout_dir and FASTAout_dir to EEfilter_dir()')
    
    for pth in [FASTQout_dir, FASTAout_dir]:
        if pth: reptools.cautious_mkdir(pth,overwrite=overwrite)
    
    filetypes = reptools.select_filetypes(filetype)
    infiles=[fn for fn in os.listdir(indir) if os.path.splitext(fn.lower())[1] in filetypes]
    
    if len(infiles)==0:
        print('No fastq files found.\n')
        return
    
    for fn in infiles:
        FASTQout = reptools.make_unpaired_filepaths(FASTQout_dir, os.path.splitext(fn)[0])
        FASTAout = reptools.make_unpaired_filepaths(FASTAout_dir, os.path.splitext(fn)[0],'fas')
        _ = reptools.EEfilter_file(
                      os.path.join(indir,fn),
                      FASTQout = FASTQout,
                      FASTAout = FASTAout,
                      maxee = maxee
                      )
    
    return(FASTQout_dir)


def EEfilter_file(infile,FASTAout=False,FASTQout=False,maxee=1):
    from reptools import dummy_context_mgr as dummy
    if not FASTAout and not FASTQout:
        raise ValueError('Please supply one or both of FASTAout and FASTQout to EEfilter()')
    with open(infile) as inhandle:
        with open(FASTAout,'w') if FASTAout else dummy() as outfasta_handle:
            with open(FASTQout,'w') if FASTQout else dummy() as outfastq_handle:
                for title,seq,qual in reptools.FASTQparser(inhandle):
                    if reptools.calculate_EE(qual) <= maxee:
                        outfasta_handle.write('>{}\n{}\n'.format(title,seq))
                        outfastq_handle.write('@{}\n{}\n+\n{}\n'.format(title,seq,qual))
    return(reptools.removeemptyfile(FASTQout)) #returns None if the file was empty (and has been removed), else the fn


def calculate_EE(qual):
    """
    Calculate expect error given a Sanger-encoded qual string.
    """
    phreds = [ord(c)-33 for c in qual]
    probs = [10**(-float(Q)/10) for Q in phreds]
    EE = sum(probs)
    return(EE)