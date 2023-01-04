import tempfile
import os
import string
import re
import csv
import reptools
import shutil
import collections
import errno
import shutil

def fileSystemSafeMove(source,dest):
    #monkey patch shutil._copyxattr to avoid permissions error
    #  as workaround for https://bugs.python.org/issue38633
    orig_copyxattr = shutil._copyxattr
    def patched_copyxattr(src, dst, *, follow_symlinks=True):
        try:
            orig_copyxattr(src, dst, follow_symlinks=follow_symlinks)
        except OSError as ex:
            if ex.errno != errno.EACCES: raise
    shutil._copyxattr = patched_copyxattr
    try:
        shutil.move(source,dest)
    except OSError as e:
        if e.errno==errno.EXDEV:
            shutil.copy(source,dest)
            os.remove(source)
        else:
            raise


class FASTQread:
    def __init__(
            self,title,seq,qual
            ):
        self.title = title
        self.seq = seq
        #
        self.qual = qual


class FASTAread:
    def __init__(
            self,title,seq,qual=None
            ):
        self.title = title
        self.seq = seq


class dummy_context_mgr():
    def write(self, data):
        pass # ignore the data
    def __enter__(self): return self
    def __exit__(*x): pass


def createTempfile(suffix=None):
    tf = tempfile.NamedTemporaryFile(delete=False,suffix=suffix)
    tf.close()
    return(tf.name)
    

def checkFiletype(infile,extrachars="",return_empty=False):
    nuc_chars = 'uatcgryswkmbdhvn-' + extrachars
    prot_chars = ''.join(set(codon_table.values())).lower()
    titlefound=False
    with open(infile) as handle:
        for row in handle:
            if len(row.strip())>0:
                if not titlefound:
                    if row.strip()[0]=='>':
                        ftype='fasta'
                        titlefound=True
                    elif row.strip()[0]=='@':
                        ftype='fastq'
                        titlefound=True
                    else:
                        raise IOError(
                                      '{} is not a FASTA or FASTQ file.  Initial character is not @ or >\n'.format(
                                                                                                                 infile
                                                                                                                 )
                                      )
                elif ftype=='fastq' and row.strip().lower()[0] in nuc_chars:
                        return(ftype)
                elif ftype=='fasta' and row.strip().lower()[0] in nuc_chars+prot_chars:
                    return(ftype)
                else:
                    raise IOError(
                                  '{} is not a FASTA or FASTQ file.  '
                                  'Initial character of second (non-empty) line\n'
                                  'is not a valid nucleotide (or valid amino acid, for FASTA only)'.format(infile)
                                  )
        if return_empty and titlefound and len(row.strip()==0):
            return(ftype)
        
        raise IOError('{} is not a FASTA or FASTQ file, or is empty.\n'.format(infile))


def check_filetypes(
    fns
    ):
    """
    Takes either a list of filenames, or a list of lists containing filename pairs.
    Checks that all files are valid FASTQ or FASTA files, and that all files supplied are of the same type.
    Returns the detected filetype ('fastq' or 'fasta').
    """
    #flatten list
    def flatten(lst):
        for x in lst:
            try:
                yield from flatten(x)
            except TypeError:
                yield x
    
    fns = list(flatten(fns))
    ftypes = set()
    for fn in fns:
        ftypes.add(checkFiletype(fn))
        if len(ftypes)>1:
            raise IOError('File structures matching both FASTQ and FASTA formats were present in input.  Pick one, '
                          'and only one.')
    
    return(ftypes.pop())


def reverse_comp(nucl_as_string,ntype='DNA'):
    return(reptools.complement(nucl_as_string,ntype)[::-1])


def complement(nucl_as_string,ntype='DNA'):
    # convert to lower case and trim whitespace at ends, and convert unicode to string
    nucl_as_string = nucl_as_string.lower().strip()
    badchars = dict.fromkeys(string.printable)
    #use translate() to get the characters in all printable characters which aren't in our list
    if ntype.upper()=='DNA':
        for sym in 'Atagcyrswmkvhdbn.--':
            badchars.pop(sym,None) #remove these from the badchars table 
        transtable = str.maketrans('uatcgryswkmbdhvn.-~','Atagcyrswmkvhdbn.--')
        #translates bases, turns ~ to -, and leaves standard symbols alone.
        #If a u is present, it will be translated to a capital A, to highlight the fact that it shouldn't be present in DNA
    elif ntype.upper()=='RNA':
        for sym in 'auAgcyrswmkvhdbn.--':
            badchars.pop(sym,None) #remove these from the badchars table 
        transtable = str.maketrans('uatcgryswkmbdhvn.-~','auAgcyrswmkvhdbn.--')
        #translates bases, turns ~ to -, and leaves standard symbols alone.
        #If a t is present, it will be translated to a capital A, to highlight the fact that it shouldn't be present in RNA
    else:
        raise ValueError('Unknown nucleic acid: %s' % ntype)
    badchars = {ord(k):v for k,v in badchars.items()} #convert badchars table keys to Unicode points
    translated = nucl_as_string.translate(transtable)
    return(translated.translate(badchars))


def select_filetypes(filetype):
    if filetype:
        filetype=filetype.lower()
        if filetype[0] != '.':
            filetype = '.' + filetype
        if filetype == '.fasta':
            filetypes = ['.fasta','.fas','.fsa','.fa']
        elif filetype == '.fastq':
            filetypes = ['.fastq','.fq']
        else:
            filetypes = [filetype]
    else:
        filetypes = ['.fas','.fasta','.fsa','.fastq','.fq','.fa']
    return(filetypes)


def cautious_mkdir(dirname,overwrite=False):
    import errno
    import time
    if overwrite: remove_dir(dirname, recursive=True)
    while True:
        try:
            os.makedirs(dirname)
        except OSError as error:
            if error.errno == errno.EEXIST:
                if os.listdir(dirname):
                    raise error
                else:
                    pass
            else:
                raise error
    
        time.sleep(0.02)
        if os.path.isdir(dirname): break


def ensure_dir(dirname):
    """
    Ensure that a named directory exists; if it does not, attempt to create it.
    @Graham Klyne in https://stackoverflow.com/questions/944536/efficient-way-of-creating-recursive-paths-python
    """
    import errno
    import time
    while True:
        try:
            os.makedirs(dirname)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        time.sleep(0.02)
        if os.path.isdir(dirname): break


def remove_dir(dirname,recursive=False,remove_if_not_empty=True,silent_fail=False):
    #import errno
    import time
    if remove_if_not_empty or recursive:
        from shutil import rmtree as rmtree
    else:
        from os import rmdir as rmtree
    while True:
        try:
            if not os.path.isdir(dirname): break
            if recursive and remove_if_not_empty: #remove everything
                shutil.rmtree(dirname)
                return(True)
            elif recursive and not remove_if_not_empty: #remove an empty tree, but not one with files
                if sum([
                        1 for path, subdirs, files in os.walk(dirname)
                        for name in files if os.path.isfile(os.path.join(path, name))
                        ]) > 0:
                    if silent_fail:
                        return(False)
                    else:
                        raise IOError('{} contains files.'.format(dirname))
                else:
                    os.removedirs(dirname)
                    return(True)
            elif not recursive and remove_if_not_empty:  #remove files, but fail if subdirs present
                if sum([1 for pth in os.listdir(dirname) if os.path.isdir(pth)]) > 0:
                    if silent_fail:
                        return(False)
                    else:
                        raise IOError('{} contains subdirectories.'.format(dirname))
                else:
                    shutil.rmtree(dirname)
                    return(True)
            elif not recursive and not remove_if_not_empty: #remove only an empty directory
                os.rmdir(dirname) #an OS Error will be raised is the dir is not empty
                return(True)
            else:
                raise ValueError('Unexpected condition.  The else fork should not have been called.')
        except OSError as error:
            #if error.errno == errno.ENOENT:
            #    pass
            #else:
                raise error
        time.sleep(0.02)



def fastqcounter(infile):
    """
    Returns the number of unique sequences in a fastq file
    """
    #check if file is derep'd using DerepCheck()
    derep = reptools.DerepCheck(infile)

    n=0
    if derep:
        with open(infile) as fn:
            for title,seq,qual in reptools.FASTQparser(fn):
                n+=reptools.DerepCount(title)
    else:
        with open(infile) as fn:
            for title,seq,qual in reptools.FASTQparser(fn):
                n+=1
    return(n)


def fascounter(infile,derep=None):
    """
    Returns the number of sequences in a fasta file, incorporating Derep information if present
    """
    #check if file is derep'd using DerepCheck()
    derep = reptools.DerepCheck(infile)
    n=0
    if derep:
        with open(infile) as fn:
            for (title,seq) in reptools.FASTAparser(fn):
                n += DerepCount(title)
    else:
        with open(infile) as fn:
            for (title,seq) in reptools.FASTAparser(fn):
                n += 1
    return(n)


def DerepCheck(fn):
    try:
        with open(fn) as inhandle:
            for title,seq in reptools.FASTAparser(inhandle):
                return(is_derepFas(title))
    except ValueError:
        with open(fn) as inhandle:
            for title,seq,qual in reptools.FASTQparser(inhandle):
                return(is_derepFas(title))


def DerepCount(titleline):
    """
    Takes a title from a uSearch dep'd fasta sequence, and returns the sequence count as an integer
    """
    return int((titleline.split('size=')[1]).split(';')[0])


def is_derepFas(title):
    """
    Takes a fasta title line (with or without the initial '>'), and tests to see whether it looks
    like a usearch dereplicated title line - i.e., does it end "size=n;"?
    """
    title = title.strip(';')
    try:
        countclause = title.split(';')[-1] #try to split at ';'
        try:
            if countclause.split('=')[0]=='size': 
            #try to split at '='; and if it splits, is the first part "size"?
                derepped = True
            else:
                derepped = False
        except:
            derepped = False
    except:
        derepped = False
    return(derepped)


def ExtractGene(title,gene):
    subunits = title.split(';')
    outputlist = [subunit.split('=')[1] for subunit in subunits if subunit.split('=')[0]==gene ]
    if len(outputlist)== 1:
        return(outputlist[0])
    elif len(outputlist)== 0:
        return('NA')
    else:
        raise ValueError('More than one gene found by ExtractGene().\noutputlist = %s\ntitle = %s\n' % (outputlist,title))


def CDR3_trans(CDR3nt):
    if len(CDR3nt)%3==0: #in frame
        return(reptools.trans(CDR3nt).upper())
    else: #out of frame
        if len(CDR3nt) < 3:
            return('_')
        forwardcut = int(((len(CDR3nt))/3)/2) * 3 #translate the first ~50% of codons
        if len(CDR3nt)%3==1:  
            reversecut = forwardcut + 1
            #translate the last 50% of codons, rounding down, and skipping two bases inbetween (the "_" translation)
        elif len(CDR3nt)%3==2:
            reversecut = forwardcut + 2
            #translate the last 50% of codons, rounding down, and skipping one base inbetween (the "_" translation)
        else:
            print(('Warning: empty or malformed string supplied to CDR3_trans(): %s \n' % CDR3nt))
            return('')
    aa = reptools.trans(CDR3nt[:forwardcut]).upper()
    aa = aa + '_'
    aa = aa + reptools.trans(CDR3nt[reversecut:]).upper()
    return(aa)


def trans(input_seq,stop='*'):
    """
    translate a string
    this is fast, but note that it doesn't look for invalid characters - it will just translate a codon containing them as "X"
    """
    #import string
    codon_table['TAA'] = stop
    codon_table['TGA'] = stop
    codon_table['TAG'] = stop
    codon_table['UAA'] = stop
    codon_table['UGA'] = stop
    codon_table['UAG'] = stop
    input_seq = input_seq.upper().strip() # convert to upper case and trim whitespace at ends
    #
    if len(input_seq)%3 >0: #is the length divisible by 3 (i.e. in codons)
        input_seq += ''.join(['n']*(3-(len(input_seq)%3))) #add 'n's if it isn't
    #
    split_seq = [input_seq[i:i+3] for i in range(0, len(input_seq), 3)] #split into codons
    #
    outseq = ''.join([codon_table.get(codon, 'X') for codon in split_seq])
    #
    return(outseq)


codon_table = {
                'TTT':'F',
                'TCT':'S',
                'TAT':'Y',
                'TGT':'C',
                'TTC':'F',
                'TCC':'S',
                'TAC':'Y',
                'TGC':'C',
                'TTA':'L',
                'TCA':'S',
                'TAA':'*',
                'TGA':'*',
                'TTG':'L',
                'TCG':'S',
                'TAG':'*',
                'TGG':'W',
                'CTT':'L',
                'CCT':'P',
                'CAT':'H',
                'CGT':'R',
                'CTC':'L',
                'CCC':'P',
                'CAC':'H',
                'CGC':'R',
                'CTA':'L',
                'CCA':'P',
                'CAA':'Q',
                'CGA':'R',
                'CTG':'L',
                'CCG':'P',
                'CAG':'Q',
                'CGG':'R',
                'ATT':'I',
                'ACT':'T',
                'AAT':'N',
                'AGT':'S',
                'ATC':'I',
                'ACC':'T',
                'AAC':'N',
                'AGC':'S',
                'ATA':'I',
                'ACA':'T',
                'AAA':'K',
                'AGA':'R',
                'ATG':'M',
                'ACG':'T',
                'AAG':'K',
                'AGG':'R',
                'GTT':'V',
                'GCT':'A',
                'GAT':'D',
                'GGT':'G',
                'GTC':'V',
                'GCC':'A',
                'GAC':'D',
                'GGC':'G',
                'GTA':'V',
                'GCA':'A',
                'GAA':'E',
                'GGA':'G',
                'GTG':'V',
                'GCG':'A',
                'GAG':'E',
                'GGG':'G',
                'UUU':'F',
                'UCU':'S',
                'UAU':'Y',
                'UGU':'C',
                'UUC':'F',
                'UCC':'S',
                'UAC':'Y',
                'UGC':'C',
                'UUA':'L',
                'UCA':'S',
                'UAA':'*',
                'UGA':'*',
                'UUG':'L',
                'UCG':'S',
                'UAG':'*',
                'UGG':'W',
                'CUU':'L',
                'CCU':'P',
                'CAU':'H',
                'CGU':'R',
                'CUC':'L',
                'CUA':'L',
                'CUG':'L',
                'AUU':'I',
                'ACU':'T',
                'AAU':'N',
                'AGU':'S',
                'AUC':'I',
                'AUA':'I',
                'AUG':'M',
                'GUU':'V',
                'GCU':'A',
                'GAU':'D',
                'GGU':'G',
                'GUC':'V',
                'GUA':'V',
                'GUG':'V'
                }


def FASTAparser(in_handle):
    """
    A simple and fast parser for FASTA files.
    This parser removes whitespace characters from the end of lines, but does no other checking, so non-FASTA compliant characters will not be detected or removed.
    Takes an open file object.
    Returns a (title,sequence) tuple.
    """
    in_handle.read(0) #will throw an error if in_handle isn't a file-like object
    recordfound = False
    for row in in_handle:
        row=row.strip()
        if len(row)!=0: #skip empty lines
            if row[0]=='>':
                if not recordfound:
                    recordfound=True
                else:
                    yield(tuple([title,''.join(seq.split())]))
                title=row[1:]
                seq=''
            else:
                if recordfound:
                    seq=seq+row
                else:
                    raise ValueError('Not a fasta file (first non-empty line must start with >).')
    if recordfound:
        yield(tuple([title,''.join(seq.split())]))
    else:
        raise ValueError('Empty fasta file or not a fasta file')


def FASTQparser(in_handle):
    """
    A simple and fast parser for FASTQ files.
    This parser removes whitespace characters from the end of lines, but does no other checking, so non-FASTQ compliant characters will not be detected or removed.
    Takes an open file object.
    Returns a (title,sequence,quality) tuple.
    """
    in_handle.read(0) #will throw an error if in_handle isn't a file-like object
    stage = 0
    for row in in_handle:
        row=row.strip()
        if stage == 0: #nothing yet found
            if row[0]=='@':
                title=row[1:]
                seq = ''
                qual = ''
                stage = 1
        elif stage == 1: #reading sequence
            if row.strip() == '+' or row.strip()=='+'+title:
                seq = ''.join(seq.split())
                stage = 2
            else:
                seq=seq+row
        elif stage == 2: #reading qual
            if row[0]=='@' and len(qual)==len(seq):
                yield(tuple([title,seq,qual]))
                title=row[1:]
                seq = ''
                qual = ''
                stage = 1
            else: 
                qual = qual + ''.join(row.split())
            if len(qual)>len(seq):
                raise ValueError('Sequence/quality string length mismatch.')
        else:
            raise ValueError('Internal error: stage not in [0,1,2]')
    if stage==0:
        raise ValueError('Empty fastq file or not a fastq file in handle: {}'.format(str(in_handle)))
    if stage==2:
        if len(qual)==len(seq):
            yield(tuple([title,seq,qual]))
        else:
            raise ValueError('Sequence/quality string length mismatch.')


def nt_check(seq,specialchars=''):
    nt_chars = 'AGCTURYKMSWBDHVN-'
    nt_chars= re.compile(r'[^'+nt_chars+specialchars+']')
    if bool(nt_chars.search(seq.upper())):
        return(False)
    else:
        return(True)

        
def remove_file_if_present(fn):
    import errno
    import time
    while True:
        try:
            os.remove(fn)
        except OSError as error:
            if error.errno == errno.ENOENT:
                pass
            else:
                raise error
        time.sleep(0.1)
        if not os.path.exists(fn): break


def remove_file_if_empty(fn):
    import errno
    try:
        if os.path.getsize(fn)==0:
            os.remove(fn)
            return(True)
    except OSError as error:
        if error.errno == errno.ENOENT:
            pass
        else:
            raise error
    return(False)


def combinefiles(infiles,outfile):
    with open(outfile,'w') as out_handle:
        for infile in infiles:
            with open(infile,'r') as in_handle:
                out_handle.write(in_handle.read())


def fastq2fasta(infile, outfile="",trimstart=0,overwrite=False):
    """
    takes infile and outfile (file names)
    """
    if outfile=="":
        if infile[-5:]=='fastq':
            outfile=infile[:-5]+'fas'
        else:
            outfile=infile+'.fas'
    if os.path.isfile(outfile) and not overwrite:
        raise IOError('Output file (%s) already exists' % (outfile))
    for (title, seq, qual) in reptools.FASTQparser(open(infile,'r')):
        with open(outfile,'a') as f:
            f.write('>%s\n%s\n' % (title,seq[trimstart:]))
    return outfile


def clean_up(file_list):
    import time
    time.sleep(0.1)
    for fn in file_list:
        if fn: #so that False is OK in the file_list
            reptools.remove_file_if_present(fn)


def clean_up_dirs(dir_list):
    import time
    time.sleep(0.1)
    for dir in dir_list:
        if dir: #so that False is OK in the dir_list
            try:
                shutil.rmtree(dir)
            except FileNotFoundError:
                pass


def removeemptyfile(fn):
    if fn and os.path.getsize(fn)==0:
        os.remove(fn)
        return(None)
    else:
        return(fn)


def removeemptyfilepair(fn1,fn2):
    if (fn1 and os.path.getsize(fn1)==0) or (fn2 and os.path.getsize(fn2)==0):
        os.remove(fn1)
        os.remove(fn2)
        return((None,None))
    else:
        return((fn1,fn2))


def fill_defaults(var, keys, values, ordered=False):
    if ordered:
        outvar = collections.OrderedDict()
    else:
        outvar = {}
    if not var: var = {}
    for k,v in zip(
            keys,
            values
            ):
        if k in var:
            outvar[k] = var[k]
        else:
            outvar[k] = v
    return(outvar)


def build_path(selected,selected_dir,default_dir,outdir):
    """
    selected should be True or False
    If selected==False (or None), return False
    If selected is True:
        if selected_dir is set:
            return it directly if it is an absolute path, or append it to outdir if it is a relative path
        if selected dir is not set (or is "None"/"none"):
            if it is "False" or "false" or False, return the path to a tempdir
            if it is None, return the supplied default (appended to outdir)
    Always returns a tuple.  The first element is the directory path, or False. The second element is None, unless a 
    temporary directory has been produced, in which case both elements contain the path to the temporary directory.
    """
    if not selected:
        return(False,None)
    
    if isinstance(selected_dir, str): #if selected_dir, use it as the basis for the return value
        if selected_dir.lower() == 'false': #unless its a string version of False, in which case set to False
            selected_dir = False
        elif selected_dir.lower() == 'none': #or a string version of None, in which case set to None
            selected_dir = None
        else:
            if os.path.isabs(selected_dir):
                return(selected_dir,None)
            else:
                return(os.path.join(outdir,selected_dir),None)
    
    if selected_dir is None: #is selected_dir is None, return the default
        return(os.path.join(outdir,default_dir),None)
    elif selected_dir is False: #is selected_dir is False, return the path to a temporary directory
        td = tempfile.mkdtemp()
        return(td,td)
    else:
        raise ValueError('selected_dir must be False, None, or a string')


def build_fn(selected_fn,default_fn,outdir,maketemp=True):
    """
    If selected is not False:
        if selected_fn is set:
            return it directly if it is an absolute path, or append it to outdir if it is a relative path
        if selected_fn dir is not set (or is "None"/"none"):
            if it is "None" or "none" or [N/B/ string, not type] or None [type], return the supplied default_fn
            (appended to outdir)
            if it is False, return the path to a tempfile
    """
    
    if isinstance(selected_fn, str): #if selected_fn, use it as the basis for the return value
        if selected_fn.lower() == 'none': #unless its a string version of None, in which case set to NoneType
            selected_fn = None
        elif selected_fn.lower() == 'false': #or  a string version of None, in which case set to False
            selected_fn = False
        else:
            if os.path.isabs(selected_fn):
                return(selected_fn)
            else:
                return(os.path.join(outdir,selected_fn))
    
    if selected_fn is None: #is selected_fn is None, return the default
        return(os.path.join(outdir,default_fn))
    elif selected_fn is False: #is selected_fn is False, return the path to a temporary directory
        if maketemp:
            return(reptools.createTempfile())
        else:
            return(False)
    else:
        raise ValueError('selected_fn must be False, None, or a string')


def make_paired_filepaths(folder,stem,pairsuffixes):
    if folder:
        fn1 = os.path.join(folder,'{}{}.fastq'.format(stem,pairsuffixes[0]))
        fn2 = os.path.join(folder,'{}{}.fastq'.format(stem,pairsuffixes[1]))
        return(fn1,fn2)
    else:
        return(False,False)


def make_unpaired_filepaths(folder,stem,ext='fastq'):
    if folder:
        return(os.path.join(folder,'{}.{}'.format(stem,ext)))
    else:
        return(False)

