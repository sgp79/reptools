import subprocess
import os
import tempfile
import multiprocessing
import csv
import string
import subprocess
import shutil
import time
import math
import reptools

def sswCall(
            infile,
            db,
            hitsOut=False,
            sswPath = 'ssw_test',
            minscore=24,
            **kwargs
    ):
    call_list = [sswPath]
    call_list = call_list + ['-f', str(minscore)]
    call_list = call_list + ['-r']
    call_list = call_list + ['-s']
    call_list = call_list + ['-c']
    call_list = call_list + [db]
    call_list = call_list + [infile]
    if hitsOut:
        with open(hitsOut,'w') as hitshandle, open(os.devnull,'w') as errhandle:
            subprocess.call(call_list, stdout=hitshandle, stderr=errhandle)
    else:
        with open(os.devnull,'w') as hitshandle, open(os.devnull,'w') as errhandle:
            subprocess.call(call_list, stdout=hitshandle, stderr=errhandle)


def nhmmerCall(
    infile,
    db,
    hitsOut = False,
    threads = False,
    evalue=0.001,
    makehmmerdbPath = 'makehmmerdb',
    nhmmerPath = 'nhmmer',
    matchedfq = False,
    notmatchedfq = False,
    **kwargs
    ):
    if not threads: 
        threads = multiprocessing.cpu_count()-1
    else:
        if threads > multiprocessing.cpu_count()-1:
            threads = multiprocessing.cpu_count()-1

    todelete = []
    
    #check file type - if fasta, make sure that no fastq output options are selected
    filetype = reptools.checkFiletype(infile)
    if filetype=='fasta':
        if matchedfq or notmatchedfq:
            raise IOError(
                'Cannot output fastq files when the input is fasta.  '
                'Set matchedfq and notmatchedfq to False.\n')
    #check file type - import appropriate Parser
    if filetype == 'fastq':
        from reptools import FASTQparser as Parser
    else:
        from reptools import FASTAparser as Parser
    #check file type - if fastq, make a fasta tempfile
    if filetype=='fastq':
        tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
        tf.close()
        reptools.fastq2fasta(infile,tf.name,overwrite=True)
        workinginfile = tf.name
        todelete.append(tf.name)
    else: 
        workinginfile = infile
    
    #make a hmmer database file
    #first, copy the input db to a tempfile, as this is the easiest way to ensure that there are no spaces in the path
    #t_in = tempfile.NamedTemporaryFile(suffix='.fas')
    #t_in.close()
    #shutil.copy(db,t_in.name)
    hmmer_db =  tempfile.NamedTemporaryFile(delete='False')
    hmmer_db.close()
    todelete.append(hmmer_db.name)
    
    call_list = [makehmmerdbPath]
    call_list = call_list + ['--informat','fasta'] 
    call_list = call_list + ['--sa_freq','4']
    call_list = call_list + [db]
    call_list = call_list + [hmmer_db.name]
    p = subprocess.check_call(call_list)
    
    call_list = [nhmmerPath]
    if hitsOut: call_list = call_list + ['--tblout', hitsOut]
    call_list = call_list + ['--notextw']
    call_list = call_list + ['--pextend','0.1']
    call_list = call_list + ['-E',str(evalue)]
    call_list = call_list + ['--qfasta']
    call_list = call_list + ['--tformat','fmindex']
    call_list = call_list + ['--cpu',str(threads)]
    call_list = call_list + ['--dna']
    call_list = call_list + [workinginfile]
    call_list = call_list + [hmmer_db.name]
    
    p = subprocess.call(call_list)
    
    reptools.clean_up(todelete)


def parasailCall(
            infile,
            db,
            id=0.95,
            hitsOut=False,
            **kwargs
    ):
    import parasail
    if threads>multiprocessing.cpu_count(): threads = multiprocessing.cpu_count()
    
    todelete = []


def tblastxCall(infile,
                db,
                id=0.95,
                matched=False,
                matchedfq=False,
                hitsOut=False,
                strand='both',
                mincols=False,
                blastPath='tblastx',
                makeblastdbPath = 'makeblastdb',
                notmatched=False,
                notmatchedfq=False,
                return_counts=False,
                evalue=0.001,
                verbose=True,
                threads=10, 
                wordlength=7,
                gapopen=3,
                gapextend=3,
                reward=1,
                penalty=-3,
                blastdb_version=5,
                **kwargs
                ):
    if threads>multiprocessing.cpu_count(): threads = multiprocessing.cpu_count()
    #from sys import platform
    
    todelete = []
    
    #check file type - if fastq, make a fasta version
    filetype = reptools.checkFiletype(infile)
    if filetype=='fasta':
        if matchedfq or notmatchedfq:
            raise IOError(
                'Cannot output fastq files when the input is fasta.  '
                'Set matchedfq and notmatchedfq to False.\n')
    #check file type - import appropriate Parser
    if filetype == 'fastq':
        from reptools import FASTQparser as Parser
    else:
        from reptools import FASTAparser as Parser
    #check file type - if fastq, make a fasta tempfile
    if filetype=='fastq':
        tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
        tf.close()
        reptools.fastq2fasta(infile,tf.name,overwrite=True)
        workinginfile = tf.name
        todelete.append(tf.name)
    else: 
        workinginfile = infile
    #

    #make a BLAST formatted db file
    #first, copy the input db to a tempfile, as this is the easiest way to ensure that there are no spaces in the path
    blast_db,td = reptools.makeblastdb(db,makeblastdbPath,blastdb_version=blastdb_version)
    
    #if hitsOut isn't set, make a temporary file to store results
    if not hitsOut:
        tf = tempfile.NamedTemporaryFile(delete=False,suffix='.b6')
        tf.close()
        hitsOut = tf.name
        todelete.append(tf.name)    
    
    #call BLAST+
    print(('BLAST+ call with %s\n' % infile))
    call_list = [blastPath]
    call_list.extend(['-db',blast_db])
    call_list.extend(['-query',workinginfile])
    call_list.extend(['-out',hitsOut])
    call_list.extend(['-evalue',str(evalue)])
    call_list.extend(['-num_threads',str(threads)])
    call_list.extend(['-outfmt','6']) #b6 output
    #call_list.extend(['-word_size',str(wordlength)])
    #call_list.extend(['-perc_identity',str(id*100)])
    #call_list.extend(['-gapopen',str(gapopen)])
    #call_list.extend(['-gapextend',str(gapextend)])
    #call_list.extend(['-reward',str(reward)])
    #call_list.extend(['-penalty',str(penalty)])
    if strand.lower() == 'plus' or strand.lower() == 'forward' or strand.lower() == 'f':
        call_list.extend(['-strand','plus'])
    elif strand.lower() == 'minus' or strand.lower() == 'reverse' or strand.lower() == 'r':
        call_list.extend(['-strand','minus'])
    else:
        call_list.extend(['-strand','both'])
    call_output = subprocess.run(call_list,text=True)
    if verbose:
        print(call_output.stdout)
    if call_output.stderr and len(call_output.stderr)>0:
        if 'warning' in call_output.stderr.lower():
            print(('\nBLAST+ reports: %s\n\n' % call_output))
        else:
            raise IOError('\nAfter call:\n%s\n\nBLAST+ reports: %s\n\n' % (call_list,call_output.stderr))
            
    reptools.clean_up(todelete)
    shutil.rmtree(td)


def makeblastdb(
                db,
                makeblastdbPath,
                blastdb_version=5,
                verbose=True
                ):    
    t_in = reptools.createTempfile(suffix='.fas')
    shutil.copy(db,t_in)
    td = tempfile.mkdtemp()
    blast_db = os.path.join(td,os.path.split(db)[1])
    call_list = [makeblastdbPath]
    call_list = call_list + ['-in',t_in] 
    call_list = call_list + ['-dbtype','nucl']
    call_list = call_list + ['-parse_seqids']
    call_list = call_list + ['-out',blast_db]
    call_list = call_list + ['-blastdb_version',str(blastdb_version)]
    
    call_output = subprocess.run(call_list, text=True, capture_output=True)
    
    if verbose:
        if call_output.stdout: print(call_output.stdout)
    if call_output.stderr and len(call_output.stderr)>0 :
        if 'error' in call_output.stderr.lower():
            raise IOError('\nAfter call:\n%s\n\nnmakeblastdb reports: %s\n\n' % (call_list,call_output.stderr))
        elif 'warning' in call_output.stderr.lower():
            print(('\nmakeblastdb reports: %s\n\n' % call_output.stderr))
        else:
            raise IOError('\nAfter call:\n%s\n\nnmakeblastdb reports: %s\n\n' % (call_list,call_output.stderr))
    
    os.remove(t_in)
    time.sleep(0.2)
    return(blast_db,td)


def blastCall(infile,
                db,
                id=False,
                matched=False,
                matchedfq=False,
                hitsOut=False,
                strand='both',
                mincols=False,
                blastPath='blastn',
                makeblastdbPath = 'makeblastdb',
                notmatched=False,
                notmatchedfq=False,
                return_counts=False,
                evalue=0.001,
                verbose=True,
                threads=False,
                wordlength=7,
                gapopen=3,
                gapextend=3,
                reward=1,
                penalty=-3,
                max_target_seqs=20,
                blastdb_version=5,
                **kwargs
                ):
    """
    Takes a fasta or fastq file and a fasta database, and calls BLAST.
    """
    if not threads: threads = multiprocessing.cpu_count()
    if threads > multiprocessing.cpu_count(): threads = multiprocessing.cpu_count()
    #from sys import platform
    #if minLength has been supplied, but not mincols, fill mincols from minLength
    if not mincols and 'minLength' in kwargs: 
        mincols=kwargs['minLength']
    todelete = []
    #check file type - if fasta, make sure that no fastq output options are selected
    filetype = reptools.checkFiletype(infile)
    if filetype=='fasta':
        if matchedfq or notmatchedfq:
            raise IOError(
                'Cannot output fastq files when the input is fasta.  '
                'Set matchedfq and notmatchedfq to False.\n')
    #check file type - import appropriate Parser
    if filetype == 'fastq':
        from reptools import FASTQparser as Parser
    else:
        from reptools import FASTAparser as Parser
    #check file type - if fastq, make a fasta tempfile
    if filetype=='fastq':
        tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
        tf.close()
        reptools.fastq2fasta(infile,tf.name,overwrite=True)
        workinginfile = tf.name
        todelete.append(tf.name)
    else: 
        workinginfile = infile
    
    #make a BLAST formatted db file
    blast_db,td = reptools.makeblastdb(db,makeblastdbPath,blastdb_version=blastdb_version)
    
    #if hitsOut isn't set, make a temporary file to store results
    if not hitsOut:
        tf = tempfile.NamedTemporaryFile(delete=False,suffix='.b6')
        tf.close()
        hitsOut = tf.name
        todelete.append(tf.name)    
    
    #call BLAST+
    print(('BLAST+ call with %s\n' % infile))
    call_list = [blastPath]
    call_list.extend(['-db',blast_db])
    call_list.extend(['-query',workinginfile])
    call_list.extend(['-out',hitsOut])
    call_list.extend(['-evalue',str(evalue)])
    call_list.extend(['-num_threads',str(threads)])
    call_list.extend(['-outfmt','6']) #b6 output
    call_list.extend(['-word_size',str(wordlength)])
    call_list.extend(['-gapopen',str(gapopen)])
    call_list.extend(['-gapextend',str(gapextend)])
    call_list.extend(['-reward',str(reward)])
    call_list.extend(['-penalty',str(penalty)])
    if max_target_seqs: call_list.extend(['-max_target_seqs',str(max_target_seqs)])
    if id: call_list.extend(['-perc_identity',str(id*100)])
    if strand.lower() == 'plus' or strand.lower() == 'forward' or strand.lower() == 'f':
        call_list.extend(['-strand','plus'])
    elif strand.lower() == 'minus' or strand.lower() == 'reverse' or strand.lower() == 'r':
        call_list.extend(['-strand','minus'])
    else:
        call_list.extend(['-strand','both'])
    call_output = subprocess.run(call_list,text=True,capture_output=True)
    if verbose:
        print(call_output.stdout)
        print(call_output.stderr)
    if call_output.stderr and len(call_output.stderr)>0 :
        if 'warning' in call_output.stderr.lower():
            print(('\nBLAST+ reports: %s\n\n' % call_output.stderr))
        else:
            raise IOError('\nAfter call:\n%s\n\nBLAST+ reports: %s\n\n' % (call_list,call_output.stderr))
    
    if matched or matchedfq or notmatched or notmatchedfq:
        hits = set([q.id for q in reptools.hitsParser(hitsOut,filetype='blast')])
        
        with open(infile) as in_handle:
            with open(matched,'w') if matched else reptools.dummy_context_mgr() as match_target,\
             open(matchedfq,'w') if matchedfq else reptools.dummy_context_mgr() as matchfq_target,\
             open(notmatched,'w') if notmatched else reptools.dummy_context_mgr() as unmatch_target,\
             open(notmatchedfq,'w') if notmatchedfq else reptools.dummy_context_mgr() as unmatchfq_target:
                for seqtuple in Parser(in_handle):
                    if seqtuple[0].split(' ')[0] in hits:
                        match_target.write('>%s\n%s\n' % seqtuple[0:2])
                        matchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple[0:3])
                    else:
                        unmatch_target.write('>%s\n%s\n' % seqtuple[0:2])
                        unmatchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple[0:3])
    #raise ValueError('Infile={}\nBLASTout={}\nBLASTdb={}'.format(workinginfile,blast_db,hitsOut))
    reptools.clean_up(todelete)
    shutil.rmtree(td)
    if return_counts:
        if matched:
            return(fascounter(matched))
        if matchedfq:
            return(reptools.fastqcounter(matchedfq))


def uSearchcall(infile,
                db,
                id=0.9,
                matched=False,
                matchedfq=False,
                hitsOut=False,
                alnout=False,
                fastapairs=False,
                strand='both',
                evalue=0.001,
                top_hits_only=True,
                top_hit_only=False,
                accel=False,
                maxaccepts='0',
                maxrejects='0',
                maxhits=False,
                mincols=False,
                usearchPath='usearch',
                return_counts=False,
                verbose = True,
                output_no_hits = False,
                notmatched=False,
                notmatchedfq=False,
                targetcov=None,
                minhsp=False,
                wordlength=False,
                fulldp=False,
                type='local',
                **kwargs
                ):
    #if minLength has been supplied, but not mincols, fill mincols from minLength
    if not mincols and minLength in kwargs: 
        mincols=minLength
    
    print(('uSearchcall with %s\n' % infile))
    call_list = [usearchPath]
    if type =='local': call_list=call_list+['-usearch_local',infile]
    if type =='global':    call_list=call_list+['-usearch_global',infile]
    if type == 'ublast': 
        call_list=call_list+['-ublast',infile]
        if accel: call_list=call_list+['-accel',str(accel)]
        if not evalue: call_list=call_list+['-evalue',str(0.05)]
    call_list=call_list+['-db',db]
    call_list=call_list+['-strand',strand]
    call_list=call_list+['-id',str(id)]
    call_list=call_list+['-qmask','none']
    call_list=call_list+['-dbmask','none']
    if alnout: call_list=call_list+['-alnout',alnout]
    if hitsOut:
        call_list=call_list+['-userout',hitsOut]
        call_list=call_list + [
                    '-userfields',
                    'query+target+ql+qlo+qhi+tlo+thi+alnlen+ids+id+qstrand+tstrand+bits+evalue'
                    ]
    if matchedfq: call_list=call_list+['-matchedfq',matchedfq]
    if matched: call_list=call_list+['-matched',matched]
    if top_hits_only: call_list=call_list+['-top_hits_only']
    if top_hit_only: call_list=call_list+['-top_hit_only']
    if maxaccepts is not False: call_list=call_list+['-maxaccepts',str(maxaccepts)]
    if maxrejects is not False: call_list=call_list+['-maxrejects',str(maxrejects)]
    if maxhits is not False: call_list=call_list+['-maxhits',str(maxhits)]
    if fastapairs: call_list=call_list+['-fastapairs',fastapairs]
    if evalue: call_list=call_list+['-evalue',str(evalue)]
    if mincols: call_list=call_list+['-mincols',str(mincols)]
    if output_no_hits: call_list=call_list+['-output_no_hits']
    if notmatched: call_list=call_list+['-notmatched',notmatched]
    if notmatchedfq: call_list=call_list+['-notmatchedfq',notmatchedfq]
    if targetcov: call_list=call_list+['-target_cov',str(targetcov)]
    if minhsp:  call_list=call_list+['-minhsp',str(minhsp)]
    if wordlength:  call_list=call_list+['-wordlength',str(wordlength)]
    if fulldp:  call_list=call_list+['-fulldp']
    subprocess.crun(call_list,stdout=output,stderr=err)
    if verbose:
        print(output)
        print(err)
    if return_counts:
        if matchedfq:
            return(reptools.fastqcounter(matchedfq))
        elif matched:
            return(reptools.fastacounter(matched))
        else:
            pass


def vsearchCall(infile,
                db,
                id=0.9,
                matched=False,
                matchedfq=False,
                hitsOut=False,
                alnout=False,
                fastapairs=False,
                strand='both',
                top_hits_only=True,
                maxaccepts='0',
                maxrejects='0',
                maxhits=False,
                mincols=False,
                vsearchPath='vsearch',
                return_counts=False,
                verbose = True,
                output_no_hits = False,
                notmatched=False,
                notmatchedfq=False,
                targetcov=None,
                minhsp=False,
                wordlength=False,
                fulldp=False,
                type='global',
                threads=False,
                gapopen = '20',
                endgapopen = '2', 
                gapext = '2',
                endgapext = '1',
                minwordmatches = False,
                **kwargs
                ):
    """
    Takes a fasta or fastq file and a fasta database, and calls the vsearch aligner.
    Outputs any combination of u14 file (of hits), and fasta/fastq files of hits or misses.
    """
    if threads:
        if threads>multiprocessing.cpu_count(): threads = multiprocessing.cpu_count()
    
    #if minLength has been supplied, but not mincols, fill mincols from minLength
    if not mincols and minLength in kwargs: 
        mincols=minLength
    
    gapopenstring = '{}I/{}E'.format(gapopen,endgapopen)
    gapextstring = '{}I/{}E'.format(gapext,endgapext)
    
    todelete = []

    if type.lower() !='global': 
        raise ValueError(
                'As of v2.10.4, vsearch only performs global searches: local and ublast have not been implemented.\n'
                        )
    print(('vSearchcall with %s\n' % infile))
    
    #if input is fastq, need to make a temporary fasta version
    filetype = reptools.checkFiletype(infile)
    if filetype=='fasta':
        if matchedfq or notmatchedfq:
            raise IOError(
                'Cannot output fastq files when the input is fasta.  '
                'Set matchedfq and notmatchedfq to False.\n')
        workinginfile=infile
    else:
        tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
        tf.close()
        todelete.append(tf.name)
        reptools.fastq2fasta(infile,tf.name,overwrite=True)
        workinginfile = tf.name
    
    #for fastq output, we need a hits file to parse
    if not hitsOut and filetype == 'fastq' and (matchedfq or notmatchedfq): 
        tf = tempfile.NamedTemporaryFile(delete=False,suffix='.u14')
        tf.close()
        hitsOut = tf.name
        todelete.append(tf.name)
    
    call_list = [vsearchPath]
    call_list=call_list+['--usearch_global',workinginfile]
    
    call_list=call_list+['--db',db]
    call_list=call_list+['--strand',strand]
    call_list=call_list+['--id',str(id)]
    call_list=call_list+['--qmask','none']
    call_list=call_list+['--dbmask','none']
    call_list=call_list+['--fasta_width','0']
    call_list=call_list+['--gapopen',gapopenstring]
    call_list=call_list+['--gapext',gapextstring]
    if alnout: call_list=call_list+['--alnout',alnout]
    if hitsOut:
        call_list=call_list+['--userout',hitsOut]
        call_list=call_list + [
                    '--userfields',
                    'query+target+ql+qlo+qhi+tlo+thi+alnlen+ids+id+qstrand+tstrand+bits+evalue'
                    ]
    if threads: call_list=call_list+['--threads',str(threads)]
    if matched: call_list=call_list+['--matched',matched]
    if fastapairs: call_list=call_list+['--fastapairs',fastapairs]
    if top_hits_only: call_list=call_list+['--top_hits_only']
    if maxaccepts is not False: call_list=call_list+['--maxaccepts',str(maxaccepts)]
    if maxrejects is not False: call_list=call_list+['--maxrejects',str(maxrejects)]
    if maxhits is not False: call_list=call_list+['--maxhits',str(maxhits)]
    if mincols: call_list=call_list+['--mincols',str(mincols)]
    if output_no_hits: call_list=call_list+['--output_no_hits']
    if notmatched: call_list=call_list+['--notmatched',notmatched]
    if targetcov: call_list=call_list+['--target_cov',str(targetcov)]
    if minhsp:  call_list=call_list+['--minhsp',str(minhsp)]
    if wordlength:  call_list=call_list+['--wordlength',str(wordlength)]
    if minwordmatches:  call_list=call_list+['--minwordmatches',str(minwordmatches)]
    call_output = subprocess.run(call_list,text=True)
    if verbose:
        print(call_output.stdout)
        print(call_output.stderr)
    #fastq output,if required (by parsing hitsfile and iterating through input file)
    if matchedfq or notmatchedfq:
        hits = set([q.id for q in reptools.hitsParser(hitsOut,filetype='stellar',evalue=evalue)])
        
        with open(infile) as in_handle:
            with open(matchedfq,'w') if matchedfq else reptools.dummy_context_mgr() as matchfq_target,\
             open(notmatchedfq,'w') if notmatchedfq else reptools.dummy_context_mgr() as unmatchfq_target:
                for seqtuple in FASTQparser(in_handle):
                    if seqtuple[0] in hits:
                        matchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple)
                    else:
                        unmatchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple)
    
    reptools.clean_up(todelete)
    
    if return_counts:
        if matchedfq:
            return(reptools.fastqcounter(matchedfq))
        elif matched:
            return(reptools.fastacounter(matched))
        else:
            pass


def swipeCall(infile,
                db,
                id=0.95,
                matched=False,
                matchedfq=False,
                hitsOut=False,
                strand='both',
                mincols=False,
                swipePath='swipe',
                makeblastdbPath = 'makeblastdb',
                notmatched=False,
                notmatchedfq=False,
                return_counts=False,
                evalue=0.001,
                verbose=True,
                threads=False,
                gapopen=2,
                gapextend=2,
                reward=1,
                penalty=-3,
                max_target_seqs=20,
                **kwargs
                ):
    """
    Takes a fasta or fastq file and a fasta database, and calls the SWIPE aligner.
    """
    
    if not threads: threads = multiprocessing.cpu_count()
    if threads > multiprocessing.cpu_count(): threads = multiprocessing.cpu_count()

    #if minLength has been supplied, but not mincols, fill mincols from minLength
    if not mincols and 'minLength' in kwargs: 
        mincols=kwargs['minLength']
    todelete = []
    #check file type - if fasta, make sure that no fastq output options are selected
    filetype = reptools.checkFiletype(infile)
    if filetype=='fasta':
        if matchedfq or notmatchedfq:
            raise IOError(
                'Cannot output fastq files when the input is fasta.  '
                'Set matchedfq and notmatchedfq to False.\n')
    #check file type - import appropriate Parser
    if filetype == 'fastq':
        from reptools import FASTQparser as Parser
    else:
        from reptools import FASTAparser as Parser
    #check file type - if fastq, make a fasta tempfile
    if filetype=='fastq':
        tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
        tf.close()
        reptools.fastq2fasta(infile,tf.name,overwrite=True)
        workinginfile = tf.name
        todelete.append(tf.name)
    else:
        workinginfile = infile
    #
    #make a BLAST formatted db file
    blast_db,td = reptools.makeblastdb(db,makeblastdbPath,blastdb_version=4)
    
    #while not os.path.exists(blast_db): time.sleep(0.2) #required for SWIPE to see the file
    
    #if hitsOut isn't set, make a temporary file to store results
    if not hitsOut:
        tf = tempfile.NamedTemporaryFile(delete=False,suffix='.tsv')
        tf.close()
        hitsOut = tf.name
        todelete.append(tf.name)    
    
    #call SWIPE
    print(('SWIPE call with %s\n' % infile))
    call_list = [swipePath]
    call_list.extend(['='.join(['--db',blast_db])])
    call_list.extend(['='.join(['--query',workinginfile])])
    call_list.extend(['--symtype=0'])
    call_list.extend(['--outfmt=8'])
    call_list.extend(['='.join(['--out',hitsOut])])
    call_list.extend(['='.join(['--num_threads',str(threads)])])
    call_list.extend(['='.join(['--penalty',str(penalty)])])
    call_list.extend(['='.join(['--reward',str(reward)])])
    call_list.extend(['='.join(['--gapopen',str(gapopen)])])
    call_list.extend(['='.join(['--gapextend',str(gapextend)])])
    if max_target_seqs: call_list.extend(['='.join(['--num_alignments',str(max_target_seqs)])])
    if evalue: call_list.extend(['='.join(['--evalue',str(evalue)])])
    if strand.lower() == 'plus' or strand.lower() == 'forward' or strand.lower() == 'f':
        call_list.extend(['--strand=1'])
    if strand.lower() == 'minus' or strand.lower() == 'reverse' or strand.lower() == 'r':
        call_list.extend(['--strand=2'])
    #print(call_list)
    #print(td)
    call_output = subprocess.run(call_list, text=True, capture_output=True)
    if verbose:
        if call_output.stdout: print(call_output.stdout)
    if call_output.stderr:
        raise IOError('\nSWIPE reports: %s\n' % call_output.stderr)
    
    #select and output hits
    hits = set([q.id for q in reptools.hitsParser(hitsOut,filetype='swipe',mincols=mincols)])
    
    with open(infile) as in_handle:
        with open(matched,'w') if matched else reptools.dummy_context_mgr() as match_target,\
         open(matchedfq,'w') if matchedfq else reptools.dummy_context_mgr() as matchfq_target,\
         open(notmatched,'w') if notmatched else reptools.dummy_context_mgr() as unmatch_target,\
         open(notmatchedfq,'w') if notmatchedfq else reptools.dummy_context_mgr() as unmatchfq_target:
            for seqtuple in Parser(in_handle):
                if seqtuple[0].split(' ')[0] in hits:
                    match_target.write('>%s\n%s\n' % seqtuple[0:2])
                    matchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple[0:3])
                else:
                    unmatch_target.write('>%s\n%s\n' % seqtuple[0:2])
                    unmatchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple[0:3])
    
    reptools.clean_up(todelete)
    shutil.rmtree(td)
    if return_counts:
        if matched:
            return(reptools.fastqcounter(matched))
        if matchedfq:
            return(reptools.fastacounter(matchedfq))


def stellarCall(infile,
                db,
                id=0.95,
                matched=False,
                matchedfq=False,
                hitsOut=False,
                alnout=None,
                fastapairs=False,
                strand='both',
                mincols=False,
                stellarPath='stellar',
                notmatched=False,
                notmatchedfq=False,
                xDrop=3,
                return_counts=False,
                alphabet='dna5',
                numMatches=10,
                evalue=0.001,
                verbose=True,
                **kwargs
                ):
    """
    Takes a fasta or fastq file and a fasta database, and calls the STELLAR aligner.
    Outputs any combination of gff file (of hits), and fasta/fastq files of hits or misses.
    Note that STELLAR delimits title lines at spaces, and the portion before the first spaces
    must be unique.
    """
    
    todelete = []
    err = 1-id #set err from id
    #check file type - if fasta, make sure that no fastq output options are selected
    filetype = reptools.checkFiletype(infile)
    if filetype=='fasta':
        if matchedfq or notmatchedfq:
            raise IOError(
                'Cannot output fastq files when the input is fasta.  '
                'Set matchedfq and notmatchedfq to False.\n')
    #check file type - import appropriate Parser
    if filetype == 'fastq':
        from reptools import FASTQparser as Parser
    else:
        from reptools import FASTAparser as Parser
    #check file type - if fastq, make a fasta tempfile
    if filetype=='fastq':
        tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
        tf.close()
        reptools.fastq2fasta(infile,tf.name,overwrite=True)
        workinginfile = tf.name
        todelete.append(tf.name)
    else:
        if os.path.splitext(infile)[1] not in ['.fasta','.fa']: #STELLAR is  fussy about file extensions
            tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
            tf.close()
            shutil.copy(infile,tf.name)
            workinginfile = tf.name
            todelete.append(tf.name)
        else:
            workinginfile=infile
    if os.path.splitext(db)[1] not in ['.fasta','.fa']: #STELLAR is fussy about file extensions, need to get db extension correct
        tf = tempfile.NamedTemporaryFile(suffix='.fasta',delete=False)
        tf.close()
        shutil.copy(db,tf.name)
        db = tf.name
        todelete.append(tf.name)
    if hitsOut:
        if os.path.splitext(hitsOut)[1].lower()!='.gff':
            #STELLAR is fussy about file extensions
            raise TypeError('Target for stellar output (hitsOut) must have the extension .gff')
    else:
        tf = tempfile.NamedTemporaryFile(delete=False,suffix='.gff') #STELLAR is fussy about file extensions
        tf.close()
        hitsOut = tf.name
        todelete.append(tf.name)
    print(('STELLAR with %s\n' % infile))
    call_list = [stellarPath]
    call_list=call_list + [db]
    call_list=call_list + [workinginfile]
    call_list=call_list + ['-e',str(err)]
    if mincols:
        call_list=call_list + ['--minLength',str(mincols)]
    call_list=call_list + ['--alphabet',alphabet]
    call_list=call_list + ['-numMatches',str(numMatches)]
    call_list=call_list + ['--verbose']
    if strand.lower() == 'plus' or strand.lower() == 'forward' or strand.lower() == 'f':
        call_list=call_list + ['-f']
    if strand.lower() == 'minus' or strand.lower() == 'reverse' or strand.lower() == 'r':
        call_list=call_list + ['-r']
    call_list=call_list + ['--out',hitsOut]
    #print(call_list)
    call_output = subprocess.run(call_list,text=True)
    if verbose:
        print(call_output.stdout)
        print(call_output.stderr)
    
    hits = set([q.id for q in reptools.hitsParser(hitsOut,filetype='stellar',evalue=evalue)])
    
    with open(infile) as in_handle:
        with open(matched,'w') if matched else reptools.dummy_context_mgr() as match_target,\
         open(matchedfq,'w') if matchedfq else reptools.dummy_context_mgr() as matchfq_target,\
         open(notmatched,'w') if notmatched else reptools.dummy_context_mgr() as unmatch_target,\
         open(notmatchedfq,'w') if notmatchedfq else reptools.dummy_context_mgr() as unmatchfq_target:
            for seqtuple in Parser(in_handle):
                if seqtuple[0].split(' ')[0] in hits:
                    match_target.write('>%s\n%s\n' % seqtuple[0:2])
                    matchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple[0:3])
                else:
                    unmatch_target.write('>%s\n%s\n' % seqtuple[0:2])
                    unmatchfq_target.write('@%s\n%s\n+\n%s\n' % seqtuple[0:3])
    
    reptools.clean_up(todelete)
    if return_counts:
        if matched:
            return(reptools.fastqcounter(matched))
        if matchedfq:
            return(reptools.fastacounter(matchedfq))


def hitsParser(
                hitFile,filetype,evalue=False,mincols=False,newline=None,
                queryidStrip=[],hitidStrip=[],hitidStripstart=[]
                ):
    """
    SWIPE tsv files seem to sometimes have query id followed by CR then TAB and then hit info (TAB-seperated)
    followed by LF
    vsearch u14 files seem to add LF to the end of some hit IDs
    the "*" bits are used only for u14 files, but seem harmless with others
    u14 files are created by usearch with a userdefined output of the format
     -userfields query+target+ql+qlo+qhi+tlo+thi+alnlen+ids+id+qstrand+tstrand+bits+evalue
    b6 files are created by BLAST+ with outfmt=6
    """
    if len(queryidStrip) > 0: raise ValueError('queryidStrip is obsolete')
    if len(hitidStrip) > 0: raise ValueError('hitidStrip is obsolete. It has been replaced by hitidStripstart')
    if filetype.lower() == 'swipe':
        hitidStripstart.append('lcl|')
        row_sep='\t'
        newline = '\n'
        commentstart = None
        from reptools import swipeHsp as Hsp
    elif filetype.lower() in ['b6','blast','blastn']:
        hitidStripstart.append('lcl|')
        row_sep='\t'
        commentstart = None
        from reptools import b6Hsp as Hsp
    elif filetype.lower() in ['usearch','ublast','local','global']:
        row_sep='\t'
        commentstart = None
        from reptools import u14Hsp as Hsp
    elif filetype.lower()=='vsearch':
        row_sep='\t'
        newline = '\n'
        commentstart = None
        from reptools import u14Hsp as Hsp
    elif filetype.lower() == 'stellar':
        row_sep='\t'
        commentstart = None
        from reptools import gffHsp as Hsp
    elif filetype.lower() == 'nhmmer':
        row_sep=None
        commentstart = '#'
        from reptools import nhmmerHsp as Hsp
    else:
        raise ValueError(
                        'Unknown filetype, must specify one of "swipe", "stellar", "b6", "blast", "blast", "nhmmer", '
                        '"usearch", "ublast", "local", "global" (latter two from usearch).\n'
                        )
    currenthsp = None
    workinghit = None
    if os.path.getsize(hitFile)>0:
        with open(hitFile,newline=newline) as in_handle: #newline specified
            for line in in_handle:
                if line.strip()[0]==commentstart or line[0]=='*': continue #skip comment and empty lines
                row = line.strip().split(row_sep)
                query_id = row[0].strip()
                hit_id = row[1].strip()
                for strip_string in hitidStripstart:
                    if hit_id.startswith(strip_string):
                        hit_id = hit_id[len(strip_string):]
                
                if hit_id == '*': continue #skip empty lines (just in case it was missed earlier)
                previoushsp = currenthsp
                currenthsp = Hsp(query_id,hit_id,*row[2:])
                
                if workinghit == None: #first line actions
                    hits_list = []
                    workinghit = Hit(currenthsp.query_id,currenthsp.hit_id,currenthsp.query_len)
                    workinghit.hsps.conditional_append(currenthsp,evalue,mincols)
                elif currenthsp.query_id != previoushsp.query_id: #new query
                    if len(workinghit.hsps)>0: hits_list.append(workinghit)
                    if len(hits_list)>0: yield(SearchResult(previoushsp.query_id,hits_list,previoushsp.query_len))
                    hits_list = []
                    workinghit = Hit(currenthsp.query_id,currenthsp.hit_id,currenthsp.query_len)
                    workinghit.hsps.conditional_append(currenthsp,evalue,mincols)
                elif currenthsp.hit_id != previoushsp.hit_id: #new hit (same query)
                    if len(workinghit.hsps)>0: hits_list.append(workinghit)
                    workinghit = Hit(currenthsp.query_id,currenthsp.hit_id,currenthsp.query_len)
                    workinghit.hsps.conditional_append(currenthsp,evalue,mincols)
                else: #same hit, same query (just a new hsp)
                    workinghit.hsps.conditional_append(currenthsp,evalue,mincols)
            
            if len(workinghit.hsps)>0: hits_list.append(workinghit)
            if len(hits_list)>0: yield(SearchResult(currenthsp.query_id,hits_list,currenthsp.query_len)) 


class gffHsp:
    """
    Takes a gff file produced by Stellar
    Attributes:
    DNAME stellar eps-matches DBEGIN DEND PERCID DSTRAND . ATTRIBUTES
    
     ATTRIBUTES   A list of attributes in the format <tag_name>=<tag>
               separated by ';'

    Attributes are: 

      ID=          Name of the query sequence
      seq2Range=   Begin and end position of the alignment in the query
                   sequence (counting from 1)
      cigar=       Alignment description in cigar format
      mutations=   Positions and bases that differ from the database sequence
    with respect to the query sequence (counting from 1)
    """
    def __init__(
                 self,
                 DNAME,
                 source,
                 feature,
                 DBEGIN,
                 DEND,
                 PERCID,
                 DSTRAND,
                 frame,
                 ATTRIBUTES,
                 dbsize=None,
                 usearch_evalues=False
                 ):
        attributes = ATTRIBUTES.split(';')
        
        self.query_id = attributes[0]
        self.hit_id = DNAME
        
        self.query_start = int([a[10:] for a in attributes if a.startswith('seq2Range=')][0].split(',')[0])
        self.query_end = int([a[10:] for a in attributes if a.startswith('seq2Range=')][0].split(',')[1])
        
        self.evalue = float(next(a for a in attributes if a.startswith('eValue='))[7:])
        
        self.query_len = int(self.query_end-self.query_start)
        
        self.hit_start = int(DBEGIN)
        self.hit_end = int(DEND)
        
        self.ident_fract = float(PERCID)/100
        cigar_info = cigarScorer([a[6:] for a in attributes if a[0:6]=='cigar='][0],float(PERCID))
        self.bitscore = float(cigar_info[0])
        if DSTRAND == '+':
            self.strand = 1
        elif DSTRAND == '-':
            self.strand = -1
            #if strand is -1, make sure that hit start and end are in the correct order
            if self.hit_end>self.hit_start:
                self.hit_start,self.hit_end = self.hit_end,self.hit_start 
        else:
            raise ValueError(
                'DSTRAND value should be + or -.  '
                'It is reported as %s for query %s.\n' % (DSTRAND,self.query_id)
            )
        
        self.aln_span = abs(int(self.hit_end-self.hit_start))


def cigarParser(cigarstring):
    """
    Generator function which returns a tuple:
        first value is count
        second value is length in query 
        third value encodes type:
            match or N = 0
            gap = 1
    """
    import string
    countstrings = []
    gap = None
    count = 0
    querycount = 0
    if cigarstring.strip()[0] not in string.digits:
        raise ValueError('CIGAR string should start with an integer.')
    for char in cigarstring.upper():
        if char in string.whitespace:
            pass
        elif char in string.digits:
            countstrings.append(char)
        elif char in ['M','N','D','I']:
            if char in ['M','N']: #check these
                if gap:
                    yield(count,querycount,1)
                    count = int(''.join(countstrings))
                    querycount = count
                else:
                    count = count + int(''.join(countstrings))
                    querycount = count
                gap = False
            elif char in ['D','I']: #check these
                if not gap and gap is not None:
                    yield(count,querycount,0)
                    count = int(''.join(countstrings))
                    if char == 'I':
                        querycount = count
                else:
                    count = count + int(''.join(countstrings))
                    if char == 'I':
                        querycount = querycount + int(''.join(countstrings))
                gap = True
            countstrings = []
        else:
            raise ValueError('Unexpected character in CIGAR string.\n') 
    if char in string.digits:    
        raise ValueError('CIGAR string should not end with a numeric value.')
    else:
        if gap:
            yield(count,querycount,1)
        else:
            yield(count,querycount,0)

#have overall pct_ident.  This gives error rate with (100-pct_ident)/100 = err_rate.
#err_rate is defined as % error, including both insertions and deletions
#if this is correct, err_rate/100*count gives num_err, and subtracting both insertions and deletions should give mismatch count


def cigarScorer(cigarString,pct_id,match=1,mismatch=-2,opengap=-10,extendgap=-1):
    aligned = 0
    cigarLen = 0
    cigarGap = 0
    gapscore = 0
    for t in cigarParser(cigarString):
        cigarLen += t[0]
        if t[2]==0: 
            aligned+=t[0]
        else:
            gapscore += (opengap + extendgap*(t[0]-1))
            cigarGap += t[0]
    num_err = cigarLen-(pct_id/100.*cigarLen)
    num_mismatched = num_err-cigarGap
    num_matched = aligned-num_mismatched
    score = num_matched*match + num_mismatched*mismatch + gapscore
    if gapscore!=0:
        gapped=True
    else:
        gapped=False
    return(score,gapped)


class hspList(list):
    """
    Subset the list type in order to define an additional method which appends only where evalue and/or mincols
    thresholds are met.
    """
    def conditional_append(self, element,evalue=False,mincols=False):
        if evalue:
            if mincols:
                if element.evalue <= evalue and element.aln_span >= mincols:
                    self.append(element)
            else:
                if element.evalue <= evalue:
                    self.append(element)
        elif mincols:
            if element.aln_span >= mincols:
                self.append(element)
        else:
            self.append(element)


class SearchResult:
    """
    Attributes:
        hits = an empty list (to put hits into)
    """
    def __init__(self,query_id,hits_list,query_len=None):
        self.id = query_id
        self.hits = hits_list
        self.query_len = query_len


class Hit:
    """
    Attributes:
        id = hit ID
        hsps = an empty list (to put hits into)
    """
    def __init__(self,query_id,hit_id,query_len=None):
        self.id = hit_id
        self.query_id = query_id
        self.query_len = query_len
        self.hsps = hspList([])
        

class u14Hsp:
    """
    Attributes:
    """
    def __init__(self,query_id,hit_id,query_len,query_start,query_end,hit_start,hit_end,aln_span,ident_num,ident_fract,query_strand,hit_strand,bitscore,evalue):
        self.query_id = query_id
        self.hit_id = hit_id
        #
        try:
            self.query_len = int(query_len)
        except:
            self.query_len = None
        #
        try:
            self.query_start = int(query_start)
        except:
            self.query_start = None
        #
        try:
            self.query_end = int(query_end)
        except:
            self.query_end = None
        #
        try:
            self.hit_start = int(hit_start)
        except:
            self.hit_start = None
        #
        try:
            self.hit_end = int(hit_end)
        except:
            self.hit_end = None
        #
        try:
            self.aln_span = int(aln_span)
        except:
            self.aln_span = None
        #
        try:
            self.ident_num = int(ident_num)
        except:
            self.ident_num = None
        #
        try:
            self.ident_fract = float(ident_fract)
        except:
            self.ident_fract = None
        #
        if query_strand == '+':
            query_strand_int = 1
        elif query_strand == '-':
            query_strand_int = -1
        else:
            query_strand_int = None
        #
        if hit_strand == '+':
            hit_strand_int = 1
        elif hit_strand == '-':
            hit_strand_int = -1
        else:
            hit_strand_int = None
        #
        self.strand = query_strand_int * hit_strand_int
        #check hit start and end are in correct order:
        if self.strand == 1:
            if self.query_start>self.query_end:
                self.query_start,self.query_end = self.query_end,self.query_start
            if self.hit_start>self.hit_end:
                self.hit_start,self.hit_end = self.hit_end,self.hit_start
        elif self.strand == -1:
            if self.query_end>self.query_start and self.hit_end>self.hit_start:
                self.hit_start,self.hit_end = self.hit_end,self.hit_start
        else:
            pass
        try:
            self.bitscore = float(bitscore)
        except:
            self.bitscore = None
        #
        try:
            if float(evalue)<=0: #because vsearch reports -1 for all nucleotide evalues
                self.evalue = None
            else:
                self.evalue = float(evalue)
        except:
            self.evalue = None


class b6Hsp:
    """
    Attributes:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    """
    
    def __init__(
                self,query_id,hit_id,pct_ident,aln_span,
                mismatch_count,gapopen_count,query_start,
                query_end,hit_start,hit_end,evalue,bitscore
                ):
        self.query_id = query_id
        self.hit_id = hit_id
        #
        self.query_len = None
        #
        try:
            self.query_start = int(query_start)
        except:
            self.query_start = None
        #
        try:
            self.query_end = int(query_end)
        except:
            self.query_end = None
        #
        try:
            self.hit_start = int(hit_start)
        except:
            self.hit_start = None
        #
        try:
            self.hit_end = int(hit_end)
        except:
            self.hit_end =  None
        #
        try:
            self.aln_span = int(aln_span)
        except:
            self.aln_span = None
        #
        try:
            self.ident_num = int(ident_num)
        except:
            self.ident_num = None
        #
        try:
            self.ident_fract = float(pct_ident)/100
        except:
            self.ident_fract = None
        #
        #if query is reported backwards, but hit is forwards, switch all the start and ends (probably not necessary)
        if self.query_start > self.query_end and self.hit_start < self.hit_end: 
            self.query_start,self.query_end = self.query_end,self.query_start
            self.hit_start,self.hit_end = self.hit_end,self.hit_start
        #
        # if hit is on reverse strand, hit_end will now definitely be lower than hit_start: set strand accordingly
        if self.hit_start > self.hit_end:
            self.strand = -1
        else:
            self.strand = 1
        #  
        try:
            self.bitscore = float(bitscore)
        except:
            self.bitscore = None
        #
        try:
            self.evalue = float(evalue)
        except:
            self.evalue = None


class swipeHsp:
    """
    Attributes:
    query_id,hit_id,ident_pct,aln_span,mismatch_num,gapopen_count,query_start,query_end,hit_start,hit_end,evalue,bitscore
    """
    def __init__(self,query_id,hit_id,ident_pct,aln_span,mismatch_num,gapopen_count,query_start,query_end,hit_start,hit_end,evalue,bitscore):
        self.query_id = query_id
        self.hit_id = hit_id
        #
        #try:
        #    self.query_len = int(query_len)
        #except:
        self.query_len = None
        #
        try:
            self.query_start = int(query_start)
        except:
            self.query_start = None
        #
        try:
            self.query_end = int(query_end)
        except:
            self.query_end = None
        #
        try:
            self.hit_start = int(hit_start)
        except:
            self.hit_start = None
        #
        try:
            self.hit_end = int(hit_end)
        except:
            self.hit_end =  None
        #
        try:
            self.aln_span = int(aln_span)
        except:
            self.aln_span = None
        #
        try:
            self.ident_num = self.aln_span - int(mismatch_num)
        except:
            self.ident_num = None
        #
        try:
            self.ident_fract = float(ident_pct)/100
        except:
            self.ident_fract = None
        #
        if self.hit_start > self.hit_end: #if on reverse strand
            self.strand = -1
            #self.hit_start,self.hit_end = self.hit_end,self.hit_start
        elif self.hit_start < self.hit_end:
            self.strand = 1
        else:
            self.strand = None
        #  
        try:
            self.bitscore = float(bitscore)
        except:
            self.bitscore = None
        #
        try:
            self.evalue = float(evalue)
        except:
            self.evalue = None


#not fully implemented
class nhmmerHsp:
    """
    Attributes:
    target_name target_accession query_name query_accession hmmfrom hmmto alifrom alito envfrom envto sqlen strand E-value score bias description_of_target
    """
    def __init__(
                self,hit_id,hit_accession,query_id,
                query_accession,query_start,query_end,
                hit_start,hit_end,env_start,
                env_end,target_length,strand,
                evalue,bitscore,bias,description
                ):
        self.query_id = query_id
        self.hit_id = hit_id
        #
        self.query_len = None
        #
        try:
            self.query_start = int(query_start)
        except:
            self.query_start = None
        #
        try:
            self.query_end = int(query_end)
        except:
            self.query_end = None
        #
        try:
            self.hit_start = int(hit_start)
        except:
            self.hit_start = None
        #
        try:
            self.hit_end = int(hit_end)
        except:
            self.hit_end =  None
        #
        #try:
        #    self.aln_span = int(aln_span)
        #except:
        self.aln_span = None
        #
        #try:
        #    self.ident_num = ident_num
        #except:
        self.ident_num = None
        #
        #try:
        #    self.ident_fract = float(pct_ident)/100
        #except:
        self.ident_fract = None
        #
        #if query is reported backwards, but hit is forwards, switch all the start and ends
        if self.query_start > self.query_end and self.hit_start < self.hit_end: 
            self.query_start,self.query_end = self.query_end,self.query_start
            self.hit_start,self.hit_end = self.hit_end,self.hit_start
        #
        if strand == '-':
            self.strand = -1
        else:
            self.strand = 1
        #  
        try:
            self.bitscore = float(bitscore)
        except:
            self.bitscore = None
        #
        try:
            self.evalue = float(evalue)
        except:
            self.evalue = None




