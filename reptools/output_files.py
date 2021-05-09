import os
import itertools as it
import reptools
import collections


def make_seq2clone(
                   CDR3file,
                   cluster_files,
                   outfile
                  ):
    if reptools.checkFiletype(CDR3file,return_empty=True) == 'fastq':
        from reptools import FASTQparser as Parser
    else:
        from reptools import FASTAparser as Parser
    
    #read the cluster files
    clusters = reptools.read_cluster_files(cluster_files)
    
    #iterate through the CDR3 fastx file
    #for each line, output V, J and CDR3, and also check clusters to see if more than one sequence is tied to this clone
    #  if so, output lines for any extra sequences as well
    with open(CDR3file) as CDR3_handle, open(outfile,'w') as out_handle:
        out_handle.write('seqID\tV\tJ\tCDR3\n')
        for seq_tuple in Parser(CDR3_handle):
            title_split = seq_tuple[0].split(';')
            seqID = title_split[0]
            V = [_ for _ in title_split if _[0:2].upper()=='V='][0].split('=')[1]
            J = [_ for _ in title_split if _[0:2].upper()=='J='][0].split('=')[1]
            out_handle.write('{}\t{}\t{}\t{}\n'.format(seqID,V,J,seq_tuple[1]))
            #see if any more sequences are linked to this one
            if seqID in clusters:
                for title in clusters[seqID]:
                    out_handle.write('{}\t{}\t{}\t{}\n'.format(title,V,J,seq_tuple[1]))


def read_cluster_files(cluster_files):
    """
    cluster_files to be a list of cluster files from reptools denoise, IN THE ORDER IN WHICH THEY WERE GENERATED
    """
    clusters = collections.defaultdict(list)
    
    for fn in cluster_files:
        with open(fn) as cluster_handle:
            for row in cluster_handle:
                row_split=row.strip().split('\t')
                if len(row_split)>1: #anything in this line?
                    recipient = row_split[0]
                    #add donor to recipient entry in clusters, creating it if necessary
                    clusters[recipient].extend(row_split[1:]) #slicing makes a deep copy
                    #if any of the donors have their own dictionary entries from a previous round, add the contents to
                    #  the new recipient sequence, removing the existing entry
                    for donor in row_split[1:]:
                        clusters[recipient].extend(clusters.pop(donor,[]))
    
    return(clusters)
                    
    

def make_VDJtools_dir(
                      indir,
                      outdir=False,
                      VDJout_dir=None,
                      genes=False,
                      emptycols=['D'],
                      overwrite=False,
                      filetype='fastq',
                      ):
    if not genes:
        genes = {'V':'V','J':'J','C':'C'}
    if not outdir:
        outdir = indir
    
    VDJout_dir,_ = reptools.build_path(True, VDJout_dir, 'VDJtools', outdir)

    filetypes = reptools.select_filetypes(filetype)
    typefiles = [fn for fn in os.listdir(indir) if os.path.splitext(fn)[1] in filetypes]
    if len(typefiles)==0:
        print('No files of specified type found.\n')
        return
    
    if overwrite:
        reptools.remove_dir(VDJout_dir,recursive=True)
    reptools.cautious_mkdir(VDJout_dir)
    
    for fn in typefiles:
        outfn = os.path.splitext(fn)[0]+'.tab'
        reptools.make_VDJtools(os.path.join(indir,fn),os.path.join(VDJout_dir,outfn),genes,emptycols,filetype)


def make_VDJtools(CDR3file,outfile,genes=False,emptycols=['D'],filetype='fastq'):
    #in order to reduce memory usage, count the reads in the whole file first (will increase run-time almost 2-fold, though)
    
    genes = reptools.fill_defaults(genes, ['V','J','C'], ['V','J','C'], ordered=True)
    
    if filetype == 'fastq':
        from reptools import FASTQparser as Parser
        from reptools import fastqcounter as Counter
    else:
        from reptools import FASTAparser as Parser
        from reptools import fascounter as Counter
    
    filecount = float(Counter(CDR3file))
    
    with open(CDR3file) as CDR3_handle, open(outfile,'w') as out_handle:
        titlestring = 'count\tfrequency\tCDR3nt\tCDR3aa\t%s\n' % '\t'.join(it.chain(genes.values(),emptycols))
        out_handle.write(titlestring)
        
        for seq_tuple in Parser(CDR3_handle):
            count = reptools.DerepCount(seq_tuple[0])
            frequency = '%.2e' % (count/filecount)
            genecounts = '\t'.join([','.join(reptools.ExtractGene(seq_tuple[0],_gene).split('+')) for _gene in genes.values()])
            emptycols = '\t'.join(['.' for _ in emptycols])
            aa = reptools.CDR3_trans(seq_tuple[1])
            outstring = '%s\t%s\t%s\t%s\t%s\t%s\n' % (count,frequency,seq_tuple[1],aa,genecounts,emptycols)
            out_handle.write(outstring)


def counts_csv(csvfile,paireddirs=(),unpaireddirs=(),pairsuffixes=('_1','_2'),filetype=None,overwrite=False,basename=False):
    if not overwrite and os.path.exists(csvfile):
        raise IOError('Target file already exists. To overwrite, set overwrite=True.')
    
    filetypes = reptools.select_filetypes(filetype)
    counts={}
    for dir in paireddirs:
        filelist = [os.path.join(dir,fn) for fn in os.listdir(dir) if os.path.splitext(fn)[1] in filetypes]
        counts[dir] = {}
        for fn in filelist:
            root = os.path.splitext(os.path.split(fn)[1])[0][:-len(pairsuffixes[0])]
            if root not in list(counts.keys()):
                counts[dir][root]=0
            if os.path.getsize(fn)>0:
                if os.path.splitext(fn)[1] in ['.fastq','.fq']:
                    counts[dir][root] += reptools.fastqcounter(fn)
                else:
                    counts[dir][root] += reptools.fascounter(fn)
    for dir in unpaireddirs:
        counts[dir] = {}
        filelist = [os.path.join(dir,fn) for fn in os.listdir(dir) if os.path.splitext(fn)[1] in filetypes]
        for fn in filelist:
            root = os.path.splitext(os.path.split(fn)[1])[0]
            if os.path.getsize(fn)>0:
                if os.path.splitext(fn)[1] in ['.fastq','.fq']:
                    counts[dir][root] = reptools.fastqcounter(fn)
                else:
                    counts[dir][root] = reptools.fascounter(fn)
            else:
                counts[dir][root] = 0
    #
    allroots = sorted(set([k for dict in counts for k in counts[dict] ]))
    
    if basename:
        titles = [os.path.split(path)[1] for path in paireddirs] + [os.path.split(path)[1] for path in unpaireddirs]
    else:
        titles = paireddirs+unpaireddirs
    
    with open(csvfile,'wb') as out_handle:
        out_handle.write('root,%s\n' % (','.join(titles) ))
        for root in allroots:
            out_handle.write('%s,%s\n' % (root,','.join([str(counts[dir][root]) if root in counts[dir] else '0' for dir in paireddirs+unpaireddirs] ) ) )

