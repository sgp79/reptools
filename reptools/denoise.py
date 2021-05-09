from numba import vectorize
from numba import njit
import numpy as np
import collections
import os
import tempfile
import reptools

#seqlens = False #if False, process all
#weight_by_qual = True
#threshold = 10
#indel_threshold = 100
#infile = 'C:/Users/Stephen Preston/Sync/TCRSeq data/MiSeq48/RD13_Lib4/CDR3/56.fastq'
#FASTQout = False
#FASTAout = False
#change_logs = False
#outdir = False

#seqs_np,gene_ids_sets = denoise(infile,change_logs=None,seqlens=[58,59,60])
#denoise(infile,change_logs=False,seqlens=[58,59,60])
#denoise(infile,change_logs=False)

def denoise_dir(
    indir,
    outdir = False,
    weight_by_qual = True,
    threshold = 10,
    indel_threshold = 100,
    FASTQout_dir = None,
    #FASTAout_dir = False,
    subs = True,
    indels = True,
    deambig = True,
    filetype = 'fastq',
    overwrite = False
    ):
    
    filetypes = reptools.select_filetypes(filetype)
    infiles=[fn for fn in os.listdir(indir) if os.path.splitext(fn.lower())[1] in filetypes]
    
    if len(infiles)==0:
        print('No fastq files found.\n')
        return
    
    if not outdir:
        outdir = indir
    
    FASTQout_dir,_ = reptools.build_path(True, FASTQout_dir, 'denoisedCDR3', outdir)
    FASTAout_dir = False
    #FASTAout_dir = build_path(True, FASTAout_dir, 'denoisedCDR3_fasta', outdir)
    
    #make output directory, deleting pre-existing data is overwrite is set
    for pth in [
                FASTQout_dir,
                FASTAout_dir
                ]:
        if pth: reptools.reptools.cautious_mkdir(pth,overwrite=overwrite)
    
    for fn in infiles:
        FASTQout = reptools.make_unpaired_filepaths(FASTQout_dir, os.path.splitext(fn)[0])
        FASTAout = reptools.make_unpaired_filepaths(FASTAout_dir, os.path.splitext(fn)[0],'fas')
        _ = reptools.denoise_file(
                      os.path.join(indir,fn),
                      weight_by_qual = weight_by_qual,
                      threshold = threshold,
                      indel_threshold = indel_threshold,
                      FASTQout = True,
                      FASTAout = False,
                      FASTQout_fn = FASTQout,
                      change_logs = False,
                      subs = subs,
                      indels = indels,
                      deambig = deambig,
                      overwrite = overwrite
                    )[0]
    
    return(FASTQout_dir)


def denoise_file(
    infile,
    seqlens = False, #if False, process all
    weight_by_qual = True,
    threshold = 10,
    indel_threshold = 100,
    FASTQout = True,
    FASTAout = False,
    FASTQout_fn = None,
    FASTAout_fn = None,
    change_logs = None,
    outdir = False,
    subs = True,
    indels = True,
    deambig = True,
    overwrite = False
):
    """
    """
    inpath = os.path.split(infile)[0]
    inbase = os.path.split(infile)[1]
    filestem = os.path.splitext(inbase)[0]
    
    if not outdir: outdir = inpath
    
    if FASTQout: FASTQout = reptools.build_fn(FASTQout_fn, '{}_denoised.fastq'.format(filestem),outdir)
    if FASTAout: FASTAout = reptools.build_fn(FASTAout_fn, '{}_denoised.fas'.format(filestem),outdir)
    
    if not change_logs:
        if change_logs is None:
            change_logs = [
                            os.path.join(outdir,'{}_dereplicated.clust'.format(filestem)),
                            os.path.join(outdir,'{}_denoised_subs.clust'.format(filestem)),
                            os.path.join(outdir,'{}_denoised_indels.clust'.format(filestem)),
                            os.path.join(outdir,'{}_denoised_genesimplify.clust'.format(filestem))
                          ]
        else:
            change_logs = [False]*4
    elif type(change_logs) not in [list,tuple] or len(change_logs)!=4:
        raise ValueError(
                         'change_logs should either be None (for automatic naming), '
                         'or a list/tuple with four items, one for deplications, one for denoise changes, one for '
                         'indel removal changes, and one for gene simplification changes.  change_logs=False will '
                         'omit log production.'
                         )
    else:
        pass
    
    #load and derep
    seqs,gene_ids =  reptools.derep_FASTQ(infile,change_logs[0])
    
    #produce sets of gene segments
    gene_ids_sets, gene_labels = reptools.process_geneids(gene_ids)
    del(gene_ids)
    
    #TODO = in derep, remove any sequences with ambiguous bases (best to just keep agctAGCT)
    
    #convert to numpy by seqlen
    seqs_np = {}
    for seqlen in list(seqs): #list(seqs) forces a copy of the keys, so that it can cope with on-the-fly deletions
        seqs_np[seqlen] = reptools.build_arrays(seqs[seqlen])
        del(seqs[seqlen]) #to free up memory
    
    del(seqs) 
    
    #each value in seqlen is a list of numpy arrays:
    #   (seq_array, gene_id_array, prob_array, counts_array, seq_names)
    
    reptools.update_counts_dtype(seqs_np) #to make sure the counts arrays all have an appropriate dtype
    
    if subs:
        print('denoising substitutions in {}'.format(inbase))
        reptools.denoise_substitutions(
                                 seqs_np,
                                 gene_ids_sets,
                                 threshold = threshold,
                                 clust_file = change_logs[1],
                                 weight_by_qual = weight_by_qual,
                                 seqlens = seqlens
                               )
    
    if indels:
        print('denoising indels in {}'.format(inbase))
        reptools.denoise_indelonly(
                             seqs_np,
                             gene_ids_sets,
                             threshold = indel_threshold,
                             clust_file = change_logs[2],
                             seqlens = seqlens
                           )
    
    if deambig:
        print('simplifying gene segments in {}'.format(inbase))
        gene_ids_sets = reptools.simplify_genes(
                                        seqs_np,
                                        gene_ids_sets,
                                        clust_file = change_logs[3],
                                        seqlens = seqlens
                                       )
    
    if FASTQout or FASTAout:
        reptools.saveFASTX(seqs_np,gene_ids_sets,FASTAout,FASTQout,gene_labels)
    
    return(FASTQout,FASTAout,change_logs)