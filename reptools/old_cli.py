import argparse
import reptools
import os
import sys
from reptools import build_fn
from reptools import build_path
import platform
from glob import glob

#TODO = sort out help display, which will be odd right now, due to migration from old cli structure
def parse_args():
    parser = argparse.ArgumentParser(
                         description='\n'
                          '\n',
                          usage='reptools [function(s)] -i <input> [...]\n',
                          prog='reptools',
                          formatter_class=argparse.RawDescriptionHelpFormatter
                                )
    
    parser.add_argument(
                          '-V',
                          '--version',
                          action='version',
                          version="{prog}s ({version})".format(prog="%(prog)", version=reptools.__version__)
                        )
    
    #parent parser, to supply shared options
    reptools_functions = ['call','denoise','EEfilter','VDJtools']
    #parent_parser = argparse.ArgumentParser(add_help=False,usage=None,description=None)
    
    requiredFuns = parser.add_argument_group(
                                   'Function(s) - enter one or more:\n'
                                   '          call  Call gene segments and extract CDR3 (takes paired fastq\n'
                                   '                files.  Trims adapters by default.\n'
                                   '       denoise  Denoise and dereplicate CDR3 sequences.  If run without\n'
                                   '                the call function, takes CDR3 fastq file(s).\n'
                                   '      EEfilter  Filter CDR3 sequences by expected error rate.  If run without\n'
                                   '                the call or denoise function, takes CDR3 fastq file(s).\n'
                                   '      VDJtools  Output VDJtools format tab-seperated files.  If run alone,\n'
                                   '                takes denoised/dereplicated and segment-labelled CDR3 fastq\n'
                                   '                file(s). \n'
                                            )
    requiredFuns.add_argument(
                              'functions',
                              nargs='+',
                              choices=reptools_functions
                              )
    
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument('-i', '--input', nargs='+', required=True)
    
    requiredAdaptertrimming = parser.add_argument_group('Arguments required unless --notrim is set')
    requiredAdaptertrimming.add_argument('--adapters', nargs='+', default=False)
    
    
    optionalArgs = parser.add_argument_group('Optional arguments')
    optionalArgs.add_argument(
                              '-d',
                              '--databases',
                              nargs='+',
                              default=False,
                              dest='dbs',
                              help='To call genes, the paths to 3 fasta databases must be supplied: V, J & C, in the '
                              'format V=<path> J=<path> C=<path>\n'
                              'If --databaseDir is set, only the filenames need to specified here.'
                              )
    optionalArgs.add_argument(
                              '-g',
                              '--genedict',
                              default=False,
                              help='To call genes, the path to a csv gene dictionary must be supplied'
                              )
    optionalArgs.add_argument(
                               '-o',
                               '--outdir',
                               nargs='?',
                               default=False,
                               help='if not supplied, defaults to a subdir of the current path called "reptools_output"'
                              )
    optionalArgs.add_argument(
                              '--notrim',
                              action='store_true',
                              help='Skip adapter trimming before calling gene segments'
                              )
    optionalArgs.add_argument('--noCDR3', action='store_true',help='Do not extract CDR3 sequences after calling genes.')
    
    denoise_types = [
                    'all',
                    'substitution',
                    'indel',
                    'genes'
    ]
    optionalArgs.add_argument(
                              '--denoiseModes',
                              nargs='+',
                              default='all',
                              choices=denoise_types,
                              help='Select which forms of denoising to perform: one or more or subsititution, indels, '
                              'and genes (the latter combines identical sequences which have gene segments in common '
                              'for all genes.\n  Defaults to "all".'
                            )
    optionalArgs.add_argument(
                              '--databaseDir',
                              default=False,
                              dest='db_dir',
                              help='If supplied, look here for the database files'
                              )
    optionalArgs.add_argument(
                               '-t',
                               '--threads',
                               default=False,
                               type=int,
                               help='Number of threads to use for alignment (defaults to number of available CPUs)'
                               )
    optionalArgs.add_argument(
                              '--titleSplit',
                              default=' ',
                              dest='title_split',
                              help='Substring seperating sequence ID from metadata in fastq title lines.  Defaults to '
                                   'a space.'
                              )
    optionalArgs.add_argument('--pairSuffixes', nargs=2, default=['_1','_2'], dest='pairsuffixes')
    optionalArgs.add_argument(
                               '--geneLabels',
                               nargs='+',
                               default=False,
                               dest = 'labels',
                               help='How to label to genes in the title line, supplied in form "--geneLabels V=XXX '
                               'J=XXX C=XXX", where XXX is the desired label.  Any or none may be supplied. Defaults '
                               'to V=V J=J C=C'
                               )
    optionalArgs.add_argument('--overwrite', action='store_true')
    optionalArgs.add_argument(
                              '--tiebreaker',
                              default = 'evalue',
                              choices=['evalue','score'], 
                              help='Value to pick best aligner hits according to: (lowest) evalue or (highest) score. '
                              'Default is evalue.'
                              )
    optionalArgs.add_argument(
                                '--aligners',
                                nargs='+',
                                default=False,
                                help='Aligner to use for each step.\nDefault is C=blastn J=swipe V=blastn C104=swipe'
                              )
    optionalArgs.add_argument(
                              '--alignerPaths', nargs='+', default=False, dest='aligner_paths',
                              help='Paths to aligners and associated files.\n'
                                   'Default is blastn=blastn swipe=swipe makeblastdb=makeblastdb bbduk=bbduk.sh'
                              )
    optionalArgs.add_argument(
                               '--denoiseThreshold', type=float, default=10, dest='threshold',
                               help='In the subsitution mode of denoise, reads are merged with those which are '
                                    'at least (hamming distance)*denoiseThreshold more common.  Default is 10.'
                              )
    optionalArgs.add_argument(
                               '--denoiseThresholdIndels', type=float, default=100, dest='indel_threshold',
                               help='In the indel mode of denoise, reads are merged with those which are '
                                    'one indel away and at least denoiseThresholdIndels more common.  Default is 100.'
                              )
    optionalArgs.add_argument(
                               '--CDR3maxee', type=float, default=1, dest='CDR3maxee',
                               help='The EEfilter will filter CDR3 to remove those with more than --CDR3maxee expected '
                                    'errors, according to their qual string.'
                            )
    optionalArgs.add_argument(
                               '--C103db_length', type=int, default=30, dest='C103db_length',
                               help="For CDR3 identification, the 5' ends of the supplied V genes are trimmed, to leave"
                               " only this number of nt (the 3' end).  Defaults to 30, but higher values may be "
                               "required, especially for genes with gene conversion and/or hypermutation."
                            )
    optionalArgs.add_argument(
                               '--titleFormat', default = 'clusterID|readID', dest='title_format',
                               help='Paired file integrity checks depend on the fastq title line format. The default '
                                    'is "clusterID|readID", which adequately represents illumina output.  For '
                                    'mixcr-tools synthetic data, use "readID|clusterID".  Use | to represent the title '
                                    'seperator, the actual identity of which is set with --titleSplit (which defaults '
                                    'to a space, but should be set to "|" for mixcr tools).'
                            )
    optionalArgs.add_argument(
                               '--blastdb_version', type=int, default=5, dest='blastdb_version', choices={4, 5},
                               help='Version 5 BLAST databases are used for BLAST calls by default, but a bug in '
                               'BLAST+ can lead to this producing an out of memory error on some systems.  If this '
                               'occurs, set "--blastdb_version 4".',
                              )
    
    alignerArgs = parser.add_argument_group(
                'Aligner arguments.\n'
                '  If supplied, these should each give one or more entries in the form:\n'
                '    XXX=[value], \n'
                '  where XXX is one of: [C, J, V, C104].\n'
                '  The --mincols argument also accepts a value for V_Fread' 
                )
    alignerArgs.add_argument('--mincols', nargs='+', default=False)
    alignerArgs.add_argument('--id', nargs='+', default=False)
    alignerArgs.add_argument('--evalue', nargs='+', default=False)
    alignerArgs.add_argument('--wordlength', nargs='+', default=False)
        
    # arguments to specify output subdirectories (these all have sensible defaults), or to suppress outputs
    dirArgs = parser.add_argument_group(
                                        'Output subdirectories (THESE DO NOT HAVE TO BE SUPPLIED):\n'
                                        'By default, reptools will create subdirectories within the output dir.\n'
                                        'Set values to override this behaviour - supply a relative path to output to\n'
                                        'a subdir of that name within the output directory, or an absolute path to\n'
                                        'output elsewhere.\n'
                                        'Set to False to omit this output entirely (it will be stored in a temporary\n'
                                        'directory, and deleted on exit).\n'
                                        )
    dirArgs.add_argument(
                          '--adaptertrimmedDir', nargs = '?', default=None, const=None, dest='adaptertrimmed_dir',
                          help='Adapter trimmed files will be saved to adaptertrimmedDir.\n'
                          'Defaults to outdir/adaptertrimmed\n'
                          'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument(
                        '--fastqDir', nargs = '?', default=None, const=None, dest='hits_dir',
                        help='Gene segment-annotated fastq files will be saved to fastqDir.\n'
                        'Defaults to outdir/fastq_hits\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument('--rawCDR3Dir', nargs = '?', default=None, const=None, dest='CDR3_dir',
                         help='Extracted CDR3 fastq files will be saved to rawCDR3Dir.\n'
                        'Defaults to outdir/rawCDR3\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument('--denoisedCDR3Dir', nargs = '?', default=None, const=None, dest='denoised_dir',
                         help='Denoised CDR3 fastq files will be saved to denoisedCDR3Dir.\n'
                        'Defaults to outdir/denoisedCDR3\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument(
                        '--filtCDR3Dir', nargs = '?', default=None, const=None, dest='filtCDR3_dir',
                        help='Expected error-filtered CDR3 fastq files will be saved to filtCDR3Dir.\n'
                        'Defaults to outdir/EEfilteredCDR3\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument(
                         '--VDJtoolsDir', nargs = '?', default=None, const=None, dest='VDJtools_dir',
                         help='VDJtools files will be saved to VDJtoolsDir.\n'
                        'Defaults to outdir/VDJtools\n'
                        'Set to False if output should not be saved.'
                         )
    
    args = parser.parse_args()
    
    #function selection logic
    if 'call' in args.functions and len(args.functions)>1 and args.noCDR3:
        raise ValueError('--noCDR3 may only be set when the call function is used on its own.\n'
                         'If subsequent functions are required, but raw CDR3 output files are not desired, set '
                         'rawCDR3Dir (or rawCDR3File)\n'
                         'to None.')
    
    if 'call' in args.functions and not args.dbs: 
        raise ValueError('If "call" is set as a function, --databases/-d is required.')
    
    if 'call' in args.functions and not args.genedict:
        raise ValueError('If "call" is set as a function, --genedict/-g is required.')
    
    #parse tiebreaker
    if args.tiebreaker == 'evalue': args.tiebreaker = ['evalue','lower']
    if args.tiebreaker == 'score': args.tiebreaker = ['score','higher']
        
    #denoiser function selection logic
    if args.denoiseModes == 'all':
        args.denoise_substitution, args.denoise_indel, args.denoise_gene_segments = True, True, True
    else:
        args.denoise_substitution, args.denoise_indel, args.denoise_gene_segments = False, False, False
        if 'substitution' in args.denoiseModes: args.denoise_substitution=True
        if 'indel' in args.denoiseModes: args.denoise_indel=True
        if 'genes' in args.denoiseModes: args.denoise_gene_segments=True
    
    #check that adapters were supplied (unless trimAdapters is to be omitted)
    if 'call' in args.functions and not args.notrim and not args.adapters:
        raise ValueError('--adapters must give the sequences of adapters to trim, unless --notrim is set.')
    
    #parse databases
    if 'call' in args.functions:
        #check length has either two or three values (V/J or V/J/C)
        #if len(args.dbs) != 3:
        #    raise ValueError('Three gene databases must be supplied: V, J & C, in the format V=<path> J=<path> C=<path>')
        args.dbs = {i.split('=')[0].upper():i.split('=')[1] for i in args.dbs}
        if 'V' not in args.dbs and 'C' not in args.dbs:
            raise ValueError(
                             "At least one of 'V' and 'C' databases must be supplied, "
                             "in the format --db V=<path> C=<path>"
                             )
        if ('V' not in args.dbs or 'J' not in args.dbs) and rawCDR3Dir:
            raise ValueError(
                             "If rawCDR3Dir is not set to None, both 'V' and 'J' databases must be supplied, "
                             "in the format --db V=<path> J=<path>"
                             )
        if ('V' not in args.dbs or 'J' not in args.dbs) and (
                                                             'denoise' in args.functions or 
                                                             'EEfilter' in args.functions or
                                                             'VDJtools' in args.functions
                                                             ):
            raise ValueError(
                             "For denoise, EEfilter, or VDJtools functions to be performed, CDR3 need to be generated. "
                             "CDR3 generation requires both 'V' and 'J' databases to be supplied, "
                             "in the format --db V=<path> J=<path>"
                             )
    #parse aligners
    args.aligners = parse_gene_setting(args.aligners,'aligners',['C','J','V','C104'],str)
    valid_aligners = ['blastn','swipe']
    if args.aligners:
        for k in args.aligners:
            if args.aligners[k] not in valid_aligners:
                raise ValueError(
                                 "Unknown aligner type supplied in aligners: {}\n"
                                 "Valid aligners are: '{}'".format(k,"','".join(valid_aligners))
                                 )
    
    #parse alignerPaths
    if args.aligner_paths: #if it's False (the default), leave as False
        valid_aligners = ['blastn','swipe','vsearch','makeblastdb','bbduk']
        alignerPaths = [i.split('=') for i in args.aligner_paths]
        args.aligner_paths = {i[0].lower():i[1] for i in alignerPaths}
        for k in args.aligner_paths:
            if k not in valid_aligners:
                raise ValueError(
                                 "Unknown aligner type supplied in alignerPaths: {}\n"
                                 "Valid aligners are: '{}'".format(k,"','".join(valid_aligners))
                                 )
    
    #parse mincols
    args.mincols = parse_gene_setting(args.mincols,'mincols',['C','J','V','V_Fread','C104'],int)
    
    #parse id
    args.id = parse_gene_setting(args.id,'id',['C','J','V','C104'],float)
    
    #parse evalue
    args.evalue = parse_gene_setting(args.evalue,'evalue',['C','J','V','C104'],float)
    
    #parse wordlength
    args.wordlength = parse_gene_setting(args.wordlength,'wordlength',['C','J','V','C104'],int)
    
    #parse labels
    args.labels = parse_gene_setting(args.labels,'labels',['C','J','V'],str)
    
    ################
    #rewritten to handle lists of input files
    #glob to unpack wildcard-conatining input
    #if platform.system() == 'Windows': #linux seems to need specific glob unpacking too
    globs = [glob(i) for i in args.input]
    args.input = [os.path.abspath(_) for g in globs for _ in g] #flatten list of glob outputs, and make abspaths
    if len(args.input) == 0:
        raise ValueError('No files found.  Have you specified the path/wildcards correctly?')
    else:
        print('Will process the following files: {}\n'.format(','.join([os.path.basename(_) for _ in args.input])))
    
    #set outdir defaults, if they were not set in arguments
    if not args.outdir:
        args.outdir = os.path.abspath('reptools_output')
    elif args.outdir.lower() == 'false':
        raise ValueError('False is not a valid entry for outdir.  If it was your intention to use an\n'
                         'output directory called "False", why not try something less confusing?\n'
                         'If it was your intention to use the default value, please omit the --outdir\n'
                         'entry entirely.')
    else:
        args.outdir = os.path.abspath(args.outdir)
    
    
    #set output subdirectory paths
    tempdirs = []
    
    if 'call' in args.functions:
        z = build_path(True, selected_dir=args.adaptertrimmed_dir, default_dir='adaptertrimmed', outdir=args.outdir)
        #This builds an abspath to a directory to store adapter-trimmed sequences in.
        #A tuple is returned, of which the first element is a path, and the second element is either None (if the 
        #  path is NOT a temporary dir), or the same as the first element (if the path IS a temporary dir).
        #   - if args.adaptertrimmed_dir is None (the default), the path will be a subdir of args.outdir called 
        #       "adapterttrimmed"
        #   - if args.adaptertrimmed_dir is False (the default), the path will be a temporary dir
        #   - if args.adaptertrimmed_dir is a absolute path, it will be used as-is
        #   - if args.adaptertrimmed_dir is a relative path, it will be appended to the main output dir
        args.adaptertrimmed_dir = z[0] #update adaptertrimmed_dir with path
        tempdirs.append(z[1]) #if the path is a temporary dir, add it to the list of dirs to remove on cleanup
        
        z = build_path(True, selected_dir=args.hits_dir, default_dir='fastq_hits', outdir=args.outdir) #as above
        args.hits_dir = z[0] #update hits_dir with path
        tempdirs.append(z[1]) #if the path is a temporary dir, add it to the list of dirs to remove on cleanup
        
        z = build_path(True, selected_dir=args.CDR3_dir, default_dir='rawCDR3', outdir=args.outdir) #as above
        args.CDR3_dir = z[0] #update hits_dir with path
        tempdirs.append(z[1]) #if the path is a temporary dir, add it to the list of dirs to remove on cleanup
    
    if 'denoise' in args.functions:
        z = build_path(True, selected_dir=args.denoised_dir, default_dir='denoisedCDR3', outdir=args.outdir) #as above
        args.denoised_dir = z[0] #update denoised_dir with path
        tempdirs.append(z[1]) #if the path is a temporary dir, add it to the list of dirs to remove on cleanup
    
    if 'EEfilter' in args.functions:
        z = build_path(True, selected_dir=args.filtCDR3_dir, default_dir='EEfilteredCDR3', outdir=args.outdir) #as above
        args.filtCDR3_dir = z[0] #update filtCDR3_dir with path
        tempdirs.append(z[1]) #if the path is a temporary dir, add it to the list of dirs to remove on cleanup
    
    if 'VDJtools' in args.functions:
        z = build_path(True, selected_dir=args.VDJtools_dir, default_dir='VDJtools', outdir = args.outdir) #as above
        args.VDJtools_dir = z[0] #update denoised_dir with path
        tempdirs.append(z[1]) #if the path is a temporary dir, add it to the list of dirs to remove on cleanup
        #not sure it CAN be a temporary path, but leave line above in for now
    
    
    #process title line format, to idenify position of cluster ID
    title_format_split = args.title_format.lower().split('|')
    try:
        args.clusterID_position = title_format_split.index('clusterid')
    except ValueError:
        raise ValueError('If --titleFormat is set, it must include clusterID')
    
    return(args,tempdirs)


def parse_gene_setting(arg,argname,valid_keys,fun):
    if arg: #if it's False (the default), leave as False
        split_list = [i.split('=') for i in arg]
        arg = {i[0]:fun(i[1]) for i in split_list}
        for k in arg:
            if k not in valid_keys:
                raise ValueError(
                                 "Unknown key supplied in {}: {}\n"
                                 "Valid keys are: '{}'".format(argname, k,"','".join(valid_keys))
                                 )
    return(arg)


def main():
    import shutil
    import tempfile
    
    args,tempdirs = parse_args()
    
    if 'call' in args.functions:
        file_pairs = reptools.assign_filepairs(args.input,args.pairsuffixes) #returns a dictionary of filepairs, keyed by file stem 
        
        call_output = reptools.call_pairslist(
                                        filepairs = file_pairs,
                                        dbs = args.dbs,
                                        genedictfile = args.genedict,
                                        db_dir = args.db_dir,
                                        notrim = args.notrim,
                                        adaptertrimmed_dir = args.adaptertrimmed_dir,
                                        hits_dir = args.hits_dir,
                                        CDR3_dir = args.CDR3_dir,
                                        Cmiss_dir = False,
                                        ChitVmiss_dir = False,
                                        ChitVhitJmiss_dir = False,
                                        Vsegmentout_dir = False,
                                        filetype = 'fastq',
                                        pairsuffixes = args.pairsuffixes,
                                        title_split = args.title_split,
                                        labels = args.labels,
                                        mincols = args.mincols,
                                        id = args.id,
                                        evalue = args.evalue,
                                        wordlength = args.wordlength,
                                        gapopen = False,
                                        gapextend = False,
                                        aligners = args.aligners,
                                        aligner_paths = args.aligner_paths,
                                        adapters = args.adapters,
                                        threads = args.threads,
                                        Vdb_length = args.C103db_length,
                                        tiebreaker = args.tiebreaker,
                                        overwrite = args.overwrite,
                                        clusterID_position = args.clusterID_position,
                                        blastdb_version = args.blastdb_version
                                        )
        args.input = call_output #pipe output paths to next function, by overwriting input paths
        
    if 'denoise' in args.functions:
        print('DENOISING')
        denoise_output = reptools.denoise_filelist(
                                                args.input,
                                                FASTQout_dir = args.denoised_dir,
                                                subs = args.denoise_substitution,
                                                indels = args.denoise_indel,
                                                deambig = args.denoise_gene_segments,
                                                weight_by_qual = True,
                                                threshold = args.threshold,
                                                indel_threshold = args.indel_threshold,
                                                overwrite = args.overwrite
                                              )
        args.input = denoise_output #pipe output paths to next function, by overwriting input paths
    
    if 'EEfilter' in args.functions:
        print('FILTERING BY EXPECTED ERROR RATE')
        EE_output = reptools.EEfilter_filelist(
                                           args.input,
                                           outdir = args.outdir,
                                           FASTQout_dir = args.filtCDR3_dir,
                                           maxee = args.CDR3maxee,
                                           overwrite = args.overwrite
                                          )
        args.input = EE_output #pipe output paths to next function, by overwriting input paths
    
    if 'VDJtools' in args.functions:
        print('WRITING VDJ FORMAT FILES')
        reptools.make_VDJtools_filelist(
                            args.input,
                            outdir = args.outdir,
                            VDJout_dir = args.VDJtools_dir,
                            genes = args.labels,
                            emptycols = ['D'],
                            overwrite = args.overwrite
                            )
    
    #remove temporary output files - TODO: put this in an "on exit" function, to handle cleanup on crashes etc.
    reptools.clean_up_dirs(tempdirs)


if __name__ == '__main__':
    sys.exit(main())