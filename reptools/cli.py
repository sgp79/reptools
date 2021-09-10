import argparse
import reptools
import os
from reptools import build_fn
from reptools import build_path
import collections

def parse_args():
    parser = argparse.ArgumentParser(
                         description='\n'
                          'mode:\n'
                          '          file  Run on the input file or filepair\n'
                          '           dir  Run on all fastq files found in the input directory\n'
                          '                    set flag -r/--recursive to search subdirs.\n'
                          '\n'
                          'function(s) - enter one or more:\n'
                          '          call  Call gene segments and extract CDR3 (takes pairs(s) of fastq\n'
                          '                files.  Trims adapters by default.\n'
                          '       denoise  Denoise and dereplicate CDR3 sequences.  If run without\n'
                          '                the call function, takes CDR3 fastq file(s).\n'
                          '      EEfilter  Filter CDR3 sequences by expected error rate.  If run without\n'
                          '                the denoise function function, takes CDR3 fastq file(s).\n'
                          '      VDJtools  Output VDJtools format tab-seperated files.  If run alone, takes \n'
                          '                denoised/dereplicated and segment-labelled CDR3 fastq file(s). \n'
                          '                If run with other functions, output will be VDJfiles produced \n'
                          '                by the previous function to be run.\n'
                          '\n'
                          'options:\n'
                          '    -i/--input  Input file/filepair or directory(ies) (required)\n'
                          '                    - for the call function in file mode, two files must be\n'
                          '                      supplied (read_1 and its read_2 pair).\n'
                          '         [...]  Further options:\n'
                          '                    use "reptools file --help" or "reptools dir --help" to see these\n'
                          '\n',
                          usage='reptools [mode] [function(s)] -i <input> [...]\n',
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
    parent_parser = argparse.ArgumentParser(add_help=False,usage=None,description=None)
    requiredFuns = parent_parser.add_argument_group('Functions') #want to hide all help for these
    requiredFuns.add_argument('functions',nargs='+',choices=reptools_functions)
    
    requiredNamed = parent_parser.add_argument_group('Required arguments')
    requiredNamed.add_argument('-i', '--input', nargs='+', required=True)
    
    requiredAdaptertrimming = parent_parser.add_argument_group('Arguments required unless --notrim is set')
    requiredAdaptertrimming.add_argument('--adapters', nargs='+', default=False)
    
    
    denoise_types = [
                    'all',
                    'substitution',
                    'indel',
                    'genes'
    ]
    
    optionalArgs = parent_parser.add_argument_group('Optional arguments')
    optionalArgs.add_argument(
                              '-d',
                              '--databases',
                              nargs='+',
                              default=False,
                              dest='dbs',
                              help='To call genes, the paths to 3 fasta databases must be supplied: V, J & C, in the format V=<path> J=<path> C=<path>'
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
                               nargs='+',
                               default=False,
                               help='if not supplied, defaults to --input (or its parent directory, in file mode)'
                              )
    optionalArgs.add_argument(
                              '--notrim',
                              action='store_true',
                              help='Skip adapter trimming before calling gene segments'
                              )
    optionalArgs.add_argument('--noCDR3', action='store_true',help='Do not extract CDR3 sequences after calling genes.')
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
    
    alignerArgs = parent_parser.add_argument_group(
                'Aligner arguments: if supplied, these should each give one or more entries in the form XXX=[value], \n'
                'where XXX is one of: [C, J, V, C104].  --mincols additionally accepts V_Fread' 
                )
    alignerArgs.add_argument('--mincols', nargs='+', default=False)
    alignerArgs.add_argument('--id', nargs='+', default=False)
    alignerArgs.add_argument('--evalue', nargs='+', default=False)
    alignerArgs.add_argument('--wordlength', nargs='+', default=False)
    
    #subparsers
    #subparsers = parser.add_subparsers(help='test',dest='mode')
    
    subparsers = parser.add_subparsers(help=argparse.SUPPRESS,dest='mode')
    #subparsers.required = True
    
    # create the parser for the "file" command
    parser_file = subparsers.add_parser('file', parents = [parent_parser], help='file help')
    fileArgs = parser_file.add_argument_group('file mode-specific arguments.\nDefaults will create paths, set to '
                                              'False to avoid producing files.')
    fileArgs.add_argument(
                          '--adaptertrimmedFilestem', default=None, dest='adaptertrimmed_filestem',
                          help='Adapter trimmed files will be saved to outdir as '
                          '[adaptertrimmedFilestem][pairsuffix].fastq'
                          'Defaults to [input filestem (without pairsuffix)]_trimmedfastq[pairsuffix].fastq\n'
                          'Set to False if output should not be saved.'
                          )
    fileArgs.add_argument(
                          '--fastqFilestem', default=None, dest='hits_filestem',
                          help='Gene segment-annotated fastq files will be saved to outdir as '
                          '[fastqFilestem]_F.fastq and [fastqFilestem]_R.fastq'
                          'Defaults to [input filestem (without pairsuffix)]_F/R.fastq\n'
                          'Set to False if output should not be saved.'
                          )
    fileArgs.add_argument(
                          '--rawCDR3File', default=None, dest='CDR3_file',
                          help='Extracted CDR3 will be saved to outdir as rawCDR3File.\n'
                          'Defaults to [input filestem (without pairsuffix)]_CDR3.fastq\n'
                          'Set to False if output should not be saved.'
    )
    fileArgs.add_argument('--denoisedCDR3File', default=None, dest='denoised_file',
                          help='Denoised CDR3 will be saved to outdir as denoisedCDR3File.\n'
                               'Defaults to [input filestem (without pairsuffix)]_denoised.fastq\n'
                               'Set to False if output should not be saved.'
                          )
    fileArgs.add_argument('--filtCDR3File', default=None, dest='filtCDR3_file',
                          help='EEfiltered CDR3 will be saved to outdir as filtCDR3File.\n'
                               'Defaults to [input filestem (without pairsuffix)]_EEfilt.fastq\n'
                              'Set to False if output should not be saved.'
                          )
    fileArgs.add_argument('--VDJtoolsFile', default=None, dest='VDJtools_file',
                          help='VDJtools output will be saved to outdir as VDJtoolsFile.\n'
                               'Defaults to [input filestem (without pairsuffix)]_VDJtools.tab'
                          )
                          
    fileArgs.add_argument('--seq2cloneFile', default=False, dest='seq2clone_file',
                          help='VDJOutput tying each assigned input sequence to its assigned V gene, J gene and CDR3 '
                               'will be saved to outdir as seq2cloneFile.\n'
                               'Default = not produced'
                          )

    # create the parser for the "dir" command
    parser_dir = subparsers.add_parser('dir', parents = [parent_parser], help='dir help')
    dirArgs = parser_dir.add_argument_group('dir mode-specific path arguments.\nDefaults will create paths, set to '
                                            ' False to avoid producing files.')
    dirArgs.add_argument(
                          '--adaptertrimmedDir', nargs = '+', default=None, dest='adaptertrimmed_dir',
                          help='Adapter trimmed files will be saved to adaptertrimmedDir.\n'
                          'Defaults to outdir/adaptertrimmed\n'
                          'Supply an absolute path, or a relative path to append to outdir\n'
                          'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument(
                        '--fastqDir', nargs = '+', default=None, dest='hits_dir',
                        help='Gene segment-annotated fastq files will be saved to fastqDir.\n'
                        'Defaults to outdir/fastq_hits\n'
                        'Supply an absolute path, or a relative path to append to outdir\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument('--rawCDR3Dir', nargs = '+', default=None, dest='CDR3_dir',
                         help='Extracted CDR3 fastq files will be saved to rawCDR3Dir.\n'
                        'Defaults to outdir/rawCDR3\n'
                        'Supply an absolute path, or a relative path to append to outdir\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument('--denoisedCDR3Dir', nargs = '+', default=None, dest='denoised_dir',
                         help='Denoised CDR3 fastq files will be saved to denoisedCDR3Dir.\n'
                        'Defaults to outdir/denoisedCDR3\n'
                        'Supply an absolute path, or a relative path to append to outdir\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument(
                        '--filtCDR3Dir', nargs = '+', default=None, dest='filtCDR3_dir',
                        help='Expected error-filtered CDR3 fastq files will be saved to filtCDR3Dir.\n'
                        'Defaults to outdir/EEfilteredCDR3.\n'
                        'Supply an absolute path, or a relative path to append to outdir\n'
                        'Set to False if output should not be saved.'
                        )
    dirArgs.add_argument(
                         '--VDJtoolsDir', nargs = '+', default=None, dest='VDJtools_dir',
                         help='VDJtools files will be saved to VDJtoolsDir.\n'
                        'Defaults to outdir/VDJtools.\n'
                        'Supply an absolute path, or a relative path to append to outdir\n'
                         )
    
    args = parser.parse_args()
    
    #catch missing subparser (mode argument) - this rather than have it set to required in order to provide more useful
    #help 
    if not args.mode:
        print()
        print('Error: Mode (dir or file) must be supplied.\n')
        parser.print_help()
        # parser.print_usage() # for just the usage line
        parser.exit()
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
        alignerPaths = [i.split('=') for i in args.alignerPaths]
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
    
    #build directory lists for dir mode:
    def set_list_lens(outlist,required_len,name):
        if outlist is None:
            return([None]*required_len)
        elif not outlist:
            return([False]*required_len)
        else:
            if len(outlist)==1 and outlist[0].lower()=='false':
                return([False]*required_len)
            elif len(outlist)==1 and outlist[0].lower()=='none':
                return([None]*required_len)
            elif len(outlist)!=required_len:
                raise ValueError(
                                 'Error with {}: If more than one input directory is supplied, a matching number of '
                                 'each output directory is also required (or output paths should be omitted, for '
                                 'default output, or set to None, where output is not required).'.format(name)
                                 )
            else:
                outlist = [i if i.lower()!='false' else False for i in outlist]
                return(outlist)
    
    args.tempdirs = []
    if args.mode=='dir':
        required_len = len(args.input)
        args.outdir = set_list_lens(args.outdir,required_len,'args.outdir')
        for n,outdir in enumerate(args.outdir): 
            #set each outdir to the matching indir, if it isn't already set
            if not outdir:
                args.outdir[n] = args.input[n]
            #and make each outdir an absolute path
            args.outdir[n]  = os.path.abspath(args.outdir[n])
        if 'call' in args.functions:
            args.adaptertrimmed_dir = set_list_lens(args.adaptertrimmed_dir, required_len,'args.adaptertrimmed_dir')
            z = [
                 build_path(
                            True, subdir, 'adaptertrimmed', outdir
                            ) for subdir,outdir in zip(args.adaptertrimmed_dir,args.outdir)
                ]
            args.adaptertrimmed_dir = [t[0] for t in z]
            args.tempdirs.extend([t[1] for t in z])
            
            args.hits_dir = set_list_lens(args.hits_dir, required_len, 'args.hits_dir')
            z = [
                 build_path(
                            True, subdir, 'fastq_hits', outdir
                            ) for subdir,outdir in zip(args.hits_dir,args.outdir)
                ]
            args.hits_dir = [t[0] for t in z]
            args.tempdirs.extend([t[1] for t in z])
            
            args.CDR3_dir = set_list_lens(args.CDR3_dir, required_len, 'args.CDR3_dir')
            z = [
                 build_path(
                            True, subdir, 'rawCDR3', outdir
                            ) for subdir,outdir in zip(args.CDR3_dir,args.outdir)
                ]
            args.CDR3_dir = [t[0] for t in z]
            args.tempdirs.extend([t[1] for t in z])
        
        if 'denoise' in args.functions:
            args.denoised_dir = set_list_lens(args.denoised_dir, required_len, 'args.denoised_dir')
            z = [
                 build_path(
                            True, subdir, 'denoisedCDR3', outdir
                            ) for subdir,outdir in zip(args.denoised_dir,args.outdir)
                ]
            args.denoised_dir = [t[0] for t in z]
            args.tempdirs.extend([t[1] for t in z])
        
        if 'EEfilter' in args.functions:
            args.filtCDR3_dir = set_list_lens(args.filtCDR3_dir, required_len, 'args.filtCDR3_dir')
            z = [
                 build_path(
                            True, subdir, 'EEfilteredCDR3', outdir
                            ) for subdir,outdir in zip(args.filtCDR3_dir,args.outdir)
                ]
            args.filtCDR3_dir = [t[0] for t in z]
            args.tempdirs.extend([t[1] for t in z])
        
        if 'VDJtools' in args.functions:
            args.VDJtools_dir = set_list_lens(args.VDJtools_dir, required_len, 'args.VDJtools_dir')
            z = [
                 build_path(
                            True, subdir, 'VDJtools', outdir
                            ) for subdir,outdir in zip(args.VDJtools_dir,args.outdir)
                ]
            args.VDJtools_dir = [t[0] for t in z]
            args.tempdirs.extend([t[1] for t in z])
    
    #check list lengths for file mode:
    if args.mode=='file':
        if args.outdir and len(args.outdir)>1:
            raise ValueError('Only one outdir may be specified in file mode.')
        elif args.outdir:
            args.outdir = args.outdir[0] #want a string, not a list for file mode
        elif args.outdir==False or args.outdir==None:
            pass 
        else:
            args.outdir=False
        if 'call' in args.functions and len(args.input)!=2:
            raise ValueError(
                             'If the call function is invoked in file mode, -i must be supplied with two file names '
                             '(read_1 and read_2).'
                             )
        if 'call' not in args.functions and len(args.input)!=1:
            raise ValueError(
                             'If the denoise, EEfilter, or VDJtools functions are invoked in file mode in the absence of '
                             'the call function, -i must be supplied with just a single file name (a CDR3 fastq file).'
                             )
        
        #for file mode, check that denoise has been called if --seq2cloneFile is set
        if 'denoise' not in args.functions and args.seq2clone_file:
            raise ValueError(
                         'For seq2cloneFile output, the denoise function must be envoked (other functions may be '
                         'invoked at the same time.'
                         )
    
    #process title line format, to idenify position of cluster ID
    title_format_split = args.title_format.lower().split('|')
    try:
        args.clusterID_position = title_format_split.index('clusterid')
    except ValueError:
        raise ValueError('If --titleFormat is set, it must include clusterID')
    return(args)


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
    
    args = parse_args()
    
    if args.mode=='file' and args.outdir:
        #to prevent trouble later with abspath detection in default filename generation
        args.outdir = os.path.abspath(args.outdir)
    
    if args.mode=='dir':
        #to prevent trouble later with abspath detection in default filename generation
        #outdirs = []
        #for outdir in args.outdir:
        #    if outdir:
        #        outdir = os.path.abspath(outdir)
        #    outdirs.append(outdir)
        
        if 'call' in args.functions:
            call_output = []
            for indir,outdir,adaptertrimmed_dir,hits_dir,CDR3_dir in zip(
                                                                  args.input,
                                                                  args.outdir,
                                                                  args.adaptertrimmed_dir,
                                                                  args.hits_dir,
                                                                  args.CDR3_dir
                                                                  ):
                call_output.append(
                                reptools.call_dir(
                                    indir = indir,
                                    dbs = args.dbs,
                                    genedictfile = args.genedict,
                                    db_dir = args.db_dir,
                                    outdir = outdir,
                                    notrim = args.notrim,
                                    adaptertrimmed_dir = adaptertrimmed_dir,
                                    hits_dir = hits_dir,
                                    CDR3_dir = CDR3_dir,
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
                                )
            args.input = call_output
        
        if 'denoise' in args.functions:
            denoise_output = []
            for indir,outdir,FASTQout_dir in zip(args.input,args.outdir,args.denoised_dir):
                denoise_output.append(
                                      reptools.denoise_dir(
                                                    indir,
                                                    outdir = outdir,
                                                    FASTQout_dir = FASTQout_dir,
                                                    subs = args.denoise_substitution,
                                                    indels = args.denoise_indel,
                                                    deambig = args.denoise_gene_segments,
                                                    weight_by_qual = True,
                                                    threshold = args.threshold,
                                                    indel_threshold = args.indel_threshold,
                                                    overwrite = args.overwrite
                                                   )
                                       )
            
            
            args.input = denoise_output
        
        if 'EEfilter' in args.functions:
            EE_output = []
            for indir,outdir,FASTQout_dir in zip(args.input,args.outdir,args.filtCDR3_dir):
                EE_output.append(
                                reptools.EEfilter_dir(
                                                indir,
                                                outdir = outdir,
                                                FASTQout_dir = FASTQout_dir,
                                                maxee = args.CDR3maxee,
                                                overwrite = args.overwrite
                                            )
                                )
            args.input = EE_output
        
        if 'VDJtools' in args.functions:
            for indir,outdir,VDJout_dir in zip(args.input,args.outdir,args.VDJtools_dir):
                reptools.make_VDJtools_dir(
                                    indir,
                                    outdir = outdir,
                                    VDJout_dir = VDJout_dir,
                                    genes = args.labels,
                                    emptycols = ['D'],
                                    overwrite = args.overwrite
                                    )
        
        #remove temporary output files
        reptools.clean_up_dirs(args.tempdirs)
    
    if args.mode=='file':
        if not args.outdir:
            args.outdir = os.path.split(os.path.realpath(os.path.expanduser(args.input[0])))[0]
        
        reptools.ensure_dir(args.outdir)
        toremove = []
        in1_stem = os.path.split(os.path.splitext(args.input[0])[0])[1]
        stem =  in1_stem.split(args.pairsuffixes[0])[0]
        if 'call' in args.functions:
            if args.db_dir:
                args.dbs = {gene:os.path.join(args.db_dir,args.dbs[gene]) for gene in args.dbs}       
            #make default file paths:
            in2_stem = os.path.split(os.path.splitext(args.input[1])[0])[1]
            if args.adaptertrimmed_filestem is None:
                adaptertrimmed1 = build_fn(None,'{}_trimmed.fastq'.format(in1_stem),args.outdir,maketemp=False)
                adaptertrimmed2 = build_fn(None,'{}_trimmed.fastq'.format(in2_stem),args.outdir,maketemp=False)
            else:
                adaptertrimmed1 = '{}{}.fastq'.format(args.adaptertrimmed_filestem, args.pairsuffixes[0])
                adaptertrimmed2 = '{}{}.fastq'.format(args.adaptertrimmed_filestem, args.pairsuffixes[1])
                adaptertrimmed1 = build_fn(adaptertrimmed1,None,args.outdir) #add the outdir if supplied is not abs path
                adaptertrimmed2 = build_fn(adaptertrimmed2,None,args.outdir)
            if args.hits_filestem is None:
                outF = build_fn(None,'{}_hits_F.fastq'.format(stem),args.outdir,maketemp=False)
                outR = build_fn(None,'{}_hits_R.fastq'.format(stem),args.outdir,maketemp=False)
            else:
                outF = '{}_F.fastq'.format(args.hits_filestem)
                outR = '{}_R.fastq'.format(args.hits_filestem)
                outF = build_fn(outF,None,args.outdir) #add the outdir if supplied is not abs path
                outR = build_fn(outR,None,args.outdir)
            
            args.CDR3_file = build_fn(args.CDR3_file,'{}_CDR3.fastq'.format(stem),args.outdir)
        
        if 'denoise' in args.functions:
            args.denoised_file = build_fn(args.denoised_file,'{}_denoised.fastq'.format(stem),args.outdir)
            if args.seq2clone_file:
                args.seq2clone_file = build_fn(args.seq2clone_file,'{}_seq2clone.tsv'.format(stem),args.outdir)
        
        if 'EEfilter' in args.functions:
            args.filtCDR3_file = build_fn(args.filtCDR3_file,'{}_EEfilt.fastq'.format(stem),args.outdir)
        
        if 'VDJtools' in args.functions:
            args.VDJtools_file = build_fn(args.VDJtools_file,'{}_VDJtools.tab'.format(stem),args.outdir)
        
        ###############
        #now run things
        ###############
        if 'call' in args.functions:
            args.input[0] = reptools.call_filepair(
                                                    in1 = args.input[0],
                                                    in2 = args.input[1],
                                                    dbs = args.dbs,
                                                    genedictfile = args.genedict,
                                                    noCDR3 = args.noCDR3,
                                                    notrim = args.notrim,
                                                    adaptertrimmed1 = adaptertrimmed1,
                                                    adaptertrimmed2 = adaptertrimmed2,
                                                    outF = outF,
                                                    outR = outR,
                                                    outCDR3 = args.CDR3_file,
                                                    Cmiss1 = False, Cmiss2 = False,
                                                    JmissF = False, JmissR = False,
                                                    VmissF = False, VmissR = False,
                                                    outVsegF = False, outVsegR = False,
                                                    title_split = args.title_split,
                                                    overwrite = args.overwrite,
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
                                                    clusterID_position = args.clusterID_position,
                                                    blastdb_version=args.blastdb_version
                                                    )
            if not args.input[0]: return(None)
        if 'denoise' in args.functions:
            args.input[0],_,changelogs = reptools.denoise_file(
                                                        args.input[0],
                                                        outdir = args.outdir,
                                                        weight_by_qual = True,
                                                        subs = args.denoise_substitution,
                                                        indels = args.denoise_indel,
                                                        deambig = args.denoise_gene_segments,
                                                        FASTQout_fn = args.denoised_file,
                                                        overwrite = args.overwrite,
                                                        threshold = args.threshold,
                                                        indel_threshold = args.indel_threshold
                                                    )
            if not args.input[0]: return(None)
        if 'EEfilter' in args.functions:
            args.input[0] = reptools.EEfilter_file(
                                            args.input[0],
                                            FASTQout = args.filtCDR3_file,
                                            maxee = args.CDR3maxee
                                        )
            if not args.input[0]: return(None)
        if 'VDJtools' in args.functions:
            reptools.make_VDJtools(
                                   CDR3file = args.input[0],
                                   outfile = args.VDJtools_file,
                                   genes = args.labels,
                                   emptycols = ['D']
                                  )
            if not args.input[0]: return(None)
        if args.seq2clone_file:
            reptools.make_seq2clone(
                                    CDR3file = args.input[0],
                                    cluster_files=changelogs,
                                    outfile=args.seq2clone_file
                                   )
            if not args.input[0]: return(None)
        
        #print(args)

        #if 'call' in arg.functions and outCDR3 is None:
        #    os.remove() #remove temporary file
        
        #if 'denoise' and args.filtCDR3_file is None:
        #    os.remove() #remove temporary file
        
#if __name__ == "__main__": 
#    main()
    
    
#py -3 "C:\Users\Stephen Preston\Sync\TCRSeq\reptools\expt_code\developing_denoiser\cli_testing.py" -h

