import argparse
import os
import sys
import reptools


def parse_args():
    parser = argparse.ArgumentParser(
                         description='main parser',
                         usage='reptools [function] -i <input> -o <outdir> [...]\n',
                         prog='reptools',
                         formatter_class=argparse.RawDescriptionHelpFormatter
                                )

    parser.add_argument(
        '-V',
        '--version',
        action='version',
        version=f"reptools ({reptools.__version__})"
    )

    reptools_functions = ['blast', 'sw', 'extract_cdr3', 'annotate', 'qualfilter', 'format_chains']
    requiredFuns = parser.add_argument_group(
        'functions',
        'Function(s) - choose one:\n'
        '         blast  Call blastn and return an json file giving info on hits, and\n'
        '                (optionally) files containing seqs with and without hits.\n'
        '            sw  Run a smith-waterman search and return an json file giving info\n'
        '                on hits, and (optionally) files containing seqs with and without\n'
        '                hits.\n'
        '  extract_cdr3  Extract CDR3 from sequences using sw alignment results.\n'
        '      annotate  Annotate sequence files with gene segment ids from json files\n'
        '    qualfilter  Filter sequences according to mean Q score.\n'
        ' format_chains  Process data from microfluidics workflow to combine alpha and\n'
        '                beta chain info.'
        )
    requiredFuns.add_argument('function', choices=reptools_functions)

    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument(
        '-i',
        '--infiles',
        nargs='+',
        required=True,
        help="Input sequence file(s).  May generally be FASTA or FASTQ (except\n"
             "for quality filtering!)."
    )
    requiredNamed.add_argument(
        '-o',
        '--outdir',
        required=True,
        help="A path to write output files to.  Will be created if it doesn't exist."
    )

    optionalArgs = parser.add_argument_group('Optional arguments')
    optionalArgs.add_argument(
        '-j',
        '--out_json',
        default=False,
        help="An output json file to hold info on search results."
    )
    optionalArgs.add_argument(
        '--v_hits',
        default=False,
        help="A json file giving info on V gene hits (produced by the blast or sw function)."
    )
    optionalArgs.add_argument(
        '--j_hits',
        default=False,
        help="A json file giving info on J gene hits (produced by the blast or sw function)."
    )
    optionalArgs.add_argument(
        '-d',
        '--database',
        help='To call genes, the path to a fasta database must be supplied'
        )
    optionalArgs.add_argument(
        '-g',
        '--gene',
        nargs='+',
        default=False,
        help='For blast and smith-waterman searches, default settings can be provided if a\n'
             'target gene is set: available options are currently V and J.  Note that the\n'
             'default V settings for Smith-Waterman alignment use only the 3 prime 80\n'
             'nucleotides in the supplied database sequences. If gene is not set, a config\n'
             'file must be supplied using --config (-c).'
        )
    optionalArgs.add_argument(
        '-c',
        '--config',
        nargs='+',
        default=False,
        help='For blast and smith-waterman searches, a config file in json format giving\n'
             'search settings.  Alternatively, --gene (-g) may be used to invoke appropriate\n'
             'defaults.'
        )
    optionalArgs.add_argument(
        '--miss',
        nargs='+',
        default=False,
        help="A path to write files containing sequences without an alignment hit to.\n"
             "Will be created if it doesn't exist."
        )
    optionalArgs.add_argument(
        '-t',
        '--threads',
        default=False,
        type=int,
        help='Number of threads to use for alignment (defaults to number of available CPUs).'
        )
    optionalArgs.add_argument(
        '--titleSplit',
        default=' ',
        dest='title_split',
        help='Substring seperating sequence ID from metadata in fastq title lines.\n'
             'Defaults to a space.'
        )
    optionalArgs.add_argument(
        '-g',
        '--genedict',
        default=False,
        help='To call genes, the path to a csv gene dictionary must be supplied'
        )
    optionalArgs.add_argument('--overwrite', action='store_true')
    optionalArgs.add_argument(
        '--aligner_path',
        nargs='+',
        default=False,
        help=('Paths to aligners and associated files.\n'
              'Defaults are blastn=blastn swipe=swipe makeblastdb=makeblastdb bbduk=bbduk.sh')
        )
    optionalArgs.add_argument(
        '--titleFormat',
        default='clusterID|readID',
        help=('Paired file integrity checks depend on the fastq title line format. The default\n'
              'is "clusterID|readID", which adequately represents illumina output.  For\n'
              'mixcr-tools synthetic data, use "readID|clusterID".  Use | to represent the\n'
              'title seperator, the actual identity of which is set with --titleSplit (which\n'
              'defaults to a space, but should be set to "|" for mixcr tools).\n'
              'used by: reptools blast/sw/\n'
              'default=5')
    )
    optionalArgs.add_argument(
        '--blastdb_version',
        type=int,
        default=5,
        choices={4, 5},
        help=('Version 5 BLAST databases are used for BLAST calls by default, but a bug in\n'
              'BLAST+ can lead to this producing an out of memory error on some systems.  If\n'
              'this occurs, set "--blastdb_version 4".\n'
              'used by: reptools blast\n'
              'default=5')
    )
    optionalArgs.add_argument(
        '--mean_qual',
        type=int,
        default=20,
        help=('The minimum mean qual value to accept.\n'
              'used by: reptools qualfilter\n'
              'default=20')
        )

    # read args
    args = parser.parse_args()

    # function selection logic - check for required arguments
    if args.functions in ['blast', 'sw']:
        for argname in ['out_json', 'database']:
            if not locals()[f"args.{argname}"]:
                print(f"To perform a blast or smith-waterman alignment, --{argname} must be set.")
                sys.exit(1)
        if not (args.gene or args.config):
            print("To perform a blast or smith-waterman alignment, either --gene or --config must\n"
                  "be set.")
            sys.exit(1)
    elif args.functions == 'extract_cdr3':
        for argname in ['v_hits', 'j_hits', 'genedict']:
            if not locals()[f"args.{argname}"]:
                print(f"To extract CDR3 sequences, --{argname} must be set.")
                sys.exit(1)
    elif args.functions == 'annotate':
        if not (args.v_hits or args.j_hits):
            print("To annotate with gene segment IDs, at least one of , either --v_hits or "
                  "--j_hits\n"
                  "must be set.")
            sys.exit(1)
    elif args.functions == 'qualfilter':
        pass  # TODO
    elif args.functions == 'format_chains':
        pass  # TODO
    else:
        # show list of functions
        parser.print_help()

    return args


def main():
    args = parse_args()


if __name__ == '__main__':
    sys.exit(main())
