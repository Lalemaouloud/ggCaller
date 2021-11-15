# imports
import argparse
from ggCaller.graph_traversal import *
import ggCaller_cpp
from functools import partial
# from memory_profiler import profile
from balrog.__main__ import *
from ggCaller.shared_memory import *
import tqdm
from panaroo_runner.set_default_args import *
from panaroo_runner.__main__ import run_panaroo
from panaroo_runner.annotate import *
import ast
import tempfile


def get_options():
    description = 'Generates ORFs from a Bifrost graph.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='ggcaller')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--graph',
                    default=None,
                    help='Bifrost GFA file generated by Bifrost build. ')
    IO.add_argument('--colours',
                    default=None,
                    help='Bifrost colours file generated by Bifrost build.  ')
    IO.add_argument('--not-ref',
                    action="store_false",
                    default=True,
                    help='If using existing graph, was not graph built exclusively with assembled genomes.  '
                         '[Default = False] ')
    IO.add_argument('--refs',
                    default=None,
                    help='List of reference genomes (one file path per line). ')
    IO.add_argument('--reads',
                    default=None,
                    help='List of read files (one file path per line). ')
    IO.add_argument('--codons',
                    default=None,
                    help='JSON file containing start and stop codon sequences. ')
    IO.add_argument('--kmer',
                    type=int,
                    default=31,
                    help='K-mer size used in Bifrost build (bp). '
                         '[Default = 31] ')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads to use. '
                         '[Default = 1] ')
    IO.add_argument('--out',
                    default='ggCaller_output',
                    help='Output directory ')
    Settings = parser.add_argument_group('Cut-off settings')
    Settings.add_argument('--max-path-length',
                          type=int,
                          default=10000,
                          help='Maximum path length during ORF finding (bp). '
                               '[Default = 10000] ')
    Settings.add_argument('--min-orf-length',
                          type=int,
                          default=90,
                          help='Minimum ORF length to return (bp). '
                               '[Default = 90] ')
    Settings.add_argument('--max-ORF-overlap',
                          type=int,
                          default=60,
                          help='Maximum overlap allowed between overlapping ORFs. '
                               '[Default = 60] ')
    Settings.add_argument('--min-path-score',
                          type=int,
                          default=100,
                          help='Minimum total Balrog score for a path of ORFs to be returned. '
                               '[Default = 100] ')
    Settings.add_argument('--min-orf-score',
                          type=int,
                          default=100,
                          help='Minimum individual Balrog score for an ORF to be returned. '
                               '[Default = 100] ')
    Settings.add_argument('--max-orf-orf-distance',
                          type=int,
                          default=10000,
                          help='Maximum distance for graph traversal during ORF connection (bp). '
                               '[Default = 10000] ')
    Settings.add_argument('--identity-cutoff',
                          type=float,
                          default=0.98,
                          help='Minimum identity at amino acid level between two ORFs for clustering. '
                               '[Default = 0.98] ')
    Settings.add_argument('--len-diff-cutoff',
                          type=float,
                          default=0.98,
                          help='Minimum ratio of length between two ORFs for clustering.  '
                               '[Default = 0.98] ')
    Algorithm = parser.add_argument_group('Settings to avoid/include algorithms')
    Algorithm.add_argument('--no-filter',
                           action="store_true",
                           default=False,
                           help='Do not filter ORF calls using Balrog. Will return all ORF calls. '
                                '[Default = False] ')
    Algorithm.add_argument('--no-write-idx',
                           action="store_false",
                           default=True,
                           help='Do not write FMIndexes to file. '
                                '[Default = False] ')
    Algorithm.add_argument('--no-write-graph',
                           action="store_false",
                           default=True,
                           help='Do not write Bifrost GFA and colours to file. '
                                '[Default = False] ')
    Algorithm.add_argument('--repeat',
                           action="store_true",
                           default=False,
                           help='Enable traversal of nodes multiple times. '
                                '[Default = False] ')
    Algorithm.add_argument('--no-clustering',
                           action="store_true",
                           default=False,
                           help='Do not cluster ORFs. '
                                '[Default = False] ')
    Panaroo_mode_opts = parser.add_argument_group('Mode')

    Panaroo_mode_opts.add_argument(
        "--clean-mode",
        dest="mode",
        help=
        ('''R|The stringency mode at which to run panaroo. Must be one of 'strict',\
    'moderate' or 'sensitive'. Each of these modes can be fine tuned using the\
     additional parameters in the 'Graph correction' section.

    strict: 
    Requires fairly strong evidence (present in  at least 5%% of genomes)\
     to keep likely contaminant genes. Will remove genes that are refound more often than\
     they were called originally.

    moderate: 
    Requires moderate evidence (present in  at least 1%% of genomes)\
     to keep likely contaminant genes. Keeps genes that are refound more often than\
     they were called originally.

    sensitive: 
    Does not delete any genes and only performes merge and refinding\
     operations. Useful if rare plasmids are of interest as these are often hard to\
     disguish from contamination. Results will likely include  higher number of\
     spurious annotations.'''),
        choices=['strict', 'moderate', 'sensitive'],
        required=False)

    Panaroo_mode_opts.add_argument(
        "--remove-invalid-genes",
        dest="filter_invalid",
        action='store_true',
        default=False,
        help=(
                "removes annotations that do not conform to the expected Prokka" +
                " format such as those including premature stop codons."))

    Panaroo_matching = parser.add_argument_group('Matching')
    Panaroo_matching.add_argument(
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        type=float)
    Panaroo_matching.add_argument("--merge_paralogs",
                                  dest="merge_paralogs",
                                  help="don't split paralogs",
                                  action='store_true',
                                  default=False)
    Panaroo_matching.add_argument("--annotation",
                                  dest="annotate",
                                  help="Annotate genes using diamond (fast) or diamond and HMMscan (sensitive)."
                                       "If not specified, no annotation done",
                                  choices=["fast", "sensitive"],
                                  default=None)
    Panaroo_matching.add_argument("--diamonddb",
                                  dest="annotation_db",
                                  help="Diamond database. Defaults are 'Bacteria' or 'Viruses'. Can also "
                                       "specify path to fasta file for custom database generation",
                                  default="Bacteria")
    Panaroo_matching.add_argument("--hmmdb",
                                  dest="hmm_db",
                                  help="HMMER hmm profile file. Default is Uniprot HAMAP. Can also"
                                       "specify path to pre-built hmm profile file generated using hmmbuild",
                                  default="default")
    Panaroo_matching.add_argument("--evalue",
                                  dest="evalue",
                                  help="Maximum e-value to return for DIAMOND and HMMER searches during annotation",
                                  default=0.001,
                                  type=float)

    Panaroo_refind = parser.add_argument_group('Refind')
    Panaroo_refind.add_argument(
        "--search_radius",
        dest="search_radius",
        help=("the distance in nucleotides surronding the " +
              "neighbour of an accessory gene in which to search for it"),
        default=5000,
        type=int)
    Panaroo_refind.add_argument(
        "--refind_prop_match",
        dest="refind_prop_match",
        help=("the proportion of an accessory gene that must " +
              "be found in order to consider it a match"),
        default=0.2,
        type=float)

    Panaroo_graph = parser.add_argument_group('Graph correction')

    Panaroo_graph.add_argument(
        "--min_trailing_support",
        dest="min_trailing_support",
        help=("minimum cluster size to keep a gene called at the " +
              "end of a contig"),
        type=int)
    Panaroo_graph.add_argument(
        "--trailing_recursive",
        dest="trailing_recursive",
        help=("number of times to perform recursive trimming of low support " +
              "nodes near the end of contigs"),
        type=int)
    Panaroo_graph.add_argument(
        "--edge_support_threshold",
        dest="edge_support_threshold",
        help=(
                "minimum support required to keep an edge that has been flagged" +
                " as a possible mis-assembly"),
        type=float)
    Panaroo_graph.add_argument(
        "--length_outlier_support_proportion",
        dest="length_outlier_support_proportion",
        help=
        ("proportion of genomes supporting a gene with a length more " +
         "than 1.5x outside the interquatile range for genes in the same cluster"
         +
         " (default=0.01). Genes failing this test will be re-annotated at the "
         + "shorter length"),
        type=float,
        default=0.1)
    Panaroo_graph.add_argument(
        "--remove_by_consensus",
        dest="remove_by_consensus",
        type=ast.literal_eval,
        choices=[True, False],
        help=
        ("if a gene is called in the same region with similar sequence a minority "
         + "of the time, remove it. One of 'True' or 'False'"),
        default=None)
    Panaroo_graph.add_argument(
        "--high_var_flag",
        dest="cycle_threshold_min",
        help=(
                "minimum number of nested cycles to call a highly variable gene " +
                "region (default = 5)."),
        type=int,
        default=5)
    Panaroo_graph.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=("minimum edge support required to call structural variants" +
              " in the presence/absence sv file"),
        type=int)
    Panaroo_graph.add_argument(
        "--all_seq_in_graph",
        dest="all_seq_in_graph",
        help=("Retains all DNA sequence for each gene cluster in the graph " +
              "output. Off by default as it uses a large amount of space."),
        action='store_true',
        default=False)
    Panaroo_graph.add_argument(
        "--no_clean_edges",
        dest="clean_edges",
        help=("Turn off edge filtering in the final output graph."),
        action='store_false',
        default=True)

    Panaroo_core = parser.add_argument_group('Gene alignment')
    Panaroo_core.add_argument(
        "--alignment",
        dest="aln",
        help=("Output alignments of core genes or all genes. Options are" +
              " 'core' and 'pan'. Default: 'None'"),
        type=str,
        choices=['core', 'pan'],
        default=None)
    Panaroo_core.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'ref' for reference-guided MSA and 'def' for default standard MSA",
        type=str,
        choices=['def', 'ref'],
        default="def")
    Panaroo_core.add_argument("--core_threshold",
                              dest="core",
                              help="Core-genome sample threshold (default=0.95)",
                              type=float,
                              default=0.95)

    # Other options
    parser.add_argument("--codon-table",
                        dest="table",
                        help="the codon table to use for translation (default=11)",
                        type=int,
                        default=11)
    parser.add_argument("--quiet",
                        dest="verbose",
                        help="suppress additional output",
                        action='store_false',
                        default=True)

    return parser.parse_args()

def main():
    # parse command line arguments for ggCaller
    options = get_options()

    # determine if references/assemblies present
    if (options.refs != None and options.reads == None) or (options.graph != None and options.colours != None
                                                            and options.not_ref):
        is_ref = True
    else:
        is_ref = False

    # define start/stop codons
    if options.codons != None:
        with open(options.codons, "r") as json_file:
            try:
                data = json.load(json_file)
                start_codons = data["codons"]["start"]
                stop_codons_for = data["codons"]["stop"]
                stop_codons_rev = [str((Seq(i)).reverse_complement()) for i in stop_codons_for]
            except:
                print("Please specify codons in the format shown in codons.json.")
                sys.exit(1)
    else:
        start_codons = ["ATG", "GTG", "TTG"]
        stop_codons_for = ["TAA", "TGA", "TAG"]
        stop_codons_rev = ["TTA", "TCA", "CTA"]

    # initialise graph
    graph = ggCaller_cpp.Graph()

    # if build graph specified, build graph and then call ORFs
    if (options.graph != None) and (options.colours != None) and (options.refs == None) and (options.reads == None):
        graph_tuple = graph.read(options.graph, options.colours, stop_codons_for, stop_codons_rev,
                                 options.threads)
    # if refs file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs != None) and (options.reads == None):
        graph_tuple = graph.build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                  options.threads, True, options.no_write_graph, "NA")
    # if reads file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs == None) and (options.reads != None):
        graph_tuple = graph.build(options.reads, options.kmer, stop_codons_for, stop_codons_rev,
                                  options.threads, False, options.no_write_graph, "NA")
    # if both reads and refs file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs != None) and (options.reads != None):
        graph_tuple = graph.build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                  options.threads, False, options.no_write_graph, options.reads)
    else:
        print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
              "- Bifrost GFA and Bifrost colours file\n"
              "- List of reference files\n"
              "- List of read files\n"
              "- A list of reference files and a list of read files.")
        sys.exit(1)

    # unpack ORF pair into overlap dictionary and list for gene scoring
    node_colour_vector, input_colours, nb_colours, overlap = graph_tuple

    # set rest of panaroo arguments
    options = set_default_args(options, nb_colours)

    # check diamond and HMMER are installed correctly
    check_diamond_install()
    check_HMMER_install()

    annotation_db = options.annotation_db
    hmm_db = options.hmm_db
    if options.annotate is not None:
        # unpack annotation database
        if annotation_db == "Bacteria" or annotation_db == "Viruses":
            db_id = annotation_db
            db_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "db", "diamond")
            annotation_db = os.path.join(db_dir, annotation_db)

            if not os.path.exists(annotation_db):
                print("Unzipping protein annotation file...")
                tar = tarfile.open(annotation_db + ".tar.gz", mode="r:gz")
                tar.extractall(db_dir)
                tar.close()

            annotation_db = os.path.join(annotation_db, db_id + ".dmnd")

        # if custom annotation database specified, then create diamond db if not present already
        else:
            if ".dmnd" not in annotation_db:
                print("Generating diamond index...")
                annotation_db = generate_diamond_index(annotation_db)

        # set-up hmm_db
        if hmm_db == "default":
            db_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "db", "hmm")
            hmm_db = os.path.join(db_dir, "HAMAP.hmm")

        if not os.path.exists(hmm_db + ".h3f"):
            print("Generating HMMER index...")
            generate_HMMER_index(hmm_db)

    # create directory if it isn't present already
    if not os.path.exists(options.out):
        os.mkdir(options.out)

    # make sure trailing forward slash is present
    output_dir = os.path.join(options.out, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=output_dir), "")

    # create numpy arrays for shared memory
    total_arr = np.array([graph])

    # load balrog models if required
    if not options.no_filter:
        print("Loading gene models...")
        model, model_tis = load_gene_models()

    else:
        model, model_tis = None, None

    # intiialise results dictionaries and lists
    high_scoring_ORFs = {}
    high_scoring_ORF_edges = {}
    cluster_id_list = None
    cluster_dict = None

    # use shared memory to generate graph vector
    print("Generating high scoring ORF calls per colour...")

    # set number of threads for graphtool and pytorch to 1
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(1)

    torch.set_num_threads(1)

    with SharedMemoryManager() as smm:
        # generate shared numpy arrays
        total_arr = np.append(total_arr, [[model], [model_tis]])
        array_shd, array_shd_tup = generate_shared_mem_array(total_arr, smm)

        # run run_calculate_ORFs with multithreading
        with Pool(processes=options.threads) as pool:
            for colour_ID, gene_dict, ORF_edges in tqdm.tqdm(pool.imap(
                    partial(run_calculate_ORFs, shd_arr_tup=array_shd_tup, repeat=options.repeat, overlap=overlap,
                            max_path_length=options.max_path_length, is_ref=is_ref,
                            no_filter=options.no_filter,
                            stop_codons_for=stop_codons_for, start_codons=start_codons,
                            min_ORF_length=options.min_orf_length,
                            max_ORF_overlap=options.max_ORF_overlap, minimum_ORF_score=options.min_orf_score,
                            minimum_path_score=options.min_path_score, write_idx=options.no_write_idx,
                            input_colours=input_colours, max_orf_orf_distance=options.max_orf_orf_distance),
                    enumerate(node_colour_vector)), total=nb_colours):
                high_scoring_ORFs[colour_ID] = gene_dict
                high_scoring_ORF_edges[colour_ID] = ORF_edges

            # generate ORF clusters
            if not options.no_clustering:
                print("Generating clusters of high-scoring ORFs...")
                cluster_id_list, cluster_dict = graph.generate_clusters(high_scoring_ORFs,
                                                                        overlap,
                                                                        options.identity_cutoff,
                                                                        options.len_diff_cutoff)

                run_panaroo(pool, array_shd_tup, high_scoring_ORFs, high_scoring_ORF_edges, cluster_id_list,
                            cluster_dict, overlap, input_colours, output_dir, temp_dir, options.verbose,
                            options.threads,
                            options.length_outlier_support_proportion, options.identity_cutoff,
                            options.family_threshold, options.min_trailing_support, options.trailing_recursive,
                            options.clean_edges, options.edge_support_threshold, options.merge_paralogs, options.aln,
                            options.alr, options.core, options.min_edge_support_sv, options.all_seq_in_graph, is_ref,
                            options.no_write_idx, overlap + 1, options.repeat, options.remove_by_consensus,
                            options.search_radius, options.refind_prop_match, options.annotate, options.evalue,
                            annotation_db, hmm_db)

    print("Finished.")

    sys.exit(0)


if __name__ == '__main__':
    main()
