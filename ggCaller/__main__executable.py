# imports
import argparse
import json
from ggCaller.graph_traversal import *
import ggCaller_cpp
from functools import partial
# from memory_profiler import profile
from balrog.__main__ import *
from ggCaller.shared_memory import *
from panaroo_runner.__main__ import run_panaroo
import math


def get_options():
    description = 'Generates ORFs from a Bifrost graph.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='ggcaller')

    IO = parser.add_argument_group('Input/options.out')
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
    IO.add_argument('--path',
                    type=int,
                    default=10000,
                    help='Maximum path length during traversal (bp). '
                         '[Default = 10000] ')
    IO.add_argument('--orf',
                    type=int,
                    default=90,
                    help='Minimum ORF length to return (bp). '
                         '[Default = 90] ')
    IO.add_argument('--maxoverlap',
                    type=int,
                    default=60,
                    help='Maximum overlap allowed between overlapping ORFs. '
                         '[Default = 60] ')
    IO.add_argument('--min-path-score',
                    type=int,
                    default=100,
                    help='Minimum total score for a path of ORFs to be returned. '
                         '[Default = 100] ')
    IO.add_argument('--min-orf-score',
                    type=int,
                    default=100,
                    help='Minimum individual score for an ORF to be returned. '
                         '[Default = 100] ')
    IO.add_argument('--no-filter',
                    action="store_true",
                    default=False,
                    help='Do not filter ORF calls using Balrog. Will return all ORF calls. '
                         '[Default = False] ')
    IO.add_argument('--no-write-idx',
                    action="store_false",
                    default=True,
                    help='Do not write FMIndexes to file. '
                         '[Default = False] ')
    IO.add_argument('--no-write-graph',
                    action="store_false",
                    default=True,
                    help='Do not write Bifrost GFA and colours to file. '
                         '[Default = False] ')
    IO.add_argument('--repeat',
                    action="store_true",
                    default=False,
                    help='Enable traversal of nodes multiple times. '
                         '[Default = False] ')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads to use. '
                         '[Default = 1] ')
    IO.add_argument('--out',
                    default='calls.fasta',
                    help='FASTA file containing ORF sequences. ')
    return parser.parse_args()


def main():
    start_codons = ["ATG", "GTG", "TTG"]
    stop_codons_for = ["TAA", "TGA", "TAG"]
    stop_codons_rev = ["TTA", "TCA", "CTA"]

    # set mimimum path score
    minimum_path_score = 100
    minimum_ORF_score = 100
    no_filter = False
    repeat = False
    max_path_length = 10000
    is_ref = True
    min_ORF_length = 90
    max_ORF_overlap = 60
    write_idx = True
    write_graph = True
    identity_cutoff = 0.98
    len_diff_cutoff = 0.98
    max_orf_orf_distance = 10000
    cluster_ORFs = True

    num_threads = 1

    # graph_tuple = ggCaller_cpp.index_existing(
    #     "/mnt/c/Users/sth19/Documents/PhD/Experiments/gene_caller_comparions/ggCaller_results/clique_119_230_372_list.gfa",
    #     "/mnt/c/Users/sth19/Documents/PhD/Experiments/gene_caller_comparions/ggCaller_results/clique_119_230_372_list.bfg_colors",
    #     stop_codons_for, stop_codons_rev, num_threads, True)

    graph = ggCaller_cpp.Graph()

    # graph_tuple = graph.build(
    #     "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group1_capsular_fa_list.txt",
    #     31, stop_codons_for, stop_codons_rev, num_threads, is_ref, write_graph, "NA")

    graph_tuple = graph.read(
        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.gfa",
        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.bfg_colors",
        stop_codons_for, stop_codons_rev, num_threads, is_ref)

    # unpack ORF pair into overlap dictionary and list for gene scoring
    node_colour_vector, input_colours, nb_colours, overlap = graph_tuple

    # panaroo options
    out_dir = "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/panaroo_temp"
    verbose = True
    length_outlier_support_proportion = 0.1
    family_threshold = 0.7
    min_trailing_support = max(2, math.ceil(0.05 * nb_colours))
    trailing_recursive = 99999999
    clean_edges = True
    edge_support_threshold = max(2, math.ceil(0.01 * nb_colours))
    merge_paralogs = True
    aln = "core"
    alr = "mafft"
    core = 0.95
    min_edge_support_sv = max(2, math.ceil(0.01 * nb_colours))
    all_seq_in_graph = False

    # create numpy arrays for shared memory
    total_arr = np.array([graph])

    # load balrog models if required
    if not no_filter:
        print("Loading gene models...")
        model, model_tis = load_gene_models()

    else:
        model, model_tis = None, None, None

    # intiialise results dictionaries and lists
    high_scoring_ORFs = {}
    high_scoring_ORF_edges = {}
    cluster_id_list = None
    cluster_dict = None

    # use shared memory to generate graph vector
    print("Generating high scoring ORF calls...")

    # set number of threads for graphtool and pytorch to 1
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(1)

    torch.set_num_threads(1)

    with SharedMemoryManager() as smm:
        # generate shared numpy arrays
        total_arr = np.append(total_arr, [[model], [model_tis]])
        array_shd, array_shd_tup = generate_shared_mem_array(total_arr, smm)

        # run run_calculate_ORFs with multithreading
        # with Pool(processes=num_threads) as pool:
        #     for colour_ID, col_true_genes in pool.map(
        #             partial(run_calculate_ORFs, shd_arr_tup=array_shd_tup, repeat=repeat, overlap=overlap,
        #                     max_path_length=max_path_length, is_ref=is_ref, no_filter=no_filter,
        #                     stop_codons_for=stop_codons_for, start_codons=start_codons, min_ORF_length=min_ORF_length,
        #                     max_ORF_overlap=max_ORF_overlap, minimum_ORF_score=minimum_ORF_score,
        #                     minimum_path_score=minimum_path_score, write_idx=write_idx,
        #                     input_colours=input_colours,
        #                     aa_kmer_set=aa_kmer_set),
        #             enumerate(node_colour_vector)):
        #         # iterate over entries in col_true_genes to generate the sequences
        #         for ORFNodeVector in col_true_genes:
        #             gene = graph.generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)
        #             if gene not in true_genes:
        #                 # create tuple to hold ORF sequence, colours and graph traversal information
        #                 empty_colours_list = ["0"] * nb_colours
        #                 true_genes[gene] = (empty_colours_list, ORFNodeVector)
        #             # update colours with current colour_ID
        #             true_genes[gene][0][colour_ID] = "1"

        # run run_calculate_ORFs with multithreading
        # with Pool(processes=num_threads) as pool:
        #     for colour_ID in pool.map(
        #             partial(run_calculate_ORFs, shd_arr_tup=array_shd_tup, repeat=repeat, overlap=overlap,
        #                     max_path_length=max_path_length, is_ref=is_ref, no_filter=no_filter,
        #                     stop_codons_for=stop_codons_for, start_codons=start_codons, min_ORF_length=min_ORF_length,
        #                     max_ORF_overlap=max_ORF_overlap, minimum_ORF_score=minimum_ORF_score,
        #                     minimum_path_score=minimum_path_score, write_idx=write_idx,
        #                     input_colours=input_colours,
        #                     aa_kmer_set=aa_kmer_set),
        #             enumerate(node_colour_vector)):
        #         # iterate over entries in col_true_genes to generate the sequences
        #         true_genes[colour_ID] = "done"

        for colour_tuple in enumerate(node_colour_vector):
            colour_ID, gene_dict, ORF_edges = run_calculate_ORFs(colour_tuple, shd_arr_tup=array_shd_tup, repeat=repeat,
                                                                 overlap=overlap,
                                                                 max_path_length=max_path_length, is_ref=is_ref,
                                                                 no_filter=no_filter,
                                                                 stop_codons_for=stop_codons_for,
                                                                 start_codons=start_codons,
                                                                 min_ORF_length=min_ORF_length,
                                                                 max_ORF_overlap=max_ORF_overlap,
                                                                 minimum_ORF_score=minimum_ORF_score,
                                                                 minimum_path_score=minimum_path_score,
                                                                 write_idx=write_idx,
                                                                 input_colours=input_colours,
                                                                 max_orf_orf_distance=max_orf_orf_distance)
            # iterate over entries in col_true_genes to generate the sequences
            # true_genes[colour_ID] = {}
            high_scoring_ORFs[colour_ID] = gene_dict
            # high_scoring_ORF_edges[colour_ID] = {}
            high_scoring_ORF_edges[colour_ID] = ORF_edges

        # cluster ORFs
        if cluster_ORFs is True:
            cluster_id_list, cluster_dict = graph.generate_clusters(high_scoring_ORFs, overlap, identity_cutoff,
                                                                    len_diff_cutoff)

            run_panaroo(graph, high_scoring_ORFs, high_scoring_ORF_edges, cluster_id_list, cluster_dict, overlap,
                        input_colours, out_dir, verbose, num_threads,
                        length_outlier_support_proportion, identity_cutoff, len_diff_cutoff,
                        family_threshold, min_trailing_support, trailing_recursive,
                        clean_edges, edge_support_threshold, merge_paralogs, aln,
                        alr, core, min_edge_support_sv, all_seq_in_graph)

    # print("Generating fasta file of gene calls...")
    # # print output to file
    # ORF_count = 1
    # with open(output, "w") as f:
    #     for gene, info_pair in true_genes.items():
    #         colour_str = "".join(info_pair[0])
    #         f.write(">" + str(ORF_count) + "_" + str(colour_str) + "\n" + gene + "\n")
    #         ORF_count += 1

    print("Finished.")

    sys.exit(0)


if __name__ == '__main__':
    main()
