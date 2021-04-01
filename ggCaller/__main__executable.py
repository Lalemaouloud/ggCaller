# imports
import argparse
import sys
import json
from Bio.Seq import Seq
import ggCaller_cpp
from ggCaller.graph_traversal import *
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from memory_profiler import profile

def get_options():
    description = 'Generates ORFs from a Bifrost graph.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='ggcaller')

    IO = parser.add_argument_group('Input/output')
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
                    help='List of reference genomes (one file path per line). ')
    IO.add_argument('--codons',
                    default=None,
                    help='JSON file containing start and stop codon sequences. ')
    IO.add_argument('--kmer',
                    default=31,
                    help='K-mer size used in Bifrost build (bp). '
                         '[Default = 31] ')
    IO.add_argument('--path',
                    default=10000,
                    help='Maximum path length during traversal (bp). '
                         '[Default = 10000] ')
    IO.add_argument('--orf',
                    default=90,
                    help='Minimum ORF length to return (bp). '
                         '[Default = 90] ')
    IO.add_argument('--maxoverlap',
                    default=60,
                    help='Maximum overlap allowed between overlapping ORFs. '
                         '[Default = 60] ')
    IO.add_argument('--min-path-score',
                    default=100,
                    help='Minimum total score for a path of ORFs to be returned. '
                         '[Default = 100] ')
    IO.add_argument('--min-orf-score',
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
                    help='Enable traversal of nodes mulitple times. '
                         '[Default = False] ')
    IO.add_argument('--threads',
                    default=1,
                    help='Number of threads for FMIndexing '
                         '[Default = 1] ')
    IO.add_argument('--out',
                    default='calls.fasta',
                    help='Output FASTA file containing ORF sequences. ')
    return parser.parse_args()


@profile
def main():
    # options = get_options()
    #
    # # parse command line arguments
    # graph_file = options.graph
    # colours_file = options.colours
    # is_ref = bool(options.not_ref)
    # refs_file = options.refs
    # reads_file = options.reads
    # codons = options.codons
    # ksize = int(options.kmer)
    # max_path_length = int(options.path)
    # min_ORF_length = int(options.orf)
    # max_ORF_overlap = int(options.maxoverlap)
    # minimum_path_score = int(options.min_path_score)
    # minimum_ORF_score = int(options.min_orf_score)
    # no_filter = bool(options.no_filter)
    # write_idx = bool(options.no_write_idx)
    # write_graph = bool(options.no_write_graph)
    # repeat = bool(options.repeat)
    # num_threads = int(options.threads)
    # output = options.out
    #
    # # define start/stop codons
    # if codons != None:
    #     with open(codons, "r") as json_file:
    #         try:
    #             data = json.load(json_file)
    #             start_codons = data["codons"]["start"]
    #             stop_codon_for = data["codons"]["stop"]
    #             stop_codon_rev = [str((Seq(i)).reverse_complement()) for i in stop_codon_for]
    #         except:
    #             print("Please specify codons in the format shown in codons.json.")
    #             sys.exit(1)
    # else:
    #     start_codons = ["ATG", "GTG", "TTG"]
    #     stop_codon_for = ["TAA", "TGA", "TAG"]
    #     stop_codon_rev = ["TTA", "TCA", "CTA"]

    start_codons = ["ATG", "GTG", "TTG"]
    stop_codons_for = ["TAA", "TGA", "TAG"]
    stop_codons_rev = ["TTA", "TCA", "CTA"]

    output = "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/group3_capsular_fa_list_integer_paths.fasta"
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

    num_threads = 4

    graph_tuple = ggCaller_cpp.index_existing(
        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.gfa",
        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.bfg_colors",
        stop_codons_for, stop_codons_rev, num_threads, True)

    # unpack ORF pair into overlap dictionary and list for gene scoring
    graph_vector, node_colour_vector, input_colours, nb_colours, overlap = graph_tuple

    if not no_filter:
        print("Loading gene scoring models...")
        model, model_tis, aa_kmer_set = load_models(num_threads)
    else:
        model, model_tis, aa_kmer_set = None

    # run run_calculate_ORFs with multithreading
    true_genes = {}
    print("Generating high scoring ORF calls...")
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for colour_ID, col_true_genes in executor.map(
                partial(run_calculate_ORFs, graph_vector=graph_vector, repeat=repeat, overlap=overlap,
                        max_path_length=max_path_length,
                        is_ref=is_ref, no_filter=no_filter, stop_codons_for=stop_codons_for, start_codons=start_codons,
                        min_ORF_length=min_ORF_length,
                        max_ORF_overlap=max_ORF_overlap, minimum_ORF_score=minimum_ORF_score,
                        minimum_path_score=minimum_path_score, write_idx=write_idx,
                        input_colours=input_colours, nb_colours=nb_colours, model=model, model_tis=model_tis,
                        aa_kmer_set=aa_kmer_set),
                enumerate(node_colour_vector)):
            for gene, ORFNodeVector in col_true_genes.items():
                gene = generate_seq(graph_vector, ORFNodeVector[0], ORFNodeVector[1], overlap)
                if gene not in true_genes:
                    # create tuple to hold ORF sequence, colours and graph traversal information
                    empty_colours_list = ["0"] * nb_colours
                    true_genes[gene] = (empty_colours_list, ORFNodeVector)
                # update colours with current colour_ID
                true_genes[gene][0][colour_ID] = "1"

    print("Generating fasta file of gene calls...")
    # print output to file
    ORF_count = 1
    with open(output, "w") as f:
        for gene, info_pair in true_genes.items():
            colour_str = "".join(info_pair[0])
            f.write(">" + str(ORF_count) + "_" + str(colour_str) + "\n" + gene + "\n")
            ORF_count += 1

    print("Finished.")

    sys.exit(0)


if __name__ == '__main__':
    main()
