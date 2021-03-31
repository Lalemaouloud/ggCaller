#imports
import argparse
import sys
import json
from Bio.Seq import Seq
import ggCaller_cpp
from ggCaller.graph_traversal import *
from multiprocessing import Pool
from functools import partial


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
                    help='List of reference genomes (one file path per line). ')
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
    IO.add_argument('--options.repeat',
                    action="store_true",
                    default=False,
                    help='Enable traversal of nodes mulitple times. '
                         '[Default = False] ')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads for FMIndexing '
                         '[Default = 1] ')
    IO.add_argument('--out',
                    default='calls.fasta',
                    help='options.out FASTA file containing ORF sequences. ')
    return parser.parse_args()

def main():
    options = get_options()

    # parse command line arguments
    options.out = options.out

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

    # if build graph specified, build graph and then call ORFs
    if (options.graph != None) and (options.colours != None) and (options.refs == None) and (options.reads == None):
        graph_tuple = ggCaller_cpp.index_existing(options.graph, options.colours, stop_codons_for, stop_codons_rev,
                                                  options.threads, True)
    # if refs file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs != None) and (options.reads == None):
        graph_tuple = ggCaller_cpp.index_build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                               options.threads, True, options.no_write_graph)
    # if reads file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs == None) and (options.reads != None):
        graph_tuple = ggCaller_cpp.index_build(options.reads, options.kmer, stop_codons_for, stop_codons_rev,
                                               options.threads, False, options.no_write_graph)
    # if both reads and refs file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs != None) and (options.reads != None):
        graph_tuple = ggCaller_cpp.index_build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                               options.threads, False, options.no_write_graph, options.reads)
    else:
        print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
              "- Bifrost GFA and Bifrost colours file\n"
              "- List of reference files\n"
              "- List of read files\n"
              "- A list of reference files and a list of read files.")
        sys.exit(1)

    # unpack ORF pair into overlap dictionary and list for gene scoring
    graph_vector, node_colour_vector, input_colours, nb_colours, overlap = graph_tuple

    # create list for high scoring ORFs to return
    true_genes = {}

    # calculate ORFs within graph
    for colour_ID, node_set in enumerate(node_colour_vector):
        print("Colour number: " + str(colour_ID))
        ORF_overlap_dict, ORF_vector = ggCaller_cpp.calculate_ORFs(graph_vector, colour_ID, node_set, options.repeat,
                                                                   overlap, options.path, options.not_ref,
                                                                   options.no_filter,
                                                                   stop_codons_for, start_codons, options.orf,
                                                                   options.maxoverlap, options.no_write_idx,
                                                                   input_colours[colour_ID])

        # if not filter specified, just append directly to true_genes
        if options.no_filter:
            for ORFNodeVector in ORF_vector:
                gene = generate_seq(graph_vector, ORFNodeVector[0], ORFNodeVector[1], overlap)
                if gene not in true_genes:
                    # create tuple to hold ORF sequence, colours and graph traversal information
                    empty_colours_list = ["0"] * nb_colours
                    true_genes[gene] = (empty_colours_list, ORFNodeVector)
                # update colours with current colour_ID
                true_genes[gene][0][colour_ID] = "1"
        else:
            # calculate scores for genes
            print("Calculating gene scores...")
            ORF_score_dict = score_genes(ORF_vector, graph_vector, options.min_orf_score, overlap, options.threads)

            # determine highest scoring genes
            print("Calculating highest scoring genes...")
            high_scoring_ORFs = call_true_genes(ORF_score_dict, ORF_overlap_dict, options.min_path_score)

            for ORF_id in high_scoring_ORFs:
                ORFNodeVector = ORF_vector[ORF_id]
                gene = generate_seq(graph_vector, ORFNodeVector[0], ORFNodeVector[1], overlap)
                if gene not in true_genes:
                    # create tuple to hold ORF sequence, colours and graph traversal information
                    empty_colours_list = ["0"] * nb_colours
                    true_genes[gene] = (empty_colours_list, ORFNodeVector)
                # update colours with current colour_ID
                true_genes[gene][0][colour_ID] = "1"

    print("Generating fasta file of gene calls...")
    # print options.out to file
    ORF_count = 1
    with open(options.out, "w") as f:
        for gene, info_pair in true_genes.items():
            colour_str = "".join(info_pair[0])
            f.write(">" + str(ORF_count) + "_" + str(colour_str) + "\n" + gene + "\n")
            ORF_count += 1

    print("Finished.")

    sys.exit(0)

if __name__ == '__main__':
    main()
