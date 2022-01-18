import graph_tool.all as gt
from balrog.__main__ import *
from ggCaller.shared_memory import *
import _pickle as cPickle

def range_overlapping(x, y):
    if x.start == x.stop or y.start == y.stop:
        return False
    return x.start <= y.stop and y.start <= x.stop

def search_graph(graph, graphfile, coloursfile, queryfile, objects_dir, output_dir, query_id, num_threads):
    # check if objects_dir present, if not exit
    if not os.path.exists(objects_dir):
        print("Please specify a ggc_data directory")
        sys.exit(1)

    objects_dir = os.path.join(objects_dir, "")

    # read in and parse sequences in queryfile
    query_vec = []
    with open(queryfile, "r") as f:
        for line in f:
            if line[0] != ">":
                query_vec.append(line.strip())

    p0 = psutil.Process()

    print("pre-graph read-in: Perc: " + str(p0.memory_percent()) + " full: " + str(p0.memory_info()))

    # read in graph object and high_scoring ORFs and node_index
    graph.data_in(objects_dir + "ggc_graph.dat", graphfile, coloursfile, num_threads)

    print("post-graph read-in: Perc: " + str(p0.memory_percent()) + " full: " + str(p0.memory_info()))

    # query the sequences in the graph
    print("Querying unitigs in graph...")
    input_colours, kmer, query_nodes = graph.search_graph(graphfile, coloursfile, query_vec, query_id, num_threads)

    # parse isolate names
    isolate_names = [
        os.path.splitext(os.path.basename(x))[0] for x in input_colours
    ]

    with open(objects_dir + "high_scoring_orfs.dat", "rb") as input_file:
        high_scoring_ORFs = cPickle.load(input_file)

    print("post-ORFs read-in: Perc: " + str(p0.memory_percent()) + " full: " + str(p0.memory_info()))

    with open(objects_dir + "node_index.dat", "rb") as input_file:
        node_index = cPickle.load(input_file)

    print("post-kmers read-in: Perc: " + str(p0.memory_percent()) + " full: " + str(p0.memory_info()))

    outfile = output_dir + "matched_queries.fasta"
    print("Matching overlapping ORFs...")
    with open(outfile, "w") as f:
        for i in range(len(query_nodes)):
            query_set = set()
            for node in query_nodes[i]:
                query_set.update(node_index[node])
                for ORF in query_set:
                    split_ID = ORF.split("_")
                    colour = int(split_ID[0])
                    ORF_ID = int(split_ID[1])
                    fasta_ID = isolate_names[colour] + "_" + str(ORF_ID).zfill(5)
                    ORF_info = high_scoring_ORFs[colour][ORF_ID]
                    seq = graph.generate_sequence(ORF_info[0], ORF_info[1], kmer - 1)
                    # add annotation if available
                    if len(ORF_info) == 8 or ORF_ID < 0:
                        ORF_annotation = ORF_info[-1]
                        f.write(
                            ">" + fasta_ID + " " + ORF_annotation[-1] + " QUERY=" + query_vec[i] + "\n" + seq + "\n")
                    else:
                        f.write(">" + fasta_ID + " QUERY=" + query_vec[i] + "\n" + seq + "\n")

    return

# @profile
def traverse_components(component, tc, component_list, edge_weights, minimum_path_score):
    # initilise high scoring ORF set to return
    high_scoring_ORFs = set()

    # generate subgraph view
    u = gt.GraphView(tc, vfilt=component_list == component)

    # initialise list of high scoring ORFs and their associated score for single component
    high_scoring_ORFs_temp = []

    high_score_temp = 0

    # iterate over edges, determine which are source and sink nodes for connected components
    vertices = u.get_vertices()
    in_degs = u.get_in_degrees(u.get_vertices())
    out_degs = u.get_out_degrees((u.get_vertices()))

    # calculate start nodes and add to list, if in-degree is 0
    start_vertices = [vertices[i] for i in range(len(in_degs)) if in_degs[i] == 0]

    # calculate end nodes and add to list, if in-degree is 0
    end_vertices = [vertices[i] for i in range(len(out_degs)) if out_degs[i] == 0]

    # iterate over start and stop vertices
    for start in start_vertices:
        # add the score of the first node
        start_score = u.vertex_properties["score"][start]

        # check if start vertex is lone node
        if start in end_vertices:
            vertex_list = [start]
            # replace high_score_ORFs with new high scoring ORF path
            if start_score > high_score_temp:
                high_score_temp = start_score
                high_scoring_ORFs_temp = vertex_list
        else:
            for end in end_vertices:
                score = start_score
                vertex_list, edge_list = gt.shortest_path(u, start, end, weights=edge_weights,
                                                          negative_weights=True, dag=True)
                # ensure a path has been found. If not, pass.
                if edge_list:
                    for e in edge_list:
                        # add the inverse of the score to get the true score
                        score -= edge_weights[e]

                    # replace high_score_ORFs with new high scoring ORF path
                    if score > high_score_temp:
                        high_score_temp = score
                        high_scoring_ORFs_temp = vertex_list

    # for highest scoring path, see if greater than cut-off, if so, add high scoring ORFs to set
    if high_score_temp >= minimum_path_score:
        ORF_ID_list = tuple([u.vertex_properties["ID"][node] for node in high_scoring_ORFs_temp])
        high_scoring_ORFs.add(ORF_ID_list)

    return high_scoring_ORFs


#@profile
def call_true_genes(ORF_score_array, ORF_overlap_dict, minimum_path_score):
    # initilise high scoring ORF set to return
    high_scoring_ORFs_all = set()

    # create empty, directed graph
    g = gt.Graph()

    # create a dictionaries assign each ORF an index in the graph to vertices
    ORF_index = {}

    # add vertexes to graph, store ORF information in ORF_index
    for index, score in np.ndenumerate(ORF_score_array):
        if score != 0:
            ORF_ID = index[0]
            v = g.add_vertex()
            ORF_index[ORF_ID] = int(v)

    # add edges and edge weights between connected ORFs using ORF_overlap_dict. ORF1 is sink, ORF2 is source
    for ORF1, overlap_dict in ORF_overlap_dict.items():
        # check that ORF1 nodes exist in graph, may have been removed due to low score
        if ORF1 in ORF_index:
            for ORF2 in overlap_dict.keys():
                # check that ORF2 nodes exist in graph as before
                if ORF2 in ORF_index:
                    # add new edge between the two ORFs, where ORF2 is the source and ORF1 is the sink
                    e = g.add_edge(g.vertex(ORF_index[ORF2]), g.vertex(ORF_index[ORF1]))

    # determine if cycles present. If so, break them by removing edge before repeated node and re-test
    cycle = True
    try:
        circuit = next(gt.all_circuits(g))
    except StopIteration:
        cycle = False

    while cycle:
        end_cycle = circuit[-1]
        start_cycle = circuit[0]
        e = g.edge(end_cycle, start_cycle)
        g.remove_edge(e)
        try:
            circuit = next(gt.all_circuits(g))
        except StopIteration:
            break

    # generate a transative closure of the graph to add all directed edges and add vertex properties
    tc = gt.transitive_closure(g)

    # clear original graph
    g.clear()

    # get label components of tc
    components = gt.label_components(tc, directed=False)[0].a

    # create vertex property map to store node IDs and scores
    vertex_ID = tc.new_vertex_property("int")
    vertex_score = tc.new_vertex_property("double")
    for ORF_ID, index in ORF_index.items():
        vertex_ID[tc.vertex(index)] = ORF_ID
        vertex_score[tc.vertex(index)] = ORF_score_array[ORF_ID]

    # add vertex IDs and scores as internal property of graph
    tc.vertex_properties["ID"] = vertex_ID
    tc.vertex_properties["score"] = vertex_score

    # create edge weights property
    edge_weights = tc.new_edge_property("double")

    # create incompatible list for edge iteration
    incomp_edges = []

    # iterate over edges and assign weights, and identify invalid overlaps to remove
    for e in tc.edges():
        # parse sink and source nodes, order same as before
        ORF1 = tc.vertex_properties["ID"][e.target()]
        ORF2 = tc.vertex_properties["ID"][e.source()]

        # parse sink ORF score calculated by Balrog
        ORF_score = ORF_score_array[ORF1]

        # check if edges are present by overlap detection. If not, set edge weight as ORF1 score
        if ORF1 not in ORF_overlap_dict:
            edge_weights[e] = -(ORF_score)
            continue
        elif ORF2 not in ORF_overlap_dict[ORF1]:
            edge_weights[e] = -(ORF_score)
            continue
        # parse overlap info, calculate edge weight
        overlap_type, abs_overlap = ORF_overlap_dict[ORF1][ORF2]

        # set penalty to 0 for "n" or "w" or "i" overlaps
        penalty = 0
        if overlap_type == "i":
            # add incompatible ORFs as tuple, added as source (1) and sink(2)
            incomp_edges.append(e)
        elif overlap_type == "u":
            penalty = unidirectional_penalty_per_base * abs_overlap
        elif overlap_type == "c":
            penalty = convergent_penalty_per_base * abs_overlap
        elif overlap_type == "d":
            penalty = divergent_penalty_per_base * abs_overlap
        # set ORF score to negative so that shortest algorithm can be applied
        edge_weights[e] = -(ORF_score - penalty)

    # add edge_weights as internal property
    tc.edge_properties["weights"] = edge_weights

    # iterate over incompatible edges to remove them
    for e in incomp_edges:
        tc.remove_edge(e)

    # iterate over components, find highest scoring path within component with multiprocessing to determine geniest path through components
    for component in set(components):
        high_scoring_ORFs = traverse_components(component, tc, components, edge_weights, minimum_path_score)
        high_scoring_ORFs_all.update(high_scoring_ORFs)

    return high_scoring_ORFs_all


# @profile
def run_calculate_ORFs(node_set_tuple, shd_arr_tup, repeat, overlap, max_path_length, is_ref, no_filter,
                       stop_codons_for, start_codons, min_ORF_length, max_ORF_overlap, minimum_ORF_score,
                       minimum_path_score, write_idx, input_colours, max_orf_orf_distance):
    # unpack tuple
    colour_ID, node_set = node_set_tuple

    # p_id = "p" + str(colour_ID)
    # p = psutil.Process()

    # print(p_id + " pre-traversal: Perc: " + str(p.memory_percent()) + " full: " + str(p.memory_info()))

    # load shared memory items
    existing_shm = shared_memory.SharedMemory(name=shd_arr_tup.name)
    shd_arr = np.ndarray(shd_arr_tup.shape, dtype=shd_arr_tup.dtype, buffer=existing_shm.buf)

    # print("Finding ORFs: " + str(colour_ID))

    # determine all ORFs in Bifrost graph
    ORF_overlap_dict, ORF_vector = shd_arr[0].findORFs(colour_ID, node_set, repeat,
                                                       overlap, max_path_length, is_ref, no_filter,
                                                       stop_codons_for, start_codons, min_ORF_length,
                                                       max_ORF_overlap, write_idx, input_colours[colour_ID])

    # print(p_id + " post-traversal: Perc: " + str(p.memory_percent()) + " full: " + str(p.memory_info()))

    # initialise return dictionaries
    gene_dict = {}

    # if no filter specified, just copy ORF_vector to gene_dict with dictionary comprehension
    if no_filter:
        gene_dict = {i: ORF_vector[i] for i in range(0, len(ORF_vector))}
        edge_list = set([(i,) for i in range(0, len(ORF_vector))])
    else:
        # print("Scoring ORFs: " + str(colour_ID))

        # calculate scores for genes
        ORF_score_array = score_genes(ORF_vector, shd_arr[0], minimum_ORF_score, overlap, shd_arr[1], shd_arr[2])

        # print(p_id + " post-scoring: Perc: " + str(p.memory_percent()) + " full: " + str(p.memory_info()))

        # print("Finding highest scoring paths: " + str(colour_ID))

        # determine highest scoring genes, stored in list of lists
        edge_list = call_true_genes(ORF_score_array, ORF_overlap_dict, minimum_path_score)

        # print(p_id + " post-high-scoring determination: Perc: " + str(p.memory_percent()) + " full: " + str(
        #     p.memory_info()))

        # generate a dictionary of all true gene info
        for entry in edge_list:
            for ORF in entry:
                gene_dict[ORF] = ORF_vector[ORF]

    # generate list of target ORFs, removing duplicates
    target_ORFs = list(set([x for f in edge_list for x in (f[0], f[-1])]))

    # print("Connecting ORFs: " + str(colour_ID))

    # determine next ORFs for each terminal ORF in edge_list
    next_ORFs = set(shd_arr[0].connect_ORFs(colour_ID, ORF_vector, target_ORFs, max_orf_orf_distance, is_ref))

    # determine redundant edges in high_scoring_ORFs
    redundant_edges = set([tuple(sorted([i[0], i[-1]])) for i in edge_list if len(i) > 1])

    # convert edge list into list
    edge_list = [i for i in edge_list]

    # remove redundant edges between high_scoring_ORFs and next_nodes
    for edge in next_ORFs:
        if edge not in redundant_edges:
            edge_list.append(edge)

    # print("Arranging highest scoring ORFs: " + str(colour_ID))

    # iterate over edge_list and append to a dictionary of high_scoring_ORF_edges for each ORF
    high_scoring_ORF_edges = {}
    for entry in edge_list:
        # work out last index of entry
        last_index = len(entry) - 1
        for i, ORF in enumerate(entry):
            # create new entry for current ORF
            if ORF not in high_scoring_ORF_edges:
                high_scoring_ORF_edges[ORF] = set()
            # if ORF not last in list, add next entry to set
            if (i != last_index):
                high_scoring_ORF_edges[ORF].add(entry[i + 1])

    # print("Finished:  " + str(colour_ID))
    # print(p_id + " post-ORF-calling: Perc: " + str(p.memory_percent()) + " full: " + str(p.memory_info()))

    return colour_ID, gene_dict, high_scoring_ORF_edges
