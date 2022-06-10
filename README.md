# ggCaller: a gene caller for Bifrost graphs

ggCaller traverses [Bifrost](https://github.com/pmelsted/bifrost) graphs constructed from bacterial genomes to identify putative protein coding sequences, known as open reading frames (ORFs). 

ggCaller incorporates [Balrog](https://github.com/salzberg-lab/Balrog) to filter ORFs to improve specificity of calls and [Panaroo](https://github.com/gtonkinhill/panaroo) for pangenome analysis and quality control.

## Installation

ggCaller is currently only available for Linux. If you are running Windows 10, it can be installed via the Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10))

### Installation via conda (recommended)

Install through [bioconda](http://bioconda.github.io/):

```conda install ggcaller```

If conda is not installed, first install [miniconda](https://docs.conda.io/en/latest/miniconda.html), then add the correct channels:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Installation from source
Required packages and versions can be found in ```environment.yml```. In addition, a C++17 compiler (e.g. gcc >=7.3) is required.

Once all required packages are installed, install ggCaller using:

```
git clone --recursive https://github.com/samhorsfield96/ggCaller 
python setup.py install
```

## Gene-call mode
In gene-call mode, ggCaller identifies ORFs within a Bifrost graph, filters them using BALROG and clusters them using Panaroo. 

ggCaller takes a list of fasta files (one file per line), or a Bifrost GFA file and Bifrost Colours file generated by ```Bifrost build```. See the [Bifrost](https://github.com/pmelsted/bifrost) repository for installation.

ggCaller additionally employs FMindexing using [SDSLv3](https://github.com/xxsds/sdsl-lite)  to remove artificial sequences generated by incorrectly phased nodes in the DBG. This is only employed for genomes designated as references. 

#### If Bifrost GFA and Colours file do not exist:

To build a new Bifrost graph using assembled genomes and/or reads, specify ONE OR BOTH:
- ```--refs <refs.txt>``` List of absolute paths to reference sequence fastas (one file per line)
- ```--reads <reads.txt>``` List of absolute paths to read fastas (one file per line)

Note: Ensure assembled genome files are exclusively passed to the ```--refs``` argument, and read files exclusively to the ```--reads```
argument. Bifrost uses kmer coverage filtering for read files to remove read errors, but does not do this for assembled genomes.

#### If Bifrost GFA and Colours file already exist:

To run ggCaller using an existing Bifrost GFA file and Colours file, specify BOTH:
- ```--graph <graph.gfa>``` Input GFA
- ```--colours <colours.bfg_colors>``` Input colours file

To inform ggCaller of which files are references, supply them as a list to the ```--refs``` argument as above.
- ```--refs <refs.txt>``` List of absolute paths to reference sequence fastas (one file per line)

If ```--refs``` is not specified, ggCaller assumes all sequences are assembled. To avoid this, i.e. when no sequences are assembled genomes, additionally specify:
- ```--not-ref```

Note: Ensure the sequences used to build the graph are in the same directories as when the graph was built.

#### Additional helpful arguments
- ```--kmer``` k-mer size for graph building. Only used for building of new graphs (default: 31 bp).
- ```--no-filter``` Do not conduct ORF filtering. ggCaller will return all ORFs present.
- ```--threads``` Number of threads (default: 1).
- ```--clean-mode {sensitive, moderate, strict}``` specify stringency for Panaroo quality control. See [Panaroo parameters](https://gtonkinhill.github.io/panaroo/#/gettingstarted/params) for details (default: 'sensitive')
- ```--annotation {fast, sensitive}``` annotate clusters. Either specify ```fast``` to use diamond, or ```sensitive``` to use diamond + HMMscan.
- ```--diamonddb``` specify path to diamond annotation database. Should be fasta format. Default is Bacterial dataset (```Bacteria```) downloaded with ggCaller. Can also specify Viral dataset (```Viruses```). Both default annotation datasets are from [Uniprot](https://www.uniprot.org/).
- ```--hmmdb``` specify HMMscan annotation database. Should be pre-trained HMMER profile database. Default is Pfam dataset from [Prokka](https://github.com/tseemann/prokka), downloaded with ggCaller.
- ```--alignment {core, pan}``` generate alignments and VCFs for core genome (```core```) or all clusters (```pan```).
- ```--aligner {def, ref}``` use with ```--alignment```, specify whether to align all genes in cluster together at once (```def```), or via reference guided approach (```ref```). ```ref``` is faster when more genomes are used to build the graph. 
- ```--out``` output directory (default: 'ggCaller_output')

## Examples
- Build graph using assembled genomes, using kmer size of 31 bp and strict clean mode. 

```ggcaller --refs refs.txt --kmer 31 --clean-mode strict --out output_dir```

- Build graph using reads using fast annotation using diamond. 

```ggcaller --reads reads.txt --annotation fast --out output_dir```

- Build graph using assembled genomes and reads, with reference-based pangenome alignment. 

```ggcaller --refs refs.txt --reads reads.txt --alignment pan --aligner ref --out output_dir```

- Use existing graph, specifying genomes which are assembled and no filtering.

```ggcaller --graph graph.gfa --colours colours.bfg_colours --refs refs.txt --no-filter --out output_dir```

- Use existing graph which was built using only reads.

```ggcaller --graph graph.gfa --colours colours.bfg_colours --not-ref --out output_dir```

Test data is available in the ```data``` directory.

### Outputs

ggCaller generates a number of outputs depending on the arguments specified.

#### Always generated
- Annotated gene calls in fasta format
- Gene presence/absence matrix in Panaroo, Roary and RTAB format
- Structural variant presence/absence matrix
- Pangenome reference, containing largest sequence for each gene cluster
- Pre/post filtered Panaroo graph
- Rarefaction curve
- Gene frequency plot
- Cluster size plot
- Neighbour joining tree based on gene presence/absence
- Summary statistics on gene frequency

#### if ```--annotation fast/sensitive``` specified
- Annotated gene calls in GFF format

#### if ```--alignment pan``` or  ```--alignment core``` specified
- Neighbour joining tree based on core genome alignment
- Core genome alignment

#### if ```--alignment pan``` specified
- Aligned sequences for each cluster (if ```--alignment pan``` specified)
- Variant calls for each cluster in VCF format

## Query mode
ggCaller can query pre-called genes within a graph. This requires ggCaller to have been run on a prior dataset, and the outputs saved.

### Example
First generate a set of gene calls and save them. This will generate a data directory within the output directory, ```output_dir/ggc_data```

```ggcaller --graph graph.gfa --colours colours.bfg_colours --save --out output_dir```

Then supply ggCaller with the graph and the data directory, along with the queries in fasta format, saving to same directory

```ggcaller --graph graph.gfa --colours colours.bfg_colours --data output_dir/ggc_data --query query.fasta --out output_dir```

### Output
ggCaller generates a fasta file detailing the genes overlapping with each query, and in which genomes they are found.

## All I/O options
```
Generates ORFs from a Bifrost graph.

optional arguments:
  -h, --help            show this help message and exit

Input/Output options:
  --graph GRAPH         Bifrost GFA file generated by Bifrost build.
  --colours COLOURS     Bifrost colours file generated by Bifrost build.
  --not-ref             If using existing graph, was not graph built exclusively with assembled genomes. [Default = False]
  --refs REFS           List of reference genomes (one file path per line).
  --reads READS         List of read files (one file path per line).
  --query QUERY         List of unitig sequences to query (either FASTA or one sequence per line)
  --codons CODONS       JSON file containing start and stop codon sequences.
  --kmer KMER           K-mer size used in Bifrost build (bp). [Default = 31]
  --save                Save graph objects for sequence querying. [Default = False]
  --data DATA           Directory containing data from previous ggCaller run generated via "--save"
  --all-seq-in-graph    Retains all DNA sequence for each gene cluster in the Panaroo graph output. Off by default as it uses a large amount of space.
  --out OUT             Output directory

ggCaller traversal and gene-calling cut-off settings:
  --max-path-length MAX_PATH_LENGTH
                        Maximum path length during ORF finding (bp). [Default = 20000]
  --min-orf-length MIN_ORF_LENGTH
                        Minimum ORF length to return (bp). [Default = 90]
  --max-ORF-overlap MAX_ORF_OVERLAP
                        Maximum overlap allowed between overlapping ORFs. [Default = 60]
  --min-path-score MIN_PATH_SCORE
                        Minimum total Balrog score for a path of ORFs to be returned. [Default = 100]
  --min-orf-score MIN_ORF_SCORE
                        Minimum individual Balrog score for an ORF to be returned. [Default = 100]
  --max-orf-orf-distance MAX_ORF_ORF_DISTANCE
                        Maximum distance for graph traversal during ORF connection (bp). [Default = 10000]
  --query-id QUERY_ID   Ratio of query-kmers to required to match in graph. [Default = 0.8]

Settings to avoid/include algorithms:
  --no-filter           Do not filter ORF calls using Balrog. Will return all ORF calls. [Default = False]
  --no-write-idx        Do not write FMIndexes to file. [Default = False]
  --no-write-graph      Do not write Bifrost GFA and colours to file. [Default = False]
  --repeat              Enable traversal of nodes multiple times. [Default = False]
  --no-clustering       Do not cluster ORFs. [Default = False]
  --no-refind           Do not refind uncalled genes [Default = False]

Gene clustering options.:
  --identity-cutoff IDENTITY_CUTOFF
                        Minimum identity at amino acid level between two ORFs for clustering. [Default = 0.98]
  --len-diff-cutoff LEN_DIFF_CUTOFF
                        Minimum ratio of length between two ORFs for clustering. [Default = 0.98]
  --family-threshold FAMILY_THRESHOLD
                        protein family sequence identity threshold (default=0.7)
  --merge-paralogs      don't split paralogs

Panaroo run mode options:
  --clean-mode {strict,moderate,sensitive}
                        R|The stringency mode at which to run panaroo. Must be one of 'strict', 'moderate' or 'sensitive'. Each of these modes can be fine tuned using the additional parameters in the 'Graph
                        correction' section. strict: Requires fairly strong evidence (present in at least 5% of genomes) to keep likely contaminant genes. Will remove genes that are refound more often than
                        they were called originally. moderate: Requires moderate evidence (present in at least 1% of genomes) to keep likely contaminant genes. Keeps genes that are refound more often than they
                        were called originally. sensitive: Does not delete any genes and only performes merge and refinding operations. Useful if rare plasmids are of interest as these are often hard to
                        disguish from contamination. Results will likely include higher number of spurious annotations.

Panaroo gene cluster annotation options:
  --annotation {fast,sensitive,ultrasensitive}
                        Annotate genes using diamond default (fast), diamond sensitive (sensitive) or diamond and HMMscan (ultrasensitive).If not specified, no annotation done
  --diamonddb ANNOTATION_DB
                        Diamond database. Defaults are 'Bacteria' or 'Viruses'. Can also specify path to fasta file for custom database generation
  --hmmdb HMM_DB        HMMER hmm profile file. Default is Uniprot HAMAP. Can alsospecify path to pre-built hmm profile file generated using hmmbuild
  --evalue EVALUE       Maximum e-value to return for DIAMOND and HMMER searches during annotation
  --truncation-threshold TRUNCATION_THRESHOLD
                        Sequences in a gene family cluster below this proportion of the length of thecentroid will be annotated as 'potential pseudogene'

Panaroo gene-refinding options:
  --search-radius SEARCH_RADIUS
                        the distance in nucleotides surronding the neighbour of an accessory gene in which to search for it
  --refind-prop-match REFIND_PROP_MATCH
                        the proportion of an accessory gene that must be found in order to consider it a match

Panaroo graph correction stringency options:
  --remove-invalid-genes
                        removes annotations that do not conform to the expected Prokka format such as those including premature stop codons.
  --min-trailing-support MIN_TRAILING_SUPPORT
                        minimum cluster size to keep a gene called at the end of a contig
  --trailing-recursive TRAILING_RECURSIVE
                        number of times to perform recursive trimming of low support nodes near the end of contigs
  --edge-support-threshold EDGE_SUPPORT_THRESHOLD
                        minimum support required to keep an edge that has been flagged as a possible mis-assembly
  --length-outlier-support-proportion LENGTH_OUTLIER_SUPPORT_PROPORTION
                        proportion of genomes supporting a gene with a length more than 1.5x outside the interquatile range for genes in the same cluster (default=0.01). Genes failing this test will be re-
                        annotated at the shorter length
  --remove-by-consensus {True,False}
                        if a gene is called in the same region with similar sequence a minority of the time, remove it. One of 'True' or 'False'
  --high-var-flag CYCLE_THRESHOLD_MIN
                        minimum number of nested cycles to call a highly variable gene region (default = 5).
  --min-edge-support-sv MIN_EDGE_SUPPORT_SV
                        minimum edge support required to call structural variants in the presence/absence sv file
  --no-clean-edges      Turn off edge filtering in the final output graph.

Gene alignment options:
  --alignment {core,pan}
                        Output alignments of core genes or all genes. Options are 'core' and 'pan'. Default: 'None'
  --aligner {def,ref}   Specify an aligner. Options:'ref' for reference-guided MSA and 'def' for default standard MSA
  --core-threshold CORE
                        Core-genome sample threshold (default=0.95)
  --no-variants         Do not call variants using SNP-sites after alignment.
  --ignore-pseduogenes  Ignore ORFs annotated as 'potential pseudogenes' in alignment

Misc. options:
  --quiet               suppress additional output
  --threads THREADS     Number of threads to use. [Default = 1]
  --version, -v         show program's version number and exit
```

## Citation

If you use this code, please cite:

Bifrost: 
Holley, G., Melsted, P. Bifrost: highly parallel construction and indexing of colored and compacted de Bruijn graphs. Genome Biol 21, 249 (2020). https://doi.org/10.1186/s13059-020-02135-8

Balrog:
Sommer MJ, Salzberg SL. Balrog: A universal protein model for prokaryotic gene prediction. PLoS Comput Biol 17(2): e1008727 (2021). https://doi.org/10.1371/journal.pcbi.1008727

Panaroo:
Tonkin-Hill, G., MacAlasdair, N., Ruis, C. et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21, 180 (2020). https://doi.org/10.1186/s13059-020-02090-4

SDSL v3:
[Succinct Data Structure Library 3.0](https://github.com/xxsds/sdsl-lite)

Edlib:
Šošić, M., Šikić, M. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. Bioinformatics 33, 9, (2017). https://doi.org/10.1093/bioinformatics/btw753

Eigen v3:
Guennebaud, G., Jacob, B. et al. Eigen v3 (2010). http://eigen.tuxfamily.org