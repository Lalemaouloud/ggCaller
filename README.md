# ggCaller: a gene caller for Bifrost graphs

ggCaller traverses Bifrost graphs constructed from bacterial genomes to identify putative protein coding sequences, known as open reading frames (ORFs). ggCaller incorporates [Balrog](https://github.com/salzberg-lab/Balrog) to filter ORFs to improve specificity of calls. 

Balrog determines a score for each ORF based on length, presence/absence of a translation initiation sites and amino acid composition using a convolutional neural network. ORF scores are then used to determine compatible high scoring ORFs based on their overlap.

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

Requirements:
-   ```python3```
-   ```biopython```
-   ```SeQan3```
-   ```bifrost```
-   ```pybind11```
-   ```cmake```
-   ```openMP```
-   ```pthreads```
-   ```zlib```
-   ```graph-tool```
-   ```pytorch```
-   ```tqdm```
-   a C++17 compiler (e.g. gcc >=7.3)

```
git clone --recursive https://github.com/samhorsfield96/ggCaller 
python setup.py install
```

## Input files for ggCaller

ggCaller can take a list of fasta files (one file per line), or a Bifrost GFA file and Bifrost Colours file generated by ```Bifrost build```. See the Bifrost repository for installation (https://github.com/pmelsted/bifrost).

## Usage

To run ggCaller from the command line, pass arguments specifying input/output and parameters for traversal.

ggCaller additionally employs FMindexing to remove artificial sequences generated by incorrectly phased nodes in the DBG. This is only employed when the input sequences are exclusively assembled genomes. 

#### If Bifrost GFA and Colours file already exist:

To run ggCaller using an existing Bifrost GFA file and Colours file, specify BOTH:
- ```--graph <graph.gfa>``` Input GFA
- ```--colours <colours.bfg_colors>``` Input colours file

To disable filtering of artificially generated sequences via FMindexing (if sequences used to generate the GFA and Colours files are not exclusively assembled genomes), additionally specify:
- ```--not-ref```

Note: Ensure the sequences used to build the graph are in the same directories as when the graph was built.

#### If Bifrost GFA and Colours file do not exist:

To build a new Bifrost graph using assembled genomes and/or reads, specify ONE OR BOTH:
- ```--refs <refs.txt>``` Input GFA
- ```--reads <reads.txt>``` Input colours file

Note: Ensure assembled genome files are exclusively passed to the ```--refs``` argument, and read files exclusively to the ```--reads```
argument. Bifrost uses kmer coverage filtering for read files to remove read errors, but does not do this for assembled genomes.

#### Additional arguments
- ```--kmer``` k-mer size for graph building. Only used for building of new graphs (default: 31 bp)
- ```--path``` maximum path length for traversal in base-pairs (default: 10000 bp)
- ```--orf``` minimum ORF length in base-pairs (default: 90 bp)
- ```--codons``` JSON file specifying START and STOP codons. See ```codons.json``` for example of format.
- ```--no-write-idx``` Don't write FMIndexes to file, otherwise in same directory as fasta files.
- ```--no-write-graph``` Don't write Bifrost GFA or Colours to file, otherwise in working directory.
- ```--repeat``` Enables traversal of nodes more than once, to detect genes with repeats.
- ```--maxoverlap``` Maximim overlap allowed between to ORFs (default: 60 bp)
- ```--min-orf-score``` Minimum individual ORF score allowed during filtering (default: 100)
- ```--min-path-score``` Minimum path score allowed for consecutive ORFs during overlap analysis in filtering (default: 100)
- ```--no-filter``` Do not conduct ORF filtering. ggCaller will return all ORFs present.
- ```--threads``` Number of threads (default: 1)
- ```--out``` output file in FASTA format (default: 'calls.fasta')

## Examples
- Build graph using assembled genomes, using kmer size of 31 bp. 

```ggcaller --refs refs.txt --kmer 31 --out calls.fasta```

- Build graph using reads, using kmer size of 31 bp. 

```ggcaller --reads reads.txt --kmer 31 --out calls.fasta```

- Build graph using assembled genomes and reads, using kmer size of 31 bp. 

```ggcaller --refs refs.txt --reads reads.txt --kmer 31 --out calls.fasta```

- Use existing graph which was built using assembled genomes only.

```ggcaller --graph graph.gfa --colours colours.bfg_colours --out calls.fasta```

- Use existing graph which was built using reads or assembled genomes and reads.

```ggcaller --graph graph.gfa --colours colours.bfg_colours --not-ref --out calls.fasta```

Test data is available in the ```data``` directory.

## Interpreting output FASTA

The output ```calls.fasta``` contains ORF sequences. Headers start with a unique ORF identifier (e.g. 1, 2, 3) followed by a colours array (e.g. 11100). The colours array describes the presence/absence of the ORF in the source genomes in the same order as they are passed in the ```refs.txt```/```reads.txt``` files. If both refs and reads are passed, then the ordering is reads followed by refs.

```
>1_11100
TTGTCATTTTATTGGTTTATTTTCTTAGCTTTGTTAGAGAGACAGAACTTGAACGTTCTTCTAG
>2_00111
ATGGGATTAATGCTGTCTTTATGGGTGTTGGTGGCAGTTTTGATGTATTATCAGGACACATTAAACGAGCTCCATTATGGATGCAAAAATTGA
>3_11100
TTGTTCAAAGGTGGTGTTACGATTTCAAGAACTCCTCTCAGTTCTGAGGACACGGTAATGATTGATGCGATAG
```

## Citation

If you use this code, please cite:

Bifrost: 
Holley G., Melsted, P. Bifrost – Highly parallel construction and indexing of colored and compacted de Bruijn graphs. (2019) bioRxiv 695338. doi: https://doi.org/10.1101/695338

SeQan3: 
Reinert, K. et al. The SeqAn C++ template library for efficient sequence analysis: A resource for programmers. (2017) Journal of biotechnology, 261, 157-168 doi: https://doi.org/10.1016/j.jbiotec.2017.07.017

Balrog:
Sommer, M.J. & Salzberg, S.L. Balrog: A universal protein model for prokaryotic gene prediction (2020) BioRxiv 2020.09.06.285304 doi: https://doi.org/10.1101/2020.09.06.285304