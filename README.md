SlidingWindowTrees
=========================


Caveats
----------
This code is a work in progress and is not yet production quality.  Please ensure you have backups of your data before you use this code.


Supported Platform Stack
----------
* **Sequencing**:  MiSeq Illumina paired-end reads with fragment sizes such that the mates overlap.
* **HyPhy**:  HYPHYMP v2.2  (MPI HyPhy)
* **FastTree**: FastTreeMP v2.1  (Multithreading via OpenMP)
* **OS**:  Linux 64-bit
* **Python**:  python v2.7.x 64-bit


Prerequisites
--------------------
* Ensure you have python installed
* Download the OpenMP version of [FastTree](www.microbesonline.org/fasttree): [FastTreeMP](http://www.microbesonline.org/fasttree/FastTreeMP)
* Download and install multithreaded HyPhy v2.2:  [HYPHYMP](https://github.com/veg/hyphy/releases/tag/2.2)


Usage
----------

* Use a reference fasta where each contig starts at an open reading frame (ORF).
* Align paired-end reads against your reference with your favourite pairwise fast read aligner to produce a SAM  file of alignments.
* Clone the SlidingWindow github repository.
* Invoke:  python sliding_window_tree.py {Parameters}
    * **Parameters:**
    * --sam:  full filepath to sam alignment file
    * --ref:  name of reference contig
    * --ref_len:  length of reference contig in nucleotide bases
    * --out_dir:  output directory in which the pipeline will write all its intermediate files
    * --map_q_thresh:  mapping quality threshold below which alignments are ignored.
    * --read_q_thresh:  read quality threshold below which bases are converted to Ns.
    * --max_N_thresh:  maximum fraction of Ns allowed in the merged paired-end read below which the paired-end read is ignored.
    * --window_size:  window size in nucleotides
    * --window_breadth_thresh:  fraction of window that merged paired-end read must cover with non-gap and non-N nucleotides.  Below this threshold, the read is omitted from the window.
    * --window_depth_thresh:  total reads that must meet the window breadth threshold.  Below this depth threshold, the window is omitted.
    * --start_nucpos:  1-based start nucleotide position in the reference contig.  The first window will start at this position.
    * --end_nucpos:  1-based end nucleotide position in the reference contig.  The last window will start at or before this position.
    * --pvalue: p-value threshold for determining selection significantly different from neutral evolution.
    * --threads_per_window:  threads allotted per window.
    * --concurrent_windows:  max number of windows to process concurrently.
    * --dnds_tsv:  full filepath of final tab-separated values file containing selection information for each codon site in the reference.
    * --hyphy_exe:  full filepath of HYPHYMP executable
    * --hyphy_basedir:  full filepath of HyPhy base directory containing template batch files
    * --fastree_exe:  full filepath of FastTreeMP executable
* The pipeline will output dN/dS values for each codin site in the tab-separated file as specified in --dnds_tsv.
* At this point, there is no pretty visualization of the dN/dS, so you will need to feed this .tsv file into R or similar for charting.

Parallelization
-----------------
Currently, the pipeline code is capable of running on multiple cores on a single node.  If you would like to run the pipeline
on multiple nodes on a cluster, you will need to do some manual or scripted massaging to kick off the pipeline for
multiple sections of the genome concurrently.  For each node, specify the genomic region that you'd like to process
using the --start_nucpos and --end_nucpos parameters.


