Overview
============================================

Umberjack unit tests require simulated population genomes and simulated sequencing.

By default, the simulated data are generating using the  linux_x64 tool binaries included in /Umberjack/test/simulations/bin.
If you are executing code on a different platform or you want to use different versions of the tools supplied, 
you will need to edit /Umberjack/test/simulations/data/umberjack_unittest/umberjack_unittest_sim.conf


Simulation Required Tools
===============================================

* [FastTree](http://www.microbesonline.org/fasttree/):  FastTree performs phylogenetic reconstruction.  
  * FastTree executable is the single threaded version of FastTree.  This is required if you want to ensure that your simulated trees are generated deterministically.
  * FastTreeMP executable is the multithreaded version of FastTree.  It will not generate deterministic trees.
* [HYPHYMP](https://github.com/veg/hyphy/releases):  This is the OpenMP multithreaded version of HyPhy.  It calculates various evolutionary statistics.
* [INDELible](http://abacus.gene.ucl.ac.uk/software/indelible/)
* [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
* [BWA](http://bio-bwa.sourceforge.net/)
* [Picard](http://broadinstitute.github.io/picard/)
* [R](http://www.r-project.org/):  R is used for displaying alignment stats for simulated reads and plotting the genome wide evolutionary stats
  * [R ggplot2](http://ggplot2.org/) package.
  * [R knitr](http://yihui.name/knitr/) package
  * [R reshape2](http://cran.r-project.org/web/packages/reshape2/index.html) package
  * [R plyr](http://cran.r-project.org/web/packages/plyr/index.html) package


Simulation Required 3rd Party Libraries, Binaries
===============================================

* HyPhy requires [OpenMP](http://openmp.org/wp/openmp-compilers) libraries.
* Picard requires [Java](http://www.oracle.com/technetwork/java/javase/overview/index.html)
* Umberjack MPI tests require [mpi4py](http://mpi4py.scipy.org/docs/usrman/index.html) and [OpenMP](http://www.open-mpi.org/) implementation of [MPI-2](http://en.wikipedia.org/wiki/Message_Passing_Interface).


How to Generate Simulated Data
===========================================

We include a pipeline for simulating populations and simulating shotgun sequence for those population.
  * /Umberjack/blob/master/test/simulations/asg_driver.py generates an ancestral selection graph.
  * [INDELible](http://abacus.gene.ucl.ac.uk/software/indelible/) simulates genomes for multiple populations based on the ancestral selection graph.  Each population has a different mutation scaling rate. 
  * /Umberjack/test/simulations/sample_genomes.py creates a new population whose codon sites are selected randomly from an INDELible population genome.
  * [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) simulates shotgun MiSeq paired-end sequencing for the new population.
  * [BWA](http://bio-bwa.sourceforge.net/) aligns the shotgun reads against the new population's consensus.
  * [Picard](http://broadinstitute.github.io/picard/) collects alignment statistics
  * [FastTree](http://www.microbesonline.org/fasttree/) creates a phylogenetic tree for the full population.
  * [HYPHYMP](https://github.com/veg/hyphy/releases) collects the evolutionary statistics for the full population. 
 
In order to execute this pipeline:

* Create a simulated data config file.  Use /Umberjack/test/simulations/data/umberjack_unittest/umberjack_unittest_sim.conf as a template.  
The simulation pipeline will generate the simulated population and simulated reads under the same directory as your config file. 
* Add /Umberjack/test to your PYTHONPATH environment variable
* Execute 

```
python /Umberjack/test/simulations/sim_pipeline.py [path to config file]
```

