Overview
============================================

Umberjack unit tests require simulated population genomes and simulated sequencing.

By default, the simulated data are generating using the  linux_x64 tool binaries included in /Umberjack/test/simulations/bin.
If you are executing code on a different platform or you want to use different version of the tools supplied, 
you will need to edit /Umberjack/test/simulations/data/umberjack_unittest/umberjack_unittest.conf.


Simulation Required Tools
===============================================

[FastTree](http://www.microbesonline.org/fasttree/)
[FastTreeMP](http://www.microbesonline.org/fasttree/)
[HyPhy](https://github.com/veg/hyphy/releases)
[BWA](http://bio-bwa.sourceforge.net/)
[INDELible](http://abacus.gene.ucl.ac.uk/software/indelible/)
[Picard](http://broadinstitute.github.io/picard/)



Simulation Required 3rd Party Libraries, Binaries
===============================================

HyPhy requires OpenMP.
Picard requires [Java](http://www.oracle.com/technetwork/java/javase/overview/index.html)
Umberjack MPI tests require MPI




