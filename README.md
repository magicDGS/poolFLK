poolFLK
=======

FLK (Bonhomme et al. 2010) is a test for selection scans based on population allele frequencies, similar to the Lewontin and Krakauer's (LK) test but using the kinship matrix
estimated from the population tree. Because it only requires allele frequencies for its computation, it is suitable to use with Pool-Seq data.

To my knowledge, only one [tool](https://qgsp.jouy.inra.fr/index.php?option=com_content&view=article&id=50&Itemid=55) implemented in R can handle allele frequencies directly and compute FLK.
Neverteless, the rooting of the NJ tree is limited to the midpoint algorithm. As suggested in the [hapFLK]() implementation, this could introduce bias in the kinship matrix and thus in the test.
Due to this bias, the authors implemented in their 1.3.0 version a maximum likelihood algorithm for rooting the tree, but this software only works with individual genotype formats.

To be able to compute FLK with the improvements developed in hapFLK, I develop this python script to perform the FLK from Pool-Seq data stored in the sync format. I use the hapFLK libary to perform
the analyses. Thus, if new versions of hapFLK are installed the script could be use with the novelties include in the library.

## Requirements

* [hapFLK](https://forge-dga.jouy.inra.fr/projects/hapflk/files) library (version >= 1.3.0 for the ML algorithm)
* [heapq](https://docs.python.org/2/library/heapq.html) library
