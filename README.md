# RS_DBSCAN
RS-DBSCAN is based on two key assumptions: (1) The data of each sufficient sample should have similar data distribution, as well as category distribution, to the entire data set; (2) the representative of each category in all sufficient samples conform to Gaussian distribution. It processes data in two stages, one is to classify data in each local sample independently, and the other is to globally classify data by assigning each point to the category of its nearest representative category center.

***********************************************************************************
The RS-DBSCAN program was compiled under Windows using matlab R2016b.
***********************************************************************************

Files
===================================================================================
These program mainly containing:

-startup code named `Huge_speed.m`.

-one main functions of C4Y named `calculateClusterLabels_Large.m`.

-The main code for running the RNN-DBSCAN algorithm named `RnnDbscan.m` and `RNN-DBSCAN tests.ipynb`, main functions of RNN-DBSCAN named `knnIndexToGraphEdges.m`, `knngragh.m` and `knnindex.m`.

-some data sets.

