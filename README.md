# fMRI-PC

Before running any of the code first download and install QUIC (a copy of the original repository is here: https://github.com/linboqiao/SICS/tree/master/QUIC and this is the original paper: http://jmlr.org/papers/volume15/hsieh14a/hsieh14a.pdf) and also PMTK3-Graphical Models (https://github.com/probml/pmtk3/tree/master/toolbox/GraphicalModels).

This code has three parts, the first two are used to run a cross-validation scheme to both examine the value of applying the adaptive graphical lasso while also, if there is value to it, estimating the optimal penalization value to apply to the structured prior (the structural connectome here) being applied through the AGL. 
