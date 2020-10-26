# mtmodels

This repository the contains the implementations of the agent-based model (AB) and the spatially-explicit model (SE) utilized in 'Universal dynamics of mitochondrial networks - A finite size scaling analysis' to simulate mitochondrial networks dynamics. Additional code and images to reproduce results from finite size saling analysis of real mitochondrial networks are also provided.

# AB model
AB.f90: Fortran code to implement the AB model. It requires random_module.90 (also provided). 

# SE model
SE.R: R code to implement the SE model.

# Real networks
CoMNStA.m (Complex Mitochondrial Network Analyzer): MATLAB script to reproduce Figure 7 from high-res microscopy images of mitochondrial networks from mouse embryonic fibroblasts (MEFs) (contained in the folder 'images').
