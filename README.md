# cylic
cyclic is a class whose instances act as substitutions for Matlab's grid tensor. 
For multi-dimensional grids, memory requirements can grow quickly, 
while the variation exists along only one of the dimensions. 
This class implements correct indexing, while storing this minimal amount of information, along with
the desctiptive metadata, namely, tensor's size and dimension of variation.
NOTE: this project is under construction