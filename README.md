# ISmethod
Matlab script implementing the iterative simulation method for recovering the diffusive coefficient from diffusive trajectories with non-harmonic confinement

This script is implementing the iterated simulations (IS) method for recovering the diffusive coefficient from a Brownian confined trajectory (original dataset) removing errors associated to non-harmonicity of the confining potential. For this, simulations are made with different diffusion coefficients that the user has to define (Dstore parameter) within a reasonable range of values.
For each element of Dstore, a set of brownian dynamics simulations is made. The  output is a plot of the fitted D from both the original dataset and the set of simulations. The output diffusion coefficient has to be identified as the input element of Dstore so that the dataset-simulation discrepancy in the fitted diffusion coefficients are minimize. An interpolation can be made to identify the ottimal D. Looking at the diffusion plot, the proper diffusion coefficient is te value of D_{in} where "D_{fit} original dataset" and "D_{fit} IS method" intersect

ORIGINAL DATASET
In the current version, the original dataset is simulated using an asymmetric potential as defined in the reference paper. 
For using an experimental dataset, replace the 'generate original dataset' section with the upload of the experimental dataset. 

For questions/comments please write to stefano.villa@ds.mpg.de
