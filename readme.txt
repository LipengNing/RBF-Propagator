This package is used to estimate the diffusion propagator using Gaussian basis functions. It includes routines to compute the return-to-the-origin probability (RTOP), mean-squared-displacement (MSD), the scalar indices for the restricted diffusivity and the generalized Kurtosis. 

Demo.m 				is the main file that illustrate the usage these functions.

run_GaussianBasis.m	is used to load the data and to call the routine GaussianBasis_Coef.

GaussianBasis_Coef.m 	is the main routine that compute the representation coefficients for Gaussian basis functions.

ComputeRTOP.m 		computes the return-to-the-origion probability (RTOP).
ComputeMSD.m		computes the mean-squared-displacement (MSD).
ComputeKurtosis.m		computes the generalized kurtosis (GK) and the non-isotropic-Gaussianity (NIG).
CompueRD.m			computes the minimum-mean-squared-displacement-error (MMSDE) and relative diffusivity (RD).
ComputeODF.m		computes the diffusion orientation-distribution-function (ODF).

TestGaussianBasis.m 	is used to check if there is bug within the code. If it returns 1, then the code is correct.




