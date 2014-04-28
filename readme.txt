This package is used to estimate the diffusion propagator using radial basis functions.

Demo.m				is the main files that illustrates the usage of all the routines.

run_GaussianBasis.m		is used to load the data and call the function GaussianBasis_Coef.
GaussianBasis_Coef.m	is used to compute the representation coefficient in the radial basis functions and diffusion tensor.

ComputeRTOP.m		computes the estimated return-to-the-origin probability (RTOP).
ComputeMSD.m			computes the estimated mean-squared-displacement (MSD).
ComputeODF.m			computes the estimated diffusion orientation-distribution-function (ODF).
ComputeRD.m			computes the minimum-mean-squared-displacement-error (MMSDE) and the relative-diffusivity (RD).
ComputeKurtosis.m		computes the generalized-kurtosis (GK) and the non-isotropic-Gaussianity (NIG).

TestGaussianBasis.m		is used to check if there is bug in the code. If it returns 1, then the code is correct.
TestData.mat			includes a set of measurements that is used in TestGaussianBasis.m to verify the code.



