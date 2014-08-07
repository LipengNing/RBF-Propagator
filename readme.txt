This package is used to estimate the diffusion propagator using radial basis functions. Other functions for NRRD/Nifti file reading and writing are also needed to run this package. Some of these functions can be found via 
https://sites.google.com/site/kittipat/mvpa-for-brain-fmri/convert_matlab_nifti



Demo_CollectIndices.m	       the main file that illustrates how to use this package to collect scalar indices.

CollectIndices.m 	       the file that computes scalar indices

run_GaussianBasis.m	       computes the representation coefficients and D0 in DTI
run_GaussianBasis_Decay.m	similar to run_GaussianBasis.m but with stronger constraints on the estimated signal

GaussianBasis_Coef.m		is used to compute the representation coefficient in the radial basis functions and diffusion tensor.
GaussianBasis_Coef_Cons.m	similar to GaussianBasis_Coef.m but with stronger constraints


ComputeRTOP.m			computes the estimated return-to-the-origin probability (RTOP).
ComputeRTAP.m			computes the estimated return-to-the-axis probability (RTAP).
ComputeRTPP.m			computes the estimated return-to-the-plane probability (RTPP).
ComputeMSD.m			computes the estimated mean-squared-displacement (MSD).
ComputeODF.m			computes the estimated diffusion orientation-distribution-function (ODF).
ComputeRD.m			computes the minimum-mean-squared-displacement-error (MMSDE) and the relative-diffusivity (RD).
ComputeKurtosis.m		computes the generalized-kurtosis (GK) and the non-isotropic-Gaussianity (NIG).
.
.
.


TestGaussianBasis.m		is used to check if there is bug in the code. If it returns 1, then the code is correct.
TestData.mat			includes a set of measurements that is used in TestGaussianBasis.m to verify the code.



