clear 
close all
%%
RunTest=1;%Test error in the code using TestData.mat
if(RunTest)
    if(~TestGaussianBasis())
        error('Code is not correct!');
    end
end
%%
if(~matlabpool('size'))
    matlabpool open;
end
%%
input_dwi='DataFolder/InputFile.nhdr';
%A Nifit file for mask is needed. 
%The data structure for the mask and the dwi file should be the same, i.e.
%the x-y-z should be the same. If not, the coordinates need to be permuted in
% run_GaussianBasis_Decay.m and run_GaussianBasis.m
mask_name='MaskFolder/MaskFile.nii'; 

%Prefix for output file names
OutputPrefix='ExportFolder/Prefix';


flag=1; %flag=0 regular constraints, flag=1 for stronger constraints; 
%Setting flag=1 usually leads to more stable results but it will take much
%longer time to compute the parameters.
CollectIndices(input_dwi,mask_name,OutputPrefix,flag);