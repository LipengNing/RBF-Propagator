clear 
close all
%%
RunTest=1;
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
mask_name='MaskFolder/MaskFile.nrrd'; %%.nii format can also be used, but make sure the coordinates are not permuted
ExportName='ExportFolder/Prefix';
flag=0; %flag=0 regular constraints, flag=1 for stronger constraints; 
CollectIndices(input_dwi,mask_name,ExportName);