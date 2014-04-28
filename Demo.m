clear all
close all

RunTest=1;

if(~matlabpool('size'))
    matlabpool open;
end


if(RunTest)
    if(~TestGaussianBasis())
        error('Code is not correct!');
    end
end
    
    

% input_dwi = '/projects/schiz/pi/yogesh/phd/dwmri/Data/Multishell/data_3x_set2_working.nhdr';
% mask_name = '/projects/schiz/pi/yogesh/phd/dwmri/Data/Multishell/mask.nhdr';

% input_dwi = '/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/case002/SparseDIffusion_InVivo_02072013/SMS2_dirs60/Ed-SMS2_60.nhdr';
% mask_name = '/projects/schiz/pi/lipeng/matlab/Code/Case3/mask.nrrd';

input_dwi = '/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/case003/dwi_3shells_b3000.nhdr';
mask_name = '/projects/schiz/pi/lipeng/matlab/Code/Case3/mask.nrrd';



[D,V,mask] = run_GaussianBasis(input_dwi,mask_name);%% change u and b inside the function

FileName='Case3_Coef';

save(FileName);

%%

%%compute return-to-the-origin probability
FileName='Case3_RTOP';
ComputeRTOP(D,V,mask,FileName);

%%compute mean-squared-displacement 
FileName='Case3_MSD';
ComputeMSD(D,V,mask,FileName);


%%compute restricted diffusion 
FileName='Case3_RD';
ComputeRD(D,V,mask,FileName);

%%compute Kurtosis: GK and NIG
FileName='Case3_Kurtosis';
ComputeKurtosis(D,V,mask,FileName);


%%compute ODF
%[ODF,gg]=ComputeODF(D,V,mask,sigma)

