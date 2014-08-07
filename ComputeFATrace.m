input_dwi='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate001/diff/dwi_2.nhdr';
mask_name='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate001/diff/dwi_mask.nii.gz';
ExportName='/projects/schiz/pi/lipeng/matlab/CaseStudy/Pediatric/Cardiac001/Indicies';

[FA,D_bar] = run_FA_trace(input_dwi,mask_name);

nii=make_nii(FA);
save_nii(nii,[ExportName,'_FA.nii']);

nii=make_nii(D_bar);
save_nii(nii,[ExportName,'_Trace.nii']);
1


clear all


input_dwi='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate002/diff/dwi-Ed.nhdr';
mask_name='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate002/diff/dwi-Ed_mask.nii.gz';
ExportName='/projects/schiz/pi/lipeng/matlab/CaseStudy/Pediatric/Cardiac002/Indicies';
[FA,D_bar] = run_FA_trace(input_dwi,mask_name);

nii=make_nii(FA);
save_nii(nii,[ExportName,'_FA.nii']);

nii=make_nii(D_bar);
save_nii(nii,[ExportName,'_Trace.nii']);

2
clear all


input_dwi='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate003/diff/dwi-Ed.nhdr';
mask_name='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate003/diff/dwi-Ed_mask.nii.gz';
ExportName='/projects/schiz/pi/lipeng/matlab/CaseStudy/Pediatric/Cardiac003/Indicies';
[FA,D_bar] = run_FA_trace(input_dwi,mask_name);

nii=make_nii(FA);
save_nii(nii,[ExportName,'_FA.nii']);

nii=make_nii(D_bar);
save_nii(nii,[ExportName,'_Trace.nii']);

3


input_dwi='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate005/diff/dwi-Ed.nhdr';
mask_name='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate005/diff/dwi-Ed_mask.nii.gz';
ExportName='/projects/schiz/pi/lipeng/matlab/CaseStudy/Pediatric/Cardiac005/Indicies';

[FA,D_bar] = run_FA_trace(input_dwi,mask_name);

nii=make_nii(FA);
save_nii(nii,[ExportName,'_FA.nii']);

nii=make_nii(D_bar);
save_nii(nii,[ExportName,'_Trace.nii']);

5
clear all

input_dwi='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate004/diff/dwi-Ed.nhdr';
mask_name='/projects/schiz/pi/yogesh/phd/dwmri/Data/HiResMultishell/Pediatric/CardiacNeonate004/diff/dwi-Ed_mask.nii.gz';
ExportName='/projects/schiz/pi/lipeng/matlab/CaseStudy/Pediatric/Cardiac004/Indicies';

[FA,D_bar] = run_FA_trace(input_dwi,mask_name);

nii=make_nii(FA);
save_nii(nii,[ExportName,'_FA.nii']);

nii=make_nii(D_bar);
save_nii(nii,[ExportName,'_Trace.nii']);


4
