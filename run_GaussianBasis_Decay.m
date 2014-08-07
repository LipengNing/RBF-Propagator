function [D,V,mask] = run_GaussianBasis_Decay(input_dwi,mask_name,sigma)
if(nargin<3)
    sigma=[0.0015 0.0008];
end

[SS, u0, b0, ~, ~, ss0] = nhdr_diff_multiB(input_dwi);

if (~isempty(strfind(mask_name,'.nii'))) %.nii format
    mask=MRIread(mask_name);
    mask=mask.vol;
    %mask=permute(mask,[2 1 3]);%the coordinates may be permuted
else
    mask = nrrdZipLoad(mask_name);%.nrrd format
end

mask(ss0<0)=0;


N=size(SS,4);
b=b0(1:N);
u=u0(1:N,:);

%%
[nx,ny,nz,N]=size(SS);
L  = nx*ny*nz;
SS = reshape(SS,1,1,L,N);
SS = SS(:,:,mask~=0,:);
Decay=[0:1000:18000]; %b values for constraints of signal
[D,V] = GaussianBasis_Coef_Cons(SS,u,b,Decay,sigma);

end

