function [FA,D_bar] = run_FA_trace(input_dwi,mask_name,IND)

[SS, u0, b0, ~, ~, ~] = nhdr_diff_multiB(input_dwi);
[nx,ny,nz,N]=size(SS);
if(~isempty(mask_name))
mask = nrrdZipLoad(mask_name);
mask(ss0<0)=0;
else
    mask=ones(nx,ny,nz);
end




SS=SS(:,:,:,IND);
u=u0(IND,:);
b=b0(IND);
N=length(IND);
SS=reshape(SS,nx*ny*nz,N);
SS=SS(mask~=0,:);
[d,D]=direct_1T(u,b,SS');
fa=zeros(1,size(SS,1));
d_bar=zeros(1,size(SS,1));

parfor i=1:size(SS,1)
    [fa_i,d_i]=tensor2fa(D(:,:,i));
    fa(i)=fa_i;
    d_bar(i)=d_i;
end
FA=zeros(nx,ny,nz);
FA(mask~=0)=fa;
D_bar=zeros(nx,ny,nz);
D_bar(mask~=0)=d_bar;
end


