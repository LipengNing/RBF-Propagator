function [D,V,mask] = run_GaussianBasis(input_dwi,mask_name,sigma)
if(nargin<3)
    sigma=[0.0015 0.0008];
end

[SS, u0, b0, voxel, Sm, ss0] = nhdr_diff_multiB(input_dwi);
mask = nrrdZipLoad(mask_name);
mask(ss0<0)=0;
id = find(mask ~= 0);

N=size(SS,4);
b=b0(1:N);
u=u0(1:N,:);

%%
[nx,ny,nz,N]=size(SS);
L  = nx*ny*nz;
SS = reshape(SS,1,1,L,N);
SS = SS(:,:,mask~=0,:);

[D,V] = GaussianBasis_Coef(SS,u,b);

end

