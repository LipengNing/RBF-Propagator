function [RTAP, RTPP]=ComputeRTAPandRTPP(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

[gg,~]=icosahedron(4); %gradient directions for ODF
L=length(gg);
rtap=zeros(nx,ny);
rtpp=zeros(nx,ny);

for ix=1:nx
    parfor iy=1:ny
        [RTAP_basis, RTPP_basis]=ConstructBasis(D(:,:,ix,iy),sigma);      
        rtap(ix,iy)=RTAP_basis*V(:,ix,iy);
        rtpp(ix,iy)=RTPP_basis*V(:,ix,iy);
    end
end
RTAP=zeros(size(mask));
RTAP(mask~=0)=rtap;
RTPP=zeros(size(mask));
RTPP(mask~=0)=rtpp;
%%
if(nargin>=4&& (~isempty(FileName)))
% sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
% mat2nhdr(RTAP, [FileName,'_RTAP'], 'custom',sd);
% mat2nhdr(RTPP, [FileName,'_RTPP'], 'custom',sd);
nii=make_nii(RTAP);
save_nii(nii,[FileName,'_RTAP.nii']);
clear nii;
nii=make_nii(RTPP);
save_nii(nii,[FileName,'_RTPP.nii']);
end



end


function [RTAP_basis, RTPP_basis]=ConstructBasis(D0,sigma)
t = 70*1e-3;
b_basis=[2 4]*1000; 
[ug,~]=icosahedron(2);
N_ug=length(ug);
ug=ug(1:N_ug/2,:);
N_ug=N_ug/2;
q_basis=(1/(2*pi))*sqrt(b_basis/t);
qg_basis=kron(diag(q_basis),eye(N_ug))*repmat(ug,[length(b_basis),1]);
N_basis=length(qg_basis);

N_column=N_basis+1;

RTAP_basis=zeros(1,N_column);
RTPP_basis=zeros(1,N_column);
[U,s,~]=svd(D0);

RTAP_basis(end)=pi*(s(2,2)*s(3,3))^(-.5);
RTPP_basis(end)=pi^(.5)*s(1,1)^(-.5);

RTAP_basis(1:N_basis)=2*pi*sigma(2)^(-1)*exp(-sigma(1)* (U(:,1)'*qg_basis').^2);
RTPP_basis(1:N_basis)=2*pi^(.5)*sigma(1)^(-.5)*exp(-sigma(2)* (U(:,2)'*qg_basis').^2-sigma(2)* (U(:,3)'*qg_basis').^2 );

end