function RTOP=ComputeRTOP(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

[gg,~]=icosahedron(4); %gradient directions for ODF
L=length(gg);
rtop=zeros(nx,ny);
for ix=1:nx
    parfor iy=1:ny
        RTOP_basis=ConstructRTOPBasis(D(:,:,ix,iy),sigma);      
        rtop(ix,iy)=RTOP_basis*V(:,ix,iy);
    end
end
RTOP=zeros(size(mask));
RTOP(mask~=0)=rtop;
%%
if(nargin>=4&& (~isempty(FileName)))
% sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
% mat2nhdr(RTOP, FileName, 'custom',sd);
nii=make_nii(RTOP);
save_nii(nii,[FileName,'_RTOP.nii']);
end

end


function RTOP_basis=ConstructRTOPBasis(D0,sigma)
t = 70*1e-3;
b_basis=[2 4]*1000; % location of basis functions
[ug,~]=icosahedron(2);
N_ug=length(ug);
ug=ug(1:N_ug/2,:);
N_ug=N_ug/2;
q_basis=(1/(2*pi))*sqrt(b_basis/t);
qg_basis=kron(diag(q_basis),eye(N_ug))*repmat(ug,[length(b_basis),1]);%points used to guarantee decay
N_basis=length(qg_basis);

N_column=N_basis+1;

RTOP_basis=zeros(1,N_column);%coefficients in the covariance matrix

RTOP_basis(end)=pi^(1.5)/sqrt(det(D0));

RTOP_basis(1:N_basis)=2*pi^(1.5)/sqrt(sigma(1)*sigma(2)^2)*ones(1,N_basis);

end
