function MSD=ComputeMSD_Qspace(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

msd=zeros(nx,ny);
for ix=1:nx
    parfor iy=1:ny
        MSD_basis=ConstructM2Basis_Qspace(D(:,:,ix,iy),sigma);
        Vec=MSD_basis*V(:,ix,iy);
        Vec(Vec<0)=0;
        msd(ix,iy)=sum(Vec);
    end
end
MSD=zeros(size(mask));
MSD(mask~=0)=msd/max(msd(:))*1000;

%%
if(nargin>=4&& (~isempty(FileName)))
% sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
% mat2nhdr(MSD, FileName, 'custom',sd);
nii=make_nii(MSD);
save_nii(nii,[FileName,'_QMSD.nii']);
end




end


function MSD_basis=ConstructM2Basis_Qspace(D0,sigma)
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
[U,~,~]=svd(D0);

MSD_basis=zeros(3,N_column);%coefficients in the covariance matrix
R=inv(D0)/2;
const=(2*pi)^(1.5)*sqrt(det(R));
MSD_basis(:,end)=const*[R(1,1) R(2,2) R(3,3)]';

D=U*diag([sigma(1) sigma(2)*[1 1]])*U';

R=inv(D)/(2);

mu=qg_basis';
const=2*(2*pi)^(1.5)*sqrt(det(R))*ones(1,N_basis);

M=repmat(R,[1,N_basis])+mu*(kron(eye(N_basis),ones(1,3)).*repmat(mu',[1 N_basis]));

M=reshape(M,3,3,N_basis);
M=reshape(M,9,N_basis);
MSD_basis(:,1:N_basis)=repmat(const,[3 1]).*M([1 5 9],:);

end
