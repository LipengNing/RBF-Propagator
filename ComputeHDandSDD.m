function [HD,SDD]=ComputeHDandSDD(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

[gg,~]=icosahedron(4); %gradient directions for ODF
L=length(gg);
rd=zeros(nx,ny);
sdd=zeros(nx,ny);

for ix=1:nx
    parfor iy=1:ny
        M2_basis=ConstructM2Basis(D(:,:,ix,iy),sigma);
        m2=M2_basis*V(:,ix,iy);
        R=m2([1  4  5;
               4  2  6; 
               5  6  3]);
        [v,d]=eig(R);d=diag(d);d(d<2e-5)=2e-5;R=v*diag(d)*v';
        R0=D(:,:,ix,iy)/(2*pi^2);
        rd(ix,iy)=real(trace(R0+R-2*sqrtm(sqrtm(R0)*R*sqrtm(R0))));
        sdd(ix,iy)=real(trace(R0\R+R\R0))/2;
    end
end

HD=zeros(size(mask));
HD(mask~=0)=rd;

SDD=zeros(size(mask));
SDD(mask~=0)=sdd;


if(nargin>=4&& (~isempty(FileName)))
% sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
% mat2nhdr(HD, strcat(FileName,'_RD'), 'custom',sd);
% mat2nhdr(SDD, strcat(FileName,'_SDD'), 'custom',sd);
nii=make_nii(HD);
save_nii(nii,[FileName,'_HD.nii']);
clear nii;
nii=make_nii(SDD);
save_nii(nii,[FileName,'_SDD.nii']);
end

end


function MSD_basis=ConstructM2Basis(D0,sigma)
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

MSD_basis=zeros(6,N_column);%coefficients in the covariance matrix
R=D0/(2*pi^2);
MSD_basis(:,end)=[R(1,1) R(2,2) R(3,3) R(1,2) R(1,3) R(2,3)]';

D=U*diag([sigma(1) sigma(2)*[1 1]])*U';

R=D/(2*pi^2);
mu=2*pi*R*qg_basis';
const=2*exp(-2*pi^2*sum((qg_basis*R).*qg_basis,2));

M=repmat(R,[1,N_basis])-mu*(kron(eye(N_basis),ones(1,3)).*repmat(mu',[1 N_basis]));

M=reshape(M,3,3,N_basis);
M=reshape(M,9,N_basis);
MSD_basis(:,1:N_basis)=repmat(const',[6 1]).*M([1 5 9 4 7 8],:);


% for n=1:N_basis
%         R=D/(2*pi^2);       
%         mu=2*pi*R*qg_basis(n,:)';
%         const=2*exp(-2*pi^2*qg_basis(n,:)*R*qg_basis(n,:)');
%         M=R-(mu*mu');
%         MSD_basis(:,n)=const*[M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3)]';
% end

end
