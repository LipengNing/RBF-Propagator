function MSD=ComputeMSD(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

msd=zeros(nx,ny);
for ix=1:nx
    parfor iy=1:ny
        MSD_basis=ConstructMSDBasis(D(:,:,ix,iy),sigma);
        Vec=MSD_basis*V(:,ix,iy);
        msd(ix,iy)=sum(Vec);
    end
end
MSD=zeros(size(mask));
MSD(mask~=0)=msd;

%%
if(nargin>=4&& (~isempty(FileName)))
% sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
% mat2nhdr(MSD, FileName, 'custom',sd);
nii=make_nii(MSD);
save_nii(nii,[FileName,'_MSD.nii']);
end

end


function MSD_basis=ConstructMSDBasis(D0,sigma)
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
MSD_basis=zeros(3,N_column);%coefficients in the covariance matrix


R=D0/(2*pi^2);
MSD_basis(:,end)=[R(1,1) R(2,2) R(3,3)]';

[U,~,~]=svd(D0);
D=U*diag([sigma(1) sigma(2)*[1 1]])*U';


% R=D/(2*pi^2);
% mu=2*pi*R*qg_basis';
% const=2*exp(-2*pi^2*sum((qg_basis*R).*qg_basis,2));
% 
% M=repmat(R,[1,N_basis])-mu*(kron(eye(N_basis),ones(1,3)).*repmat(mu',[1 N_basis]));
% 
% M=reshape(M,3,3,N_basis);
% M=reshape(M,9,N_basis);
% MSD_basis(:,1:N_basis)=repmat(const',[3 1]).*M([1 5 9],:);


for n=1:N_basis
        R=D/(2*pi^2);       
        mu=2*pi*R*qg_basis(n,:)';
        const=2*exp(-2*pi^2*qg_basis(n,:)*R*qg_basis(n,:)');
        M=R-(mu*mu');
        MSD_basis(:,n)=const*[M(1,1) M(2,2) M(3,3)]';
end



end
