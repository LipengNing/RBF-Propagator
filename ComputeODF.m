function [ODF,gg]=ComputeODF(D,V,mask,sigma)
if(nargin<4)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

[gg,~]=icosahedron(4); %gradient directions for ODF
L=length(gg);
odf=zeros(nx,ny,L);
for ix=1:nx
    parfor iy=1:ny
        ODF_basis=ConstructODFBasis(D(:,:,ix,iy),gg,sigma);
        odf(ix,iy,:)=ODF_basis*V(:,ix,iy);
    end
end
odf=reshape(odf,nx*ny,L);

[nx,ny,nz]=size(mask);
ODF=zeros(nx*ny*nz,L);
ODF(mask~=0,:)=odf;
ODF=reshape(ODF,nx,ny,nz,L);

end


function ODF_basis=ConstructODFBasis(D0,gg,sigma)
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

ODF_basis=zeros(size(gg,1),N_column);



x_u=sum((gg/D0).*gg,2);
ODF_basis(:,end)=1./(4*pi*sqrt(det(D0))*x_u.^(1.5));

D=U*diag([sigma(1) sigma(2)*[1 1]])*U';

x=sum((gg/D).*gg,2);
temp=.5*repmat(x.^(-1.5),[1 N_basis])-repmat(x.^(-2.5),[1 N_basis]).*((gg*qg_basis').^2);
ODF_basis(:,1:N_basis)=temp.*exp(-((gg*qg_basis').^2./(repmat(x,[1 N_basis]))))/(pi*sqrt(det(D)));

% for n=1:N_basis
%         x=sum((gg/D).*gg,2);
%         temp=(.5*x.^(-1.5)-x.^(-2.5).*((gg*qg_basis(n,:)').^2));
%         ODF_basis(:,n)=temp.*exp(-((gg*qg_basis(n,:)').^2./x))/(pi*sqrt(det(D)));
% end


end
