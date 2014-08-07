function signal=ComputeSignal(D,V,qg,mask,sigma)
% interportate the diffusion signal at location qg in q-space
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);
N=length(qg);
signal=zeros(nx,ny,N);

for ix=1:nx
    parfor iy=1:ny
        Signal_basis=ConstructSignalBasis(D(:,:,ix,iy),qg,sigma);
        signal(ix,iy,:)=Signal_basis*V(:,ix,iy);
    end
end
signal=reshape(signal,nx*ny,N);

[nx,ny,nz]=size(mask);
S=zeros(nx*ny*nz,N);
S(mask~=0,:)=signal;
S=reshape(S,nx,ny,nz,N);

end


function Signal_basis=ConstructSignalBasis(D0,qg,sigma)
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
N=length(qg);
Signal_basis=zeros(N,N_column);

x_gold=sum( (qg*D0).*qg,2);
Signal_basis(:,end)=exp(-x_gold);

D=U*diag([sigma(1) sigma(2)*[1 1]])*U';

qg=qg*U;
qg_basis=qg_basis*U;
X1=((repmat(qg,[1 N_basis])-repmat(vec(qg_basis')',[N,1])).^2)*kron(eye(N_basis),[sigma(1);sigma(2);sigma(2)]);
X2=((repmat(qg,[1 N_basis])+repmat(vec(qg_basis')',[N,1])).^2)*kron(eye(N_basis),[sigma(1);sigma(2);sigma(2)]);
Signal_basis(:,1:N_basis)=exp(-X1)+exp(-X2);


end
