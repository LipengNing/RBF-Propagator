function [D_store,V] = GaussianBasis_Coef(S,u,b,sigma)
% S : signal of size (nx,ny,nz,n)
% u : sample gradient directions
% b : b values
%%%
%%%%%%%% parameters  %%%%%%%%%%%%%%%%%%%
if(nargin<4)
    sigma=[0.0015 0.0008];
end

% CondNumber=1e7;
%%%%%%%%%%%%%%%%%%%%%%%%%
%%
t = 70*1e-3;
b=b(:);
q = (1/(2*pi))*sqrt(b/t);
qg=repmat(q,[1 3]).*u;

if(length(size(S))==3)
    [nx,ny,N]=size(S);
else 
    [nx,ny,nz,N]=size(S);
    S=reshape(S,[nx,ny*nz,N]);
    [nx,ny,N]=size(S);
end

[gg,~]=icosahedron(4); %gradient directions for ODF
L=length(gg);

    
D_store=zeros(3,3,nx,ny);
V=zeros(163,nx,ny);

for ix=1:nx
    parfor iy=1:ny
        E=squeeze(S(ix,iy,:)); % change (ix,iy,1) to (ix,1,iy) for MGHdata
        E=double(E);         
        e=sort(E,'descend'); %
        e0=mean(e(1:5));
        if(e0>1)
            E=E/e0; %
        end
        E(E>1)=1;
        D0=EstTensor(qg,E);
        D_store(:,:,ix,iy)=D0;
        [A,B,c]=ConstructBasisMatrix(D0,qg,sigma);
        lamb=0.00055;
        opts = optimset('Algorithm','active-set','Display','off','MaxIter',1000);
        v=quadprog(A'*A+lamb*eye(size(A,2)),-A'*E,-B,zeros(size(B,1),1),c,1,[],[],flipud(eye(size(A,2),1)),opts);    
        V(:,ix,iy)=v; 
    end
end

end



function [A,B,c]=ConstructBasisMatrix(D0,qg,sigma)

b_fine=[0:1:8]*1000; %
b_basis=[2 4]*1000; % location of basis functions

t = 70*1e-3;
N=size(qg,1);

[ug,~]=icosahedron(2);
N_ug=length(ug);% choose high resolution gradient directions used for constraints, this 
ug=ug(1:N_ug/2,:);
N_ug=N_ug/2;

q_basis=(1/(2*pi))*sqrt(b_basis/t);
qg_basis=kron(diag(q_basis),eye(N_ug))*repmat(ug,[length(b_basis),1]);%points used to guarantee decay
N_basis=length(qg_basis);

q_fine = (1/(2*pi))*sqrt(b_fine/t);
N_bfine=length(b_fine);
qg_fine=kron(diag(q_fine),eye(N_ug))*repmat(ug,[N_bfine,1]);%points used to guarantee decay
N_fine=length(qg_fine);

%%
N_column=N_basis+1;
[U,~,~]=svd(D0);
D=U*diag([sigma(1) sigma(2)*[1 1]])*U';

A=zeros(N,N_column);
B=zeros(N_fine,N_column);
c=zeros(1,N_column);

%%the Gaussian centered at zero
x_basis=sum((qg*D0).*qg,2);
A(:,end)=exp(-x_basis);
x_fine=sum((qg_fine*D0).*qg_fine,2);
B(:,end)=exp(-x_fine);   
c(end)=1;
%%other Gaussians
qg=qg*U;
qg_fine=qg_fine*U;
qg_basis=qg_basis*U;
X1=((repmat(qg,[1 N_basis])-repmat(vec(qg_basis')',[N,1])).^2)*kron(eye(N_basis),[sigma(1);sigma(2);sigma(2)]);
X2=((repmat(qg,[1 N_basis])+repmat(vec(qg_basis')',[N,1])).^2)*kron(eye(N_basis),[sigma(1);sigma(2);sigma(2)]);
A(:,1:N_basis)=exp(-X1)+exp(-X2);


X1=((repmat(qg_fine,[1 N_basis])-repmat(vec(qg_basis')',[N_fine,1])).^2)*kron(eye(N_basis),[sigma(1);sigma(2);sigma(2)]);
X2=((repmat(qg_fine,[1 N_basis])+repmat(vec(qg_basis')',[N_fine,1])).^2)*kron(eye(N_basis),[sigma(1);sigma(2);sigma(2)]);
B(:,1:N_basis)=exp(-X1)+exp(-X2);

c(1:N_basis)=2*exp(-(qg_basis.^2)*[sigma(1);sigma(2);sigma(2)]);

%%
% for n=1:N_basis
%     x1=sum(((qg-ones(N,1)*qg_basis(n,:))*D).*(qg-ones(N,1)*qg_basis(n,:)),2);
%     x2=sum(((qg+ones(N,1)*qg_basis(n,:))*D).*(qg+ones(N,1)*qg_basis(n,:)),2);
%     A(:,n)=exp(-x1)+exp(-x2);
% 
%     x1=sum(((qg_fine-ones(N_fine,1)*qg_basis(n,:))*D).*(qg_fine-ones(N_fine,1)*qg_basis(n,:)),2);
%     x2=sum(((qg_fine+ones(N_fine,1)*qg_basis(n,:))*D).*(qg_fine+ones(N_fine,1)*qg_basis(n,:)),2);
%     B(:,n)=exp(-x1)+exp(-x2);
% 
%     c(n)=2*exp(-(qg_basis(n,:)*D*qg_basis(n,:)'));
% end

Dif=kron(toeplitz(eye(1,N_bfine-1),[1 -1 zeros(1,N_bfine-2)]),eye(N_ug));%ensure signal decay along each direction
B=[B; Dif*B]; 

end


function D=EstTensor(q,E)
n=size(q,1);
E=E(:);
E(E>1)=1;
E=E(E>1e-3);
q=q(E>1e-3,:);
s=-log(E);
A=[q(:,1).^2  q(:,2).^2  q(:,3).^2  2*q(:,1).*q(:,2)  2*q(:,1).*q(:,3)  2*q(:,2).*q(:,3)];
d=pinv(A'*A)*A'*s;d=real(d);
D=[d(1) d(4) d(5);
   d(4) d(2) d(6);
   d(5) d(6) d(3)];
[U,e]=eig(D);
e=diag(e);
e(e<5e-4)=5e-4;
D=U*diag(e)*U';
end
