function [D_store,V] = GaussianBasis_Coef_Cons(S,u,b,Decay,sigma)
% S : signal of size (nx,ny,nz,n)
% u : sample gradient directions
% b : b values
%%%
%%%%%%%% parameters  %%%%%%%%%%%%%%%%%%%
if(nargin<5)
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
        ReRun=1;
        while(ReRun)  
        e=sort(E,'descend'); %
        e0=mean(e(1:5));
        if(e0>1)
            E=E/e0; %
        end

        D0=EstTensor(qg,E);
        
        D_store(:,:,ix,iy)=D0;
        
        [A,B,c]=ConstructBasisMatrix(D0,qg,Decay,sigma);
%         lambA=eig(A'*A);
%         lamb=max(lambA)/CondNumber-min(lambA);
%         lamb=max(lamb,0)
        lamb=0.00055;
        opts = optimset('Algorithm','active-set','Display','off','MaxIter',1000);
        v=quadprog(A'*A+lamb*eye(size(A,2)),-A'*E,-B,zeros(size(B,1),1),c,1,[],[],[],opts);    
        V(:,ix,iy)=v; 
        
        [~,M4_basis]=ConstructM2M4(D0,sigma);
        m4=M4_basis*V(:,ix,iy);
        M4=m4([ 1   4   5   4   10  13  5   13  11;
                4   10  13  10  6   14  13  14  15;
                5   13  11  13  14  15  11  15  8;
                4   10  13  10  6   14  13  14  15;
                10  6   14  6   2   7   14  7   12;
                13  14  15  14  7   12  15  12  9;
                5   13  11  13  14  15  11  15  8;
                13  14  15  14  7   12  15  12  9;
                11  15  8   15  12  9   8   9   3]);
          if(min(eig(M4))<-1e-11)
              E=E*0.95;
          else
              ReRun=0;
          end
        end
    end
end

end



function [A,B,c]=ConstructBasisMatrix(D0,qg,Decay,sigma)

b_fine=Decay; %
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



function [M2_basis,M4_basis]=ConstructM2M4(D0,sigma)
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

M2_basis=zeros(6,N_column);%coefficients in the covariance matrix
M4_basis=zeros(15,N_column);

R=D0/(2*pi^2);
M2_basis(:,end)=[R(1,1) R(2,2) R(3,3) R(1,2) R(1,3) R(2,3)]';

M4_basis(:,end)=[3*R(1,1)^2;            
           3*R(2,2)^2;           
           3*R(3,3)^2; 
           3*R(1,1)*R(1,2);     
           3*R(1,1)*R(1,3);    
           3*R(2,2)*R(2,1); 
           3*R(2,2)*R(2,3);
           3*R(3,3)*R(3,1);     
           3*R(3,3)*R(3,2);    
           R(1,1)*R(2,2)+2*R(1,2)^2;  
           R(1,1)*R(3,3)+2*R(1,3)^2;        
           R(2,2)*R(3,3)+2*R(2,3)^2;
           R(1,1)*R(2,3)+2*R(1,2)*R(1,3);   
           R(2,2)*R(1,3)+2*R(2,1)*R(2,3);  
           R(3,3)*R(1,2)+2*R(3,1)*R(3,2)];

D=U*diag([sigma(1) sigma(2)*[1 1]])*U';
for n=1:N_basis
        R=D/(2*pi^2);       
        mu=2*pi*R*qg_basis(n,:)';
        const=2*exp(-2*pi^2*qg_basis(n,:)*R*qg_basis(n,:)');
        
        M=R-(mu*mu');
        M2_basis(:,n)=const*[M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3)]';
        M4_basis(:,n)=const*[mu(1)^4-6*mu(1)^2*R(1,1)+3*R(1,1)^2; 
                            mu(2)^4-6*mu(2)^2*R(2,2)+3*R(2,2)^2; 
                            mu(3)^4-6*mu(3)^2*R(3,3)+3*R(3,3)^2;
                            mu(1)^3*mu(2)-3*mu(1)*mu(2)*R(1,1)-3*mu(1)^2*R(1,2)+3*R(1,1)*R(1,2);  
                            mu(1)^3*mu(3)-3*mu(1)*mu(3)*R(1,1)-3*mu(1)^2*R(1,3)+3*R(1,1)*R(1,3);
                            mu(2)^3*mu(1)-3*mu(2)*mu(1)*R(2,2)-3*mu(2)^2*R(2,1)+3*R(2,2)*R(2,1);  
                            mu(2)^3*mu(3)-3*mu(2)*mu(3)*R(2,2)-3*mu(2)^2*R(2,3)+3*R(2,2)*R(2,3); 
                            mu(3)^3*mu(1)-3*mu(3)*mu(1)*R(3,3)-3*mu(3)^2*R(3,1)+3*R(3,3)*R(3,1);  
                            mu(3)^3*mu(2)-3*mu(3)*mu(2)*R(3,3)-3*mu(3)^2*R(3,2)+3*R(3,3)*R(3,2);
                            mu(1)^2*mu(2)^2-mu(2)^2*R(1,1)-mu(1)^2*R(2,2)-4*mu(1)*mu(2)*R(1,2)+R(1,1)*R(2,2)+2*R(1,2)^2; 
                            mu(1)^2*mu(3)^2-mu(3)^2*R(1,1)-mu(1)^2*R(3,3)-4*mu(1)*mu(3)*R(1,3)+R(1,1)*R(3,3)+2*R(1,3)^2; 
                            mu(2)^2*mu(3)^2-mu(3)^2*R(2,2)-mu(2)^2*R(3,3)-4*mu(2)*mu(3)*R(2,3)+R(2,2)*R(3,3)+2*R(2,3)^2;
                            mu(1)^2*mu(2)*mu(3)-mu(2)*mu(3)*R(1,1)-2*mu(1)*mu(3)*R(1,2)-2*mu(1)*mu(2)*R(1,3)-mu(1)^2*R(2,3)+R(1,1)*R(2,3)+2*R(1,2)*R(1,3);
                            mu(2)^2*mu(1)*mu(3)-mu(1)*mu(3)*R(2,2)-2*mu(2)*mu(3)*R(2,1)-2*mu(2)*mu(1)*R(2,3)-mu(2)^2*R(1,3)+R(2,2)*R(1,3)+2*R(2,1)*R(2,3);
                            mu(3)^2*mu(1)*mu(2)-mu(1)*mu(2)*R(3,3)-2*mu(3)*mu(2)*R(3,1)-2*mu(3)*mu(1)*R(3,2)-mu(3)^2*R(1,2)+R(3,3)*R(1,2)+2*R(3,1)*R(3,2);
                            ];  
end

end
