function [GK,NIG]=ComputeKurtosis_Qspace(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

gk=zeros(nx,ny);
nig=zeros(nx,ny);

for ix=1:nx
    for iy=1:ny
        [M2_basis,M4_basis]=ConstructM2M4_Qspace(D(:,:,ix,iy),sigma);
        m2=M2_basis*V(:,ix,iy);
        m4=M4_basis*V(:,ix,iy);
        RTOP_basis=ConstructRTOPBasis(D(:,:,ix,iy),sigma);
        rtop=RTOP_basis*V(:,ix,iy);
        M2=m2([1  4  5;
               4  2  6; 
               5  6  3])/rtop;
        M4=m4([ 1   4   5   4   10  13  5   13  11;
                4   10  13  10  6   14  13  14  15;
                5   13  11  13  14  15  11  15  8;
                4   10  13  10  6   14  13  14  15;
                10  6   14  6   2   7   14  7   12;
                13  14  15  14  7   12  15  12  9;
                5   13  11  13  14  15  11  15  8;
                13  14  15  14  7   12  15  12  9;
                11  15  8   15  12  9   8   9   3])/rtop;
        [v,d]=eig(M4);d=diag(d);d(d<0)=0;M4=v*diag(d)*v';
        [v,d]=eig(M2);d=diag(d);d(d<0)=0;M2=v*diag(d)*v';
        
        R=pinv(M2);r=vec(R);
        gk(ix,iy)=r'*M4*r;
        nig(ix,iy)=trace(M4)/(trace(M2)^2);
          
    end
end

GK=zeros(size(mask));
GK(mask~=0)=gk;

NIG=zeros(size(mask));
NIG(mask~=0)=nig;


if(nargin>=4&& (~isempty(FileName)))
sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
mat2nhdr(GK, strcat(FileName,'_GK'), 'custom',sd);
mat2nhdr(NIG, strcat(FileName,'_NIG'), 'custom',sd);
end


end


function [M2_basis,M4_basis]=ConstructM2M4_Qspace(D0,sigma)
t = 70*1e-3;
b_basis=[2 4]*1000; % location of basis functions
[ug,~]=icosahedron(2);
N_ug=length(ug);
ug=ug(1:N_ug/2,:);
N_ug=N_ug/2;
q_basis=(1/(2*pi))*sqrt(b_basis/t);
qg_basis=kron(diag(q_basis),eye(N_ug))*repmat(ug,[length(b_basis),1]);
N_basis=length(qg_basis);

N_column=N_basis+1;
[U,~,~]=svd(D0);

M2_basis=zeros(6,N_column);%coefficients in the covariance matrix
M4_basis=zeros(15,N_column);


Rq=inv(D0)/2;
const=(2*pi)^(1.5)*sqrt(det(Rq));
M2_basis(:,end)=const*[Rq(1,1) Rq(2,2) Rq(3,3) Rq(1,2) Rq(1,3) Rq(2,3)]';
M4_basis(:,end)=const*[3*Rq(1,1)^2;            
                       3*Rq(2,2)^2;           
                       3*Rq(3,3)^2; 
                       3*Rq(1,1)*Rq(1,2);     
                       3*Rq(1,1)*Rq(1,3);    
                       3*Rq(2,2)*Rq(2,1); 
                       3*Rq(2,2)*Rq(2,3);
                       3*Rq(3,3)*Rq(3,1);     
                       3*Rq(3,3)*Rq(3,2);    
                       Rq(1,1)*Rq(2,2)+2*Rq(1,2)^2;  
                       Rq(1,1)*Rq(3,3)+2*Rq(1,3)^2;        
                       Rq(2,2)*Rq(3,3)+2*Rq(2,3)^2;
                       Rq(1,1)*Rq(2,3)+2*Rq(1,2)*Rq(1,3);   
                       Rq(2,2)*Rq(1,3)+2*Rq(2,1)*Rq(2,3);  
                       Rq(3,3)*Rq(1,2)+2*Rq(3,1)*Rq(3,2)];

D=U*diag([sigma(1) sigma(2)*[1 1]])*U';
for n=1:N_basis
        Rq=inv(D)/(2);       
        mu=qg_basis(n,:)';
        const=2*(2*pi)^(1.5)*sqrt(det(Rq));
        
        M=Rq+(mu*mu');
        M2_basis(:,n)=const*[M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3)]';        
        M4_basis(:,n)=const*[mu(1)^4+6*mu(1)^2*Rq(1,1)+3*Rq(1,1)^2; 
                              mu(2)^4+6*mu(2)^2*Rq(2,2)+3*Rq(2,2)^2; 
                              mu(3)^4+6*mu(3)^2*Rq(3,3)+3*Rq(3,3)^2;
                              mu(1)^3*mu(2)+3*mu(1)*mu(2)*Rq(1,1)+3*mu(1)^2*Rq(1,2)+3*Rq(1,1)*Rq(1,2);  
                              mu(1)^3*mu(3)+3*mu(1)*mu(3)*Rq(1,1)+3*mu(1)^2*Rq(1,3)+3*Rq(1,1)*Rq(1,3);
                              mu(2)^3*mu(1)+3*mu(2)*mu(1)*Rq(2,2)+3*mu(2)^2*Rq(2,1)+3*Rq(2,2)*Rq(2,1);  
                              mu(2)^3*mu(3)+3*mu(2)*mu(3)*Rq(2,2)+3*mu(2)^2*Rq(2,3)+3*Rq(2,2)*Rq(2,3); 
                              mu(3)^3*mu(1)+3*mu(3)*mu(1)*Rq(3,3)+3*mu(3)^2*Rq(3,1)+3*Rq(3,3)*Rq(3,1);  
                              mu(3)^3*mu(2)+3*mu(3)*mu(2)*Rq(3,3)+3*mu(3)^2*Rq(3,2)+3*Rq(3,3)*Rq(3,2);
                              mu(1)^2*mu(2)^2+mu(2)^2*Rq(1,1)+mu(1)^2*Rq(2,2)+4*mu(1)*mu(2)*Rq(1,2)+Rq(1,1)*Rq(2,2)+2*Rq(1,2)^2; 
                              mu(1)^2*mu(3)^2+mu(3)^2*Rq(1,1)+mu(1)^2*Rq(3,3)+4*mu(1)*mu(3)*Rq(1,3)+Rq(1,1)*Rq(3,3)+2*Rq(1,3)^2; 
                              mu(2)^2*mu(3)^2+mu(3)^2*Rq(2,2)+mu(2)^2*Rq(3,3)+4*mu(2)*mu(3)*Rq(2,3)+Rq(2,2)*Rq(3,3)+2*Rq(2,3)^2;
                              mu(1)^2*mu(2)*mu(3)+mu(2)*mu(3)*Rq(1,1)+2*mu(1)*mu(3)*Rq(1,2)+2*mu(1)*mu(2)*Rq(1,3)+mu(1)^2*Rq(2,3)+Rq(1,1)*Rq(2,3)+2*Rq(1,2)*Rq(1,3);
                              mu(2)^2*mu(1)*mu(3)+mu(1)*mu(3)*Rq(2,2)+2*mu(2)*mu(3)*Rq(2,1)+2*mu(2)*mu(1)*Rq(2,3)+mu(2)^2*Rq(1,3)+Rq(2,2)*Rq(1,3)+2*Rq(2,1)*Rq(2,3);
                              mu(3)^2*mu(1)*mu(2)+mu(1)*mu(2)*Rq(3,3)+2*mu(3)*mu(2)*Rq(3,1)+2*mu(3)*mu(1)*Rq(3,2)+mu(3)^2*Rq(1,2)+Rq(3,3)*Rq(1,2)+2*Rq(3,1)*Rq(3,2);
                              ]; 
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
R=D0/(2*pi^2);
RTOP_basis(end)=pi^(1.5)/sqrt(det(D0));

RTOP_basis(1:N_basis)=2*pi^(1.5)/sqrt(sigma(1)*sigma(2)^2)*ones(1,N_basis);

end

