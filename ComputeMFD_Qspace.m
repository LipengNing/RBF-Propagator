function MFD=ComputeMFD_Qspace(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);

mfd=zeros(nx,ny);
for ix=1:nx
    parfor iy=1:ny
        MFD_basis=ConstructM4Basis_Qspace(D(:,:,ix,iy),sigma);
        Vec=MFD_basis*V(:,ix,iy);
        Vec(Vec<0)=0;
        mfd(ix,iy)=sum(Vec);
    end
end
MFD=zeros(size(mask));
MFD(mask~=0)=mfd/max(mfd(:))*1000;

%%
if(nargin>=4&& (~isempty(FileName)))
% sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
% mat2nhdr(MFD, FileName, 'custom',sd);
nii=make_nii(MFD);
save_nii(nii,[FileName,'_QMFD.nii']);
end

end


function M4_basis=ConstructM4Basis_Qspace(D0,sigma)
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


M4_basis=zeros(9,N_column);
ID=[ 1    10    11    10     2    12    11    12     3];

Rq=inv(D0)/2;
const=(2*pi)^(1.5)*sqrt(det(Rq));
m=const*[3*Rq(1,1)^2;            
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

M4_basis(:,end)=m(ID);

               
               
D=U*diag([sigma(1) sigma(2)*[1 1]])*U';
for n=1:N_basis
        Rq=inv(D)/(2);       
        mu=qg_basis(n,:)';
        const=2*(2*pi)^(1.5)*sqrt(det(Rq));
        
        M=Rq+(mu*mu');
   
        m=const*[mu(1)^4+6*mu(1)^2*Rq(1,1)+3*Rq(1,1)^2; 
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
            M4_basis(:,n)=m(ID);
end

end