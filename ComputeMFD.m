function MFD=ComputeMFD(D,V,mask,FileName,sigma)
if(nargin<5)
    sigma=[0.0015 0.0008];
end

[~,~,nx,ny]=size(D);


mfd=zeros(nx,ny);
for ix=1:nx
    parfor iy=1:ny
        MFD_basis=ConstructMFDBasis(D(:,:,ix,iy),sigma);
        Vec=MFD_basis*V(:,ix,iy);
        mfd(ix,iy)=sum(Vec);
    end
end
MFD=zeros(size(mask));
MFD(mask~=0)=mfd;

%%
if(nargin>=4&& (~isempty(FileName)))
% sd=[2.500000,-0.000000,0.000000; -0.000000,2.500000,0.000000; -0.000000,-0.000000,2.500000]; 
% mat2nhdr(MFD, FileName, 'custom',sd);
nii=make_nii(MFD);
save_nii(nii,[FileName,'_MFD.nii']);
end




end


function MFD_basis=ConstructMFDBasis(D0,sigma)
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

MFD_basis=zeros(9,N_column);

R=D0/(2*pi^2);
ID= [1    10    11    10     2    12    11    12     3];

m=[3*R(1,1)^2;            
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
MFD_basis(:,end)=m(ID);

D=U*diag([sigma(1) sigma(2)*[1 1]])*U';
for n=1:N_basis
        R=D/(2*pi^2);       
        mu=2*pi*R*qg_basis(n,:)';
        const=2*exp(-2*pi^2*qg_basis(n,:)*R*qg_basis(n,:)');

        m=const*[mu(1)^4-6*mu(1)^2*R(1,1)+3*R(1,1)^2; 
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
        MFD_basis(:,n)=m(ID);
end

end