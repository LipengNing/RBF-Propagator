function [u,fcs] = icosahedron(level)
% ICOSAHEDRON Sampling the unit sphere  
%  ICOSAHEDRON(level) samples the unit sphere using a regular icosahedron
% tessellation of the order defined by 'level'.
%
%                     [u,fcs] = icosahedron(level)
% Input:
%             level - tessellation order
% Output:
%                 u - N-by-3 matrix of samples of the unit sphere
%               fcs - vertices-to-faces correspondence matrix
%
% written by Oleg Michailovich, February 5th, 2009

C=1/sqrt(1.25);
t=(2*pi/5)*(0:4)';
u1=C*[cos(t) sin(t) 0.5*ones(5,1)];
u2=C*[cos(t+0.2*pi) sin(t+0.2*pi) -0.5*ones(5,1)];
u=[[0 0 1]; u1; u2; [0 0 -1]];

if (level>0),
    for lev=1:level,
        fcs=convhulln(u);
        N=size(fcs,1);
        U=zeros(3*N,3);
        for k=1:N,
            A=u(fcs(k,1),:);
            B=u(fcs(k,2),:);
            C=u(fcs(k,3),:);
            U(3*k-2:3*k,:)=0.5*[A+B; B+C; A+C];
        end
        U=unique(U,'rows');
        U=U./repmat(sqrt(sum(U.^2,2)),[1 3]);
        u=[u; U]; %#ok<AGROW>
    end    
    [C,ind]=sort(u(:,3),1,'descend');
    u=u(ind,:);
    index=find(u(:,3)==0);
    v=u(index,:);
    [C,ind]=sort(v(:,2),1,'descend');
    u(index,:)=v(ind,:);
end

u=u./repmat(sqrt(sum(u.^2,2))+eps,[1 3]);
fcs=convhulln(u);