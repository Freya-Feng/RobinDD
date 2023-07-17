function [a,b,m,c,r,bbx,bby,fr,g,fn,RB] = mystokes_q1p0(xy,xyp,mv,ev,ebound)
%STOKES_Q1P0 Q1-P0 matrix generator
%   [A,B,Q,C,G,Bx,By,f,g] = stokes_q1p0(xy,xyp,mv,ev);
%   input
%          xy         Q2 nodal coordinate vector
%          xyp        Q1 nodal coordinate vector
%          mv         Q2 element mapping matrix
%          ev         element mapping matrix
%   output
%          A          Q1 vector diffusion matrix
%          B          Q1-P0 divergence matrix
%          Q          P0 mass matrix
%          C          pressure stabilization matrix
%          G          Q1 vector mass matrix
%          Bx         Q1 x-derivative matrix
%          By         Q1 y-derivative matrix
%          f          velocity rhs vector
%          g          pressure rhs vector
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function flowbc.
%   IFISS function:
% Copyright (c) 2023 Yani Feng, Qifeng Liao
nngpt=4;
x=xy(:,1); y=xy(:,2);
xp=xyp(:,1); yp=xyp(:,2);
nvtx=length(x); nu=2*nvtx; np=length(xp);
nel=length(ev(:,1)); mp=[1:nel]';
lx=max(x)-min(x); ly=max(y)-min(y);
hx=max(diff(x)); hy=max(diff(y));
% the size of all the elements (hxe,hye)
hxe=zeros(nel,1);
hye=zeros(nel,1);
for i=1:nel
    nodes=ev(i,:);
    coords=xy(nodes,:);
    hxe(i)=coords(3,1)-coords(1,1);
    hye(i)=coords(3,2)-coords(1,2);
end
fprintf('setting up Q1-P0 matrices...  ')
%
% initialise global matrices
a = sparse(nu,nu);
r = sparse(nu,nu);
bbx = sparse(nvtx,nvtx);
bby = sparse(nvtx,nvtx);
bx = sparse(np,nvtx);
by = sparse(np,nvtx);
b = sparse(np,nu);
m = sparse(np,np);
frxe= zeros(nel,4);
frye= zeros(nel,4);
frx = zeros(nvtx,1);% right hand f
fry = zeros(nvtx,1);
g = zeros(np,1);
%
% B: Robin bc block
edge=4;
for k=1:edge
    RB{k} = sparse(nu,nu);
    fn{k} = zeros(nu,1);% neumann bc f
    be{k} = zeros(nel,4,4);
    fne{k}= zeros(nel,4);
end
%
% Gauss point integration rules
if (nngpt==4)        % 2x2 Gauss points
    gpt=1.0e0/sqrt(3.0e0);
    s(1) = -gpt; t(1) = -gpt; wt(1)=1;
    s(2) =  gpt; t(2) = -gpt; wt(2)=1;
    s(3) =  gpt; t(3) =  gpt; wt(3)=1;
    s(4) = -gpt; t(4) =  gpt; wt(4)=1;
elseif (nngpt==1)   % 1x1 Gauss point
    s(1) =    0; t(1) =    0; wt(1)=4;
else
    error('Check Gauss point integration specification')
end
%
% set up 1D Gauss points
[oneg,onew]=gausspoints_oned(2);
% inner loop over elements
for ivtx = 1:4
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
end
ae = zeros(nel,4,4);
re = zeros(nel,4,4);
bbxe = zeros(nel,4,4);
bbye = zeros(nel,4,4);
bxe = zeros(nel,1,4);
bye = zeros(nel,1,4);
mpe = zeros(nel,1,1);
%
% loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    wght=wt(igpt);
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
    rhs = gauss_source_rb(sigpt,tigpt,xl_v,yl_v);
    for j = 1:4
        for i = 1:4
            ae(:,i,j)  = ae(:,i,j)  + wght*dphidx(:,i).*dphidx(:,j).*invjac(:);
            ae(:,i,j)  = ae(:,i,j)  + wght*dphidy(:,i).*dphidy(:,j).*invjac(:);
            re(:,i,j)  = re(:,i,j)  + wght*phi(:,i).*phi(:,j).*jac(:);
            bbxe(:,i,j) = bbxe(:,i,j) - wght*phi(:,i) .*dphidx(:,j);
            bbye(:,i,j) = bbye(:,i,j) - wght*phi(:,i) .*dphidy(:,j);
        end
        bxe(:,1,j) = bxe(:,1,j) - wght* dphidx(:,j);
        bye(:,1,j) = bye(:,1,j) - wght* dphidy(:,j);
        frxe(:,j) = frxe(:,j)  + wght*rhs(1:nel) .* phi(:,j) .* jac(:);
        frye(:,j) = frye(:,j)  + wght*rhs(nel+1:end) .* phi(:,j) .* jac(:);
    end
    mpe(:,1,1) = mpe(:,1,1) + wght*jac(:);
    % end of Gauss point loop
end
%
% boundary condition
for igpt=1:2
    s=oneg(igpt);
    w=onew(igpt);
    for k=1:edge
        if k==1
            [jac,invjac,phi,dphidx,dphidy] = deriv(s,-1,xl_v,yl_v);
            jac_oned=hxe/2;
            t=find(ebound(:,2)==1);
        elseif k==2
            [jac,invjac,phi,dphidx,dphidy] = deriv(1,s,xl_v,yl_v);
            jac_oned=hye/2;
            t=find(ebound(:,2)==2);
        elseif k==3
            [jac,invjac,phi,dphidx,dphidy] = deriv(s,1,xl_v,yl_v);
            jac_oned=hxe/2;
            t=find(ebound(:,2)==3);
        elseif k==4
            [jac,invjac,phi,dphidx,dphidy] = deriv(-1,s,xl_v,yl_v);
            jac_oned=hye/2;
            t=find(ebound(:,2)==4);
        end
        ebc=ebound(t,1);
        % Robin bc
        for j=1:4
            for i=1:4
                be{k}(ebc,i,j)=be{k}(ebc,i,j)+w*phi(ebc,i).*phi(ebc,j).*jac_oned(ebc);
            end
        end
        % Neumann bc \partial u/\partial n=1
        for j=1:4
            fne{k}(ebc,j) = fne{k}(ebc,j)  + w*phi(ebc,j) .* jac_oned(ebc);
        end
    end
end
% element assembly into global matrices
% component velocity matrices ...
for krow=1:4
    nrow=ev(:,krow);
    for kcol=1:4
        ncol=ev(:,kcol);
        a = a + sparse(nrow,ncol,ae(:,krow,kcol),nu,nu);
        a = a + sparse(nrow+nvtx,ncol+nvtx,ae(:,krow,kcol),nu,nu);
        r = r + sparse(nrow,ncol,re(:,krow,kcol),nu,nu);
        r = r + sparse(nrow+nvtx,ncol+nvtx,re(:,krow,kcol),nu,nu);
        bbx = bbx + sparse(nrow,ncol,bbxe(:,krow,kcol),nvtx,nvtx);
        bby = bby + sparse(nrow,ncol,bbye(:,krow,kcol),nvtx,nvtx);
        for k=1:edge
            RB{k} = RB{k} + sparse(nrow,ncol,be{k}(:,krow,kcol),nu,nu);
            RB{k} = RB{k} + sparse(nrow+nvtx,ncol+nvtx,be{k}(:,krow,kcol),nu,nu);
        end
    end
    
    for k=1:edge
        fn{k}(nrow,1) = fn{k}(nrow,1) + fne{k}(:,krow);
        fn{k}(nrow+nvtx,1)=fn{k}(nrow+nvtx,1)+fne{k}(:,krow);
    end
    %right hand f
    frx(nrow,1) = frx(nrow,1) + frxe(:,krow);
    fry(nrow,1) = fry(nrow,1) + frye(:,krow);
    kcol=1;
    ncol=mp;
    bx = bx + sparse(ncol,nrow,bxe(:,kcol,krow),np,nvtx);
    by = by + sparse(ncol,nrow,bye(:,kcol,krow),np,nvtx);
end
fr=[frx;fry];
%
% vector velocity matrices ...
b = [bx,by];
%
% pressure matrices ...
m = sparse(mp,mp,mpe(:,1,1),np,np);
%



% stabilisation matrix
mel=length(mv(:,1));
cm=zeros(mel,4,4); hm=zeros(mel,1);
c=sparse(np,np);
% loop over macroelements
for ivtx = 1:9
    xlm(:,ivtx) = x(mv(:,ivtx));
    ylm(:,ivtx) = y(mv(:,ivtx));
end
elarea=full(diag(m));
hm(:)=sum(reshape(elarea,4,mel))';
hm=hm/4;
%      le=1/he;
%      cm=0.25*hm*[ le41+le12,     -le12,         0,     -le41;
%  	               -le12, le12+le23,     -le23,         0;
%                      0,     -le23, le23+le34,     -le34;
%	 			    -le41,         0,     -le34, le34+le41];
%
xe=xlm(:,9)-xlm(:,5);ye=ylm(:,9)-ylm(:,5);le12=hm;
cm(:,1,1) =cm(:,1,1)+ le12; cm(:,1,2)=-le12;
cm(:,2,2) =cm(:,2,2)+ le12; cm(:,2,1)=-le12;
xe=xlm(:,9)-xlm(:,6);ye=ylm(:,9)-ylm(:,6);le23=hm;
cm(:,2,2) =cm(:,2,2)+ le23; cm(:,2,3)=-le23;
cm(:,3,3) =cm(:,3,3)+ le23; cm(:,3,2)=-le23;
xe=xlm(:,9)-xlm(:,7);ye=ylm(:,9)-ylm(:,7);le34=hm;
cm(:,3,3) =cm(:,3,3)+ le34; cm(:,3,4)=-le34;
cm(:,4,4) =cm(:,4,4)+ le34; cm(:,4,3)=-le34;
xe=xlm(:,9)-xlm(:,8);ye=ylm(:,9)-ylm(:,8);le41=hm;
cm(:,4,4) =cm(:,4,4)+ le41; cm(:,4,1)=-le41;
cm(:,1,1) =cm(:,1,1)+ le41; cm(:,1,4)=-le41;
%
%  macroelement assembly into global matrices
for krow=1:4
    nrow=[0:4:nel-4]+krow;
    for kcol=1:4
        ncol=[0:4:nel-4]+kcol;
        c = c + sparse(nrow,ncol,cm(:,krow,kcol),np,np);
    end
end
%
fprintf('done\n')
return
