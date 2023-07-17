function fn=myupdateg_q1p0(xy,ev,ebound,g,k0)
% compute the inner product betwwen g and basis functionc on the interface
% input 
%         xy        nodal coordinate vector 
%         ev        element vertex matrix
%         ebound    element boundary edge matrix  
%         g         Robin interface conditions
%         k0        the number of boundary (k0=1/2/3/4 for bottom/right/top/left boundary)
% output 
%         fn        the inner product betwwen g and basis function on the interface
%   IFISS function: YF; 1 May 2023.
% Copyright (c) 2023 Yani Feng, Qifeng Liao
x=xy(:,1); y=xy(:,2);
miny=min(xy(:,2));
dy=(max(xy(:,2))-min(xy(:,2)))/(length(unique(xy(:,2)))-1);
nvtx=length(x); nu=2*nvtx;
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
edge=4;
for k=1:edge
    fn{k} = zeros(nu,1);% neumann bc f
    fnxe{k}= zeros(nel,4);
    fnye{k}= zeros(nel,4);
end
% set up 1D Gauss points
[oneg,onew]=gausspoints_oned(2);
% inner loop over elements
for ivtx = 1:4
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
end
% boundary condition
for igpt=1:2
    s=oneg(igpt);
    w=onew(igpt);
    k=k0;
        if k==2
            [jac,invjac,phi,dphidx,dphidy] = deriv(1,s,xl_v,yl_v);
            jac_oned=hye/2;
            t=find(ebound(:,2)==2);
            ebc=ebound(t,1);
            [rhgx,rhgy] = g_source_rb(1,s,xl_v(ebc,:),yl_v(ebc,:),g,dy,miny);
        elseif k==4
            [jac,invjac,phi,dphidx,dphidy] = deriv(-1,s,xl_v,yl_v);
            jac_oned=hye/2;
            t=find(ebound(:,2)==4);
            ebc=ebound(t,1);
            [rhgx,rhgy] = g_source_rb(-1,s,xl_v(ebc,:),yl_v(ebc,:),g,dy,miny);
        end
        % robin bc
        for j=1:4
            fnxe{k}(ebc,j) = fnxe{k}(ebc,j)  + w*rhgx(:).*phi(ebc,j) .* jac_oned(ebc);
            fnye{k}(ebc,j) = fnye{k}(ebc,j)  + w*rhgy(:).*phi(ebc,j) .* jac_oned(ebc);
        end
end
for krow=1:4
    nrow=ev(:,krow);
    for k=1:edge
        fn{k}(nrow,1)=fn{k}(nrow,1) +fnxe{k}(:,krow);
        fn{k}(nrow+nvtx,1)=fn{k}(nrow+nvtx,1)+fnye{k}(:,krow);
    end
end

end
