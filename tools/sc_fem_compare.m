function [e,pr]=sc_fem_compare(p,xyp,pcs,xypcs)
% input 
%         p    global pressure solution  
%         xyp   nodal coordinate vector corresponding to p
%         pcs  domain decomposition pressure solution
%         xypcs nodal coordinate vector corresponding to pcs
% output 
%         e      error between global pressure solution and domain decomposition pressure solution
%         pr     domain decomposition pressure solution corresponding to xyp
%   IFISS function: YF; 1 May 2023.
% Copyright (c) 2023 Yani Feng, Qifeng Liao
nvtx=length(p);
tol=10^(-10);
pr=0*p;
for i=1:nvtx
     % find static condenstion solution wich match this point
    t=min(find(abs(xypcs(:,1)-xyp(i,1))<tol&abs(xypcs(:,2)-xyp(i,2))<tol));
    pr(i)=pcs(t);
    e(i)=abs(p(i)-pcs(t));
end
return

