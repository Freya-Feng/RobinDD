function [e,ur]= sc_fem_stokes_compare(u,xy,ucs,xycs)
% input 
%         u    global velocity solution  
%         xy   nodal coordinate vector corresponding to u
%         ucs  domain decomposition velocity solution
%         xycs nodal coordinate vector corresponding to ucs
% output 
%         e      error between global velocity solution and domain decomposition velocity solution
%         ur     domain decomposition velocity solution corresponding to xy
%   IFISS function: YF; 1 May 2023.
% Copyright (c) 2023 Yani Feng, Qifeng Liao
nvtx=length(u)/2;
tol=10^(-10);
ur=0*u;
for i=1:nvtx
     % find static condenstion solution wich match this point
    t=min(find(abs(xycs(:,1)-xy(i,1))<tol&abs(xycs(:,2)-xy(i,2))<tol));
    ur(i)=ucs(t);
    ur(i+nvtx)=ucs(t+nvtx);
    e(i)=abs(u(i)-ucs(t));
    e(i+nvtx)=abs(u(i+nvtx)-ucs(t+nvtx));
end
return

