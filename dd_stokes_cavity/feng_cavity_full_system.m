function [u,p]=feng_cavity_full_system(start_data,viscosity)
% obtain solutions for the whole domain
% input 
%         start_data  Galerkin system data
%         viscosity   incompressible flow parameter
% output 
%         u           velocity in the whole domain
%         p           pressure in the whole domain 
%   IFISS function: YF; 1 May 2023.
% Copyright (c) 2023 Yani Feng, Qifeng Liao
A=start_data.A;
B=start_data.B;
g=start_data.g;
bound=start_data.bound;
xy=start_data.xy;
nvtx=length(xy(:,1));
f0=zeros(2*nvtx,1);
[Ast,Bst,fst,gst] = flowbc(viscosity*A,B,f0,g,xy,bound);
%% compute solution
np=length(gst); nu=length(f0)/2;
qmethod=1;
if qmethod==1
    beta=1/4*viscosity;
    C=start_data.C;
    % stabilized version
%     xstz=[Ast,Bst',zeros(2*nu,1);Bst,-beta*C,ones(np,1)/np; ...
%         zeros(1,2*nu),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
%     xst=xstz(1:end-1); multiplier=xstz(end);
    xst=[Ast Bst';Bst,-beta*C]\[fst;gst];
end
u=xst(1:2*nvtx);p=xst(2*nvtx+1:end);
end