function [u2,p2,u2bc]=feng_cavity_dd_right(start_data,g2,gamma2,viscosity)
% obtain solutions for right domain
% input 
%         start_data  Galerkin system data
%         g2          Robin interface conditions
%         gamma2      Robin interface parameters
%         viscosity   incompressible flow parameter
% output 
%         u2          velocity in right domain
%         p2          pressure in right domain 
%         u2bc        velocity on the interface
%   IFISS function: YF; 1 May 2023.
% Copyright (c) 2023 Yani Feng, Qifeng Liao
A=start_data.A;
B=start_data.B;
g=start_data.g;
RB=start_data.RB;
xy=start_data.xy;
port=start_data.port;
bound=start_data.bound;
ev=start_data.ev;
ebound=start_data.ebound;
nvtx=length(xy(:,1));
k0=4;
fng=myupdateg_q1p0(xy,ev,ebound,g2,k0);
A0=viscosity*A+gamma2*RB{4};f=fng{k0};
%% boundary conditions
pb=unique([setdiff(bound,port);port(1);port(end)]);
[Ast,Bst,fst,gst] = flowbc(A0,B,f,g,xy,pb);
%np=length(gst);
qmethod=1;
if qmethod==1
    beta=1/4*viscosity;
    C=start_data.C;
    xst=[Ast,Bst';Bst,-beta*C]\[fst;gst];
end
u2=xst(1:2*nvtx);p2=xst(2*nvtx+1:end);
%% compute comunication data mu1
u2bc=[u2(port);u2(port+nvtx)];
end