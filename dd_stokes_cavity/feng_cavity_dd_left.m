function [u1,p1,u1bc]=feng_cavity_dd_left(start_data,g1,gamma1,viscosity)
% obtain solutions for left domain
% input 
%         start_data  Galerkin system data
%         g1          Robin interface conditions
%         gamma1      Robin interface parameters
%         viscosity   incompressible flow parameter
% output 
%         u1          velocity in left domain
%         p1          pressure in left domain 
%         u1bc        velocity on the interface
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
k0=2;
fng=myupdateg_q1p0(xy,ev,ebound,g1,k0);
A0=viscosity*A+gamma1*RB{2};f=fng{k0};
%% boundary conditions
pb=unique([setdiff(bound,port);port(1);port(end)]);
[Ast,Bst,fst,gst] = flowbc(A0,B,f,g,xy,pb);
%np=length(gst);nu=length(f)/2;
%% compute solution
qmethod=1;
%beta=default('stabilization parameter (default is 1/4)',1/4);
if qmethod==1
    beta=1/4*viscosity;
    C=start_data.C;
    xst=[Ast Bst';Bst,-beta*C]\[fst;gst];
end
%% compute comunication data mu1
u1=xst(1:2*nvtx);p1=xst(2*nvtx+1:end);
u1bc=[u1(port);u1(port+nvtx)];
end