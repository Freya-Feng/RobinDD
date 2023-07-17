function [u,p]=mystokes_full_system_rb(start_data,viscosity)
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
f=start_data.f;
g=start_data.g;
bound=start_data.bound;
xy=start_data.xy;
nvtx=length(xy(:,1));
[Ast,Bst,fst,gst] = myflowbc_rb(viscosity*A,B,f,g,xy,bound,1);
np=length(gst);
%% compute solution
qmethod=1;
%beta=default('stabilization parameter (default is 1/4)',1/4);
if qmethod==1
    beta=1/4*viscosity;
    C=start_data.C;
%     // ev=start_data.ev;
%     // % get jac
%     // for ivtx = 1:4
%     //     xl_v(:,ivtx) = xy(ev(:,ivtx),1);
%     //     yl_v(:,ivtx) = xy(ev(:,ivtx),2);
%     // end
%     // jac=deriv(0,0,xl_v,yl_v);
%     // %solve system with zero mean
%     // xst_larg=[Ast,Bst',zeros(2*nvtx,1);Bst,-beta*C,zeros(np,1);zeros(1,2*nvtx),jac',0]\[fst;gst;0];
%     // xst=xst_larg(1:end-1);
   xst=[Ast Bst';Bst,-beta*C]\[fst;gst];
    % else
    %     beta=0;
    %     xst=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
end
u=xst(1:2*nvtx);p=xst(2*nvtx+1:end);
end