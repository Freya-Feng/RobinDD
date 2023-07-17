function [u2,p2,u2bc]=feng_navier_dd_right_rb(start_data,h2,gamma2,viscosity,nlmethod,tol_nl)
% obtain solutions for right domain
% input 
%         start_data  Galerkin system data
%         g2          Robin interface conditions
%         gamma2      Robin interface parameters
%         viscosity   incompressible flow parameter
%         nlmethod    Picard/Newton/hybrid linearization 1/2/3 (default hybrid)
%         tol_nl      nonlinear tolerance 
% output 
%         u2          velocity in right domain
%         p2          pressure in right domain 
%         u2bc        velocity on the interface
%   IFISS function: YF; 1 May 2023.
% Copyright (c) 2023 Yani Feng, Qifeng Liao
A=start_data.A;
B=start_data.B;
fr=start_data.f;
fn=start_data.fn;
g=start_data.g;
RB=start_data.RB;
xy=start_data.xy;
port=start_data.port;
bound=start_data.bound;
ebound=start_data.ebound;
ev=start_data.ev;
%ebound=start_data.ebound;
nvtx=length(xy(:,1));
%%
% fprintf('Incompressible flow problem on step domain ...\n')
% viscosity=1/200;
% nlmethod=3;
%viscosity=default('viscosity parameter (default 1/50)',1/50);
% nlmethod=default('Picard/Newton/hybrid linearization 1/2/3 (default hybrid)',3);
nlmethod=nlmethod-1;
if nlmethod==0,
%     maxit_p=default('number of Picard iterations (default 9)',9);
    maxit_p=9;
    maxit_n=0;
elseif nlmethod==1,
    maxit_p=0;
    maxit_n=6;
%     maxit_n=default('number of Newton iterations (default 6)',6);
else
    maxit_p=4;
    maxit_n=6;
%     maxit_p=default('number of Picard iterations (default 2)',2);
%     maxit_n=default('number of Newton iterations (default 4)',4);
end
% tol_nl=default('nonlinear tolerance (default 1.e-8)',1.e-8);
%tol_nl=1.e-10;
%% initialize for nonlinear iteration: compute Stokes solution
% A0=A+gamma2*RB{4};f=fn{4}.*neu;
k0=4;
fng=myupdateg_q1p0(xy,ev,ebound,h2,k0);
A0=A+gamma2*RB{4};
f=fng{k0};
%% boundary conditions this bound include port
pb=unique([bound(find(xy(bound,2)==-1));bound(find(xy(bound,2)==1))]);
[Ast,Bst,fst,gst] = flowbc(A0,B,f,g,xy,pb);
nlres0_norm = norm([fst;gst]);
%
nv=length(fst)/2; np=length(gst);
%% compute initial solution
qmethod=1;
if qmethod==1
    beta=1/4*viscosity;
    C=start_data.C;
    xst=[Ast Bst';Bst,-beta*C]\[fst;gst];
end
% compute residual of Stokes solution
nubeta=beta/viscosity;
N = navier_q1(xy,ev,xst);
Anst = viscosity*A + gamma2*RB{4}+ [N, sparse(nv,nv); sparse(nv,nv), N];
[Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,pb);
nlres = [Anst,Bst';Bst,-nubeta*C]*xst-[fst;gst];
nlres_norm  = norm(nlres);
fprintf('\n\ninitial nonlinear residual is %e ',nlres0_norm)
fprintf('\nStokes solution residual is %e\n', nlres_norm)
flowsol = xst;
%
%
pde=4;
it_p = 0;
%
% nonlinear iteration
%%% Picard startup step
while nlres_norm>nlres0_norm*tol_nl & it_p<maxit_p,
    it_p = it_p+1;
    fprintf('\nPicard iteration number %g \n',it_p),
    % compute Picard correction and update solution
    dxns = -[Anst,Bst';Bst,-nubeta*C]\nlres;
    xns = flowsol + dxns;
    % compute residual of new solution
    N = navier_q1(xy,ev,xns);
    Anst = viscosity*A + gamma2*RB{4}+ [N, sparse(nv,nv); sparse(nv,nv), N];
    [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,pb);
    nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
    nlres_norm = norm(nlres);
    nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
    fprintf('nonlinear residual is %e',nlres_norm)
    fprintf('\n   velocity change is %e\n',soldiff)
    flowsol = xns;
    % end of Picard iteration loop
end
%%%
%
it_nl = it_p;
it_n = 0;
%%%% Newton iteration loop
while (nlres_norm > nlres0_norm*tol_nl) & (it_nl < maxit_p + maxit_n),
    it_n = it_n+1;
    it_nl = it_nl+1;
    fprintf('\nNewton iteration number %g \n',it_n),
    % compute Jacobian of current solution
     [Nxx,Nxy,Nyx,Nyy] = newton_q1(xy,ev,flowsol);
    J = viscosity*A + gamma2*RB{4}+ [N + Nxx, Nxy; Nyx, N + Nyy];
    Jnst = newtonbc(J,xy,pb);
    % compute Newton correction and update solution
    dxns = -[Jnst,Bst';Bst,-nubeta*C]\nlres;
    xns = flowsol + dxns;
    % compute residual of new solution
     N = navier_q1(xy,ev,xns);
    Anst = viscosity*A + gamma2*RB{4}+ [N, sparse(nv,nv); sparse(nv,nv), N];
    [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,pb);
    nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
    nlres_norm = norm(nlres);
    nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
    fprintf('nonlinear residual is %e',nlres_norm)
    fprintf('\n   velocity change is %e\n',soldiff)
    flowsol = xns;
    %% end of Newton iteration loop
end
%
if nlres_norm <= nlres0_norm * tol_nl,
    fprintf('\nfinished, nonlinear convergence test satisfied\n\n');
else
    fprintf('\nfinished, stopped on iteration counts\n\n');
end
%% compute comunication data mu1
u2=flowsol(1:2*nvtx);p2=flowsol(2*nvtx+1:end);
u2bc=[u2(port);u2(port+nvtx)];
end