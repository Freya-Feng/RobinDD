%% Test problem 3: the Navier-Stokes problem-flow over a backward facing step with two subdomains in section 6.2
% L_D=5
clear all
close all
gohome
system('/bin/cp ./stokes_flow/test_problems/backwardstep_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/backwardstep_bc.m ./stokes_flow/stream_bc.m');
%%
N=150; % the number of iterations by Dirichlet-Neumann method
% generate start_data
start_data_left =feng_navier_start_data_two(1);
start_data_right =feng_navier_start_data_two(2);
start_data=feng_navier_start_data_two(3);
% mesh size h
h=2/(length(start_data_left.y)-1);
nport=length(start_data_left.xyport);
% xy are arranged from bottom to up and from left to right 
% then number each point
xyport1=start_data_left.xyport; 
xyport2=start_data_right.xyport;
% port number
port1=start_data_left.port;
port2=start_data_right.port;
% mesh size h
h=2/(length(start_data.y)-1)
%gamma1=2;gamma2=64*h;
% theta=2/3; % suitably chosen relaxation parameter
theta1=1/3;
theta2=1/3;
%% initial guess for g1
g1=zeros(2*nport,1); g2=zeros(2*nport,1);
%viscosity=default('viscosity parameter (default 1/50)',1/50);
viscosity=1/100;
%nlmethod=default('Picard/Newton/hybrid linearization 1/2/3 (default hybrid)',3);
nlmethod=3;
%tol_nl=default('nonlinear tolerance (default 1.1*eps)',1.1*eps);
tol_nl=1.e-12;
%% start iterations
flag=1;i=0;
gamma1=1/2;gamma2=4*h/viscosity;
%gamma1=16*viscosity; gamma2=64*viscosity/h;
while flag
    i=i+1;
    g1_old=g1;
    g2_old=g2;
    % information on the interface.
    [u1,p1,u1bc]=feng_navier_dd_left_rb(start_data_left,g1,gamma1,viscosity,nlmethod,tol_nl);
    [u2,p2,u2bc]=feng_navier_dd_right_rb(start_data_right,g2,gamma2,viscosity,nlmethod,tol_nl);
    %parallel update g1,g2
    g1=theta1*g1+(1-theta1)*((gamma1+gamma2)*u2bc-g2);
    g2=theta2*g2+(1-theta2)*((gamma1+gamma2)*u1bc-g1);
    criter=max([abs(g1-g1_old);abs(g2-g2_old)]);
    % save data
    e(i)=norm(u1bc-u2bc);
    if criter<1.e-8
        flag=0;
    end
    percentage_count(i,N,'Processing...')
end
%% full solution right!!
[u_real,p_real]=feng_navier_full_rb(start_data,viscosity,nlmethod,tol_nl);
%% compare
u1x=u1(1:length(u1)/2);u1y=u1(length(u1)/2+1:end);u2x=u2(1:length(u2)/2);u2y=u2(length(u2)/2+1:end);
xy=[start_data_left.xy;start_data_right.xy];[xy,t]=unique(xy,'rows');
x=unique([start_data_left.x;start_data_right.x]);y=unique([start_data_left.y;start_data_right.y]);
x_gal=[u1x;u2x];x_gal=x_gal(t);
y_gal=[u1y;u2y];y_gal=y_gal(t);
xyp=[start_data_left.xyp;start_data_right.xyp];[xyp,t1]=unique(xyp,'rows');
p_gal=[p1;p2];p_gal=p_gal(t1);
u_gal=[x_gal;y_gal];
%%
[e_compare,ur]=sc_fem_stokes_compare(u_real,start_data.xy,u_gal,xy);
[p_compare,up]=sc_fem_compare(p_real,start_data.xyp,p_gal,xyp);
%% plot solution
%spc=default('uniform/exponential streamlines 1/2 (default uniform)',1);
%contourn = default('number of contour lines (default 50)',50);
spc=1; contourn=50;
xst=[ur;up];
qmethod=1;
%%
fprintf(' Mesh size: %8.3e \n',h)
fprintf(' Total iterations: %8.3e \n',i)
% plot the interface error
figure(30)
semilogy(1:i,e,'-.r')
fprintf(' relative error (u_dd-u_h)/u: %8.3e \n',sqrt((u_real-ur)'*start_data.G*(u_real-ur))/sqrt(u_real'*start_data.G*u_real))
fprintf(' relative error (p_dd-p_h)/u: %8.3e \n',sqrt((p_real-up)'*start_data.M*(p_real-up))/sqrt(p_real'*start_data.M*p_real))
%%plot
backflowplot09(qmethod,xst,start_data,contourn,spc,11,12)
backflowplot09(qmethod,[u_real;p_real],start_data,contourn,spc,21,22)