%% Test problem 1: the Stokes problem with exact solutions in section 5.1
% split the unit square Omega=(0,1)^2 into two nonoverlapping subdomains
clear all
close all
system('/bin/cp ./domain_decomposition_Robin/dd_stokes_exact/myspecific_left_flow_rb.m ./diffusion/specific_bc.m');
%%
% generate start_data
start_data_left =mydd_stokes_start_data_rb(1);%left domain 
start_data_right =mydd_stokes_start_data_rb(2);%right domain 
start_data=mydd_stokes_start_data_rb(3);%the whole domain
nport=length(start_data_left.xyport);
port1=start_data_left.port;
port2=start_data_right.port;
h=1/(length(start_data.x)-1) % mesh size
viscosity=1/20;
gamma1=16*viscosity;gamma2=viscosity*64/h;
theta1=2/3;theta2=2/3;
%theta1=0;theta2=0;
%gamma1=2;gamma2=12; 
%%
% initial guess for g1
g1=zeros(2*nport,1); g2=zeros(2*nport,1);
flag=1;i=0;
while flag
    i=i+1;
    g1_old=g1;
    g2_old=g2;
    % information on the interface.
    [u1,p1,u1bc]=mystokes_dd_left_rb(start_data_left,g1,gamma1,viscosity);
    % right part
    [u2,p2,u2bc]=mystokes_dd_right_rb(start_data_right,g2,gamma2,viscosity);  
    %parallel update g1,g2
    g1=theta1*g1+(1-theta1)*((gamma1+gamma2)*u2bc-g2);
    g2=theta2*g2+(1-theta2)*((gamma1+gamma2)*u1bc-g1);
    criter=max([abs(g1-g1_old);abs(g2-g2_old)]);
    e(i)=norm(u1bc-u2bc);
    if criter<1.e-8
        flag=0;
    end
end
fprintf(' Total iterations: %8.3e \n',i)
% plot the interface error
figure(30)
semilogy(1:i,e,'-.r')
%%
nvtx=length(u1)/2;
u1x=u1(1:nvtx);u1y=u1(nvtx+1:end);
u2x=u2(1:nvtx);u2y=u2(nvtx+1:end);
xy=[start_data_left.xy;start_data_right.xy];
[xy,t]=unique(xy,'rows');
x=unique([start_data_left.x;start_data_right.x]);
y=unique([start_data_left.y;start_data_right.y]);
% DD result
x_gal=[u1x;u2x];x_gal=x_gal(t);
y_gal=[u1y;u2y];y_gal=y_gal(t); 
u_gal=[x_gal;y_gal];
p_gal=[p1;p2];
xyp=[start_data_left.xyp;start_data_right.xyp];
[xyp,t]=unique(xyp,'rows');p_gal=p_gal(t);
xp=unique(xyp(:,1));
yp=unique(xyp(:,2));

[u,p]=mystokes_full_system_rb(start_data,viscosity);%full_system consist of two system
% zero mean
p_integral=sum(p)*h^2;
p=p-p_integral/1;
p_gal=p_gal-sum(p_gal)*h^2/1;
%%
[e_compare,ur]=sc_fem_stokes_compare(u,start_data.xy,u_gal,xy); % change u_gal into ur
[p_compare,pr]=sc_fem_compare(p,start_data.xyp,p_gal,xyp);
fprintf(' relative error (u_dd-u_h)/u: %8.3e \n',sqrt((u-ur)'*start_data.G*(u-ur))/sqrt(u'*start_data.G*u))
fprintf(' relative error (p_dd-p_h)/u: %8.3e \n',sqrt((p-pr)'*start_data.M*(p-pr))/sqrt(p'*start_data.M*p))
xy_full=start_data.xy;
x1=unique(xy_full(:,1));x2=unique(xy_full(:,2));

%% compare true value and numerical value--full solution
handle_p=@test_pressure;
handle_solvelocity=@test_velocity;
handle_velocity_gradient=@test_velocity_gradient;
%% numerical solve full system
ev=start_data.ev;
xst=[ur;pr];
u_real=handle_solvelocity(start_data.x,start_data.y)/viscosity;
xp=unique(xyp(:,1));yp=unique(xyp(:,2));
p_real=handle_p(xp,yp);
[uele_h1,uele_l2,pele_l2]=numerical_q1p0_error_rb(start_data.xy,ev,xst,handle_solvelocity,handle_p,handle_velocity_gradient,viscosity);

q1p0_l2=sqrt(sum(pele_l2));
q1p0_semi_h1=sqrt(sum(uele_h1));
fprintf(' semi-H1 norm of u-u_h: %8.3e \n',q1p0_semi_h1)
fprintf(' L2 norm of p-p_h: %8.3e \n',q1p0_l2)

