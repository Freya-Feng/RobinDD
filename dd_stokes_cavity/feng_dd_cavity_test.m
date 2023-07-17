%% stokes problem-driven cavity
clear all
close all
gohome
% regularised lid driven cavity
system('/bin/cp ./stokes_flow/test_problems/regcavity_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/zero_bc.m ./stokes_flow/stream_bc.m');
%% Robin-Robin DD
N=400; % the number of iterations by Robin-Robin method
% generate start_data
start_data_left=feng_dd_cavity_start_data(1);
start_data_right=feng_dd_cavity_start_data(2);
start_data=feng_dd_cavity_start_data(3);
% mesh size h
h=2/(length(start_data.x)-1)
% parameters
viscosity=1;gamma1=16*viscosity;gamma2=viscosity*64/h;
%viscosity=1;gamma1=16;gamma2=64*2;
theta1=2/3;theta2=2/3;
%theta1=0;theta2=0;
nport=length(start_data_left.xyport);
g1=zeros(2*nport,1); g2=zeros(2*nport,1);
%% DD iterations
flag=1;i=0;
while flag
%for i=1:N
    i=i+1;
    g1_old=g1;g2_old=g2;
    [u1,p1,u1bc]=feng_cavity_dd_left(start_data_left,g1,gamma1,viscosity);
    [u2,p2,u2bc]=feng_cavity_dd_right(start_data_right,g2,gamma2,viscosity);  
    %parallel update g1,g2
    g1=theta1*g1+(1-theta1)*((gamma1+gamma2)*u2bc-g2);
    g2=theta2*g2+(1-theta2)*((gamma1+gamma2)*u1bc-g1);
    % stopping criteri
    criter(i)=max([abs(g1-g1_old);abs(g2-g2_old)]);
    e(i)=norm(u1bc-u2bc);
    if criter(i)<1.0e-8
        flag=0;
    end
    percentage_count(i,N,'Processing...')
%end
end
%%
fprintf(' Total iterations: %8.3e \n',i)
% plot the interface error
figure(30)
semilogy(1:i,e,'-.r')
%% combine u1 and u2
nvtx=length(u1)/2;u1x=u1(1:nvtx);u1y=u1(nvtx+1:end);u2x=u2(1:nvtx);u2y=u2(nvtx+1:end);
xy=[start_data_left.xy;start_data_right.xy];[xy,t]=unique(xy,'rows');
x=unique([start_data_left.x;start_data_right.x]);y=unique([start_data_left.y;start_data_right.y]);
x_gal=[u1x;u2x];x_gal=x_gal(t);
y_gal=[u1y;u2y];y_gal=y_gal(t);
p_gal=[p1;p2];xyp=[start_data_left.xyp;start_data_right.xyp];
[xyp,t]=unique(xyp,'rows');p_gal=p_gal(t);
xp=unique(xyp(:,1));yp=unique(xyp(:,2));
u_gal=[x_gal;y_gal];
%% full solution
[u,p]=feng_cavity_full_system(start_data,viscosity);
p_integral=h^2*sum(p);
p=p-p_integral/4;
%% compare
[e_compare,ur]=sc_fem_stokes_compare(u,start_data.xy,u_gal,xy);
up_integral=h^2*sum(p_gal);
p_gal=p_gal-up_integral/4;
[~,up]=sc_fem_compare(p,start_data.xyp,p_gal,xyp);
up1=up;
fprintf(' relative error (u_dd-u_h)/u: %8.3e \n',sqrt((u-ur)'*start_data.G*(u-ur))/sqrt(u'*start_data.G*u))
fprintf(' relative error (p_dd-p_h)/p_h: %8.3e \n',sqrt((p-up1)'*start_data.M*(p-up1))/sqrt(p'*start_data.M*p))
%%
spc=1;
qmethod=1;
xst=[ur;up1];
cavity_flowplot(qmethod,[u;p],start_data.By,start_data.Bx,start_data.A,start_data.xy,start_data.xyp,start_data.x,start_data.y,start_data.bound,spc,51);
cavity_flowplot(qmethod,xst,start_data.By,start_data.Bx,start_data.A,start_data.xy,start_data.xyp,start_data.x,start_data.y,start_data.bound,spc,61);
hold on
plot([0,0],[-1,1],'--r')

