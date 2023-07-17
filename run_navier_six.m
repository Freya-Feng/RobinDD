%%Test problem 4: the Navier-Stokes problem-flow over a backward facing step with six subdomains in section 6.3
% L_D=30
%delete(gcp('nocreate'))
clear all
close all
gohome
system('/bin/cp ./stokes_flow/test_problems/backwardstep_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/backwardstep_bc.m ./stokes_flow/stream_bc.m');
%%
N=220;
num=6; %split long domain into 6 parts
%%
% generate start_data
parpool(8)
parfor i=1:num
    start_data_domain{i}=feng_dd_navier_start_data(i,num);
end
%%
start_data=feng_dd_navier_start_data(num+1,num);%full domain
%%
xyport=start_data_domain{1}.xyport;
nport=length(xyport);
%%
%viscosity=default('viscosity parameter (default 1/50)',1/50);
viscosity=1/100;
nlmethod=3;
%nlmethod=default('Picard/Newton/hybrid linearization 1/2/3 (default hybrid)',3);
%tol_nl=default('nonlinear tolerance (default 1.1*eps)',1.1*eps);
%tol_nl=1.1*eps;
tol_nl=1.0e-12;
%% start iterations
% initial guess for g1
parfor i=1:num-1
    g1{i}=zeros(2*nport,1);
end
g2=g1;
%%set parameters
theta1=1/32;theta2=1/32;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
gamma1=1/2;gamma2=0.5/viscosity;
%%
e=[]; 
flag=1;
n_iter=0;
while flag
    n_iter=n_iter+1;
    g1_old=g1;g2_old=g2;
    fprintf('\nRobin-Robin iteration number %g \n',n_iter)
    % left part
    [u{1},p{1},ubc1]=feng_navier_dd_left_rb(start_data_domain{1},g1{1},gamma1,viscosity,nlmethod,tol_nl);
    % middle part
    parfor j=2:num-1
        [u{j},p{j},ubc_l{j},ubc_r{j}]=feng_navier_dd_middle_rb(start_data_domain{j},g2{j-1},g1{j},gamma1,gamma2,viscosity,nlmethod,tol_nl);
    end
    % right part
    [u{num},p{num},ubc_f]=feng_navier_dd_right_rb(start_data_domain{num},g2{num-1},gamma2,viscosity,nlmethod,tol_nl);
    %parallel upadat left g1{1},g2{1}
    g1{1}=theta1*g1{1}+(1-theta1)*((gamma1+gamma2)*ubc_l{2}-g2{1});
    g2{1}=theta2*g2{1}+(1-theta2)*((gamma1+gamma2)*ubc1-g1{1});
    %parallel update middle g1,g2
    for j=2:num-2
        g1{j}=theta1*g1{j}+(1-theta1)*((gamma1+gamma2)*ubc_l{j+1}-g2{j});
        g2{j}=theta2*g2{j}+(1-theta2)*((gamma1+gamma2)*ubc_r{j}-g1{j});
    end
    % parallel update right g1{num-1},g2{num-1}
    g1{num-1}=theta1*g1{num-1}+(1-theta1)*((gamma1+gamma2)*ubc_f-g2{num-1});
    g2{num-1}=theta2*g2{num-1}+(1-theta2)*((gamma1+gamma2)*ubc_r{num-1}-g1{num-1});
    sum=[abs(ubc1-ubc_l{2});abs(ubc_r{num-1}-ubc_f)];
    for j=2:num-2
        sum=[sum;abs(ubc_r{j}-ubc_l{j+1})];
    end
    if max(sum)<1.0e-8
        flag=0;
    end
end
%% full solution right!!
[u_real,p_real]=feng_navier_full_rb(start_data,viscosity,nlmethod,tol_nl);
%% compare
x_gal=[];y_gal=x_gal;xy=[];xyp=[];p_gal=[];
for i=1:num
    u0=u{i};x_gal=[x_gal;u0(1:length(u0)/2)];y_gal=[y_gal;u0(1+length(u0)/2:end)];p_gal=[p_gal;p{i}];
    xy=[xy;start_data_domain{i}.xy];xyp=[xyp;start_data_domain{i}.xyp];
end
[xy,t]=unique(xy,'rows');x_gal=x_gal(t);y_gal=y_gal(t);
[xyp,t1]=unique(xyp,'rows');p_gal=p_gal(t1);
x=start_data.x;y=start_data.y;
u_gal=[x_gal;y_gal];
%%
[e_compare,ur]=sc_fem_stokes_compare(u_real,start_data.xy,u_gal,xy);
[p_compare,up]=sc_fem_compare(p_real,start_data.xyp,p_gal,xyp);
%% plot solution
%spc=default('uniform/exponential streamlines 1/2 (default uniform)',1);
%contourn = default('number of contour lines (default 50)',50);
spc=1;contourn=50;
qmethod=1;
%%
xst=[ur;up];
fprintf(' Total iterations: %8.3e \n',n_iter)
fprintf(' relative error (u_dd-u_h)/u: %8.3e \n',sqrt((u_real-ur)'*start_data.G*(u_real-ur))/sqrt(u_real'*start_data.G*u_real))
fprintf(' relative error (p_dd-p_h)/u: %8.3e \n',sqrt((p_real-up)'*start_data.M*(p_real-up))/sqrt(p_real'*start_data.M*p_real))
%%
backflowplot09(qmethod,[u_real;p_real],start_data,contourn,spc,11,12)
backflowplot09(qmethod,[ur;up],start_data,contourn,spc,21,22)
%% eddy length
fprintf(' Eddy length for global solutions: %8.3e \n',feng_eddy_length_compute(u_real(1:length(u_real)/2),start_data.xy))
fprintf(' Eddy length for domain decomposition solutions: %8.3e \n',feng_eddy_length_compute(ur(1:length(ur)/2),start_data.xy))
