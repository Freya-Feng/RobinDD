%% Test problem 4: the Navier-Stokes problem-flow over a backward facing step with 15/33/55/165/330 subdomains in section 6.3
% L_D=165
% set h=1/32, then the total degrees of freedom are 1027586
tic
clear all
close all
gohome
system('/bin/cp ./stokes_flow/test_problems/backwardstep_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/backwardstep_bc.m ./stokes_flow/stream_bc.m');
%%
% generate start_data
parpool('local',32)
N=100;
num=165*2; %split long domain into num parts
outbnd=165;
newstep_domain;%1/32
L=outbnd;
%% 
parfor i=1:num
    start_data_domain{i}=generate_navier_start_data(i,num,L);
end
%%
h=2/(length(start_data_domain{2}.y)-1)
xyport=start_data_domain{1}.xyport;
nport=length(xyport);
%%
viscosity=1/200;
nlmethod=3;
tol_nl=1.0e-12;
%% start iterations
% initial guess for g1
parfor i=1:num-1
    g1{i}=zeros(2*nport,1);
end
g2=g1;
%%set parameters
theta1=1/32;theta2=1/32;  
%gamma1=1/4;gamma2=0.5/viscosity;
gamma1=4/6;gamma2=16/viscosity;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
%gamma1=64*64;gamma2=64*16/viscosity;
%%
e=[];
flag=1;
i=0;
while flag
%for i=1:N
    i=i+1;
    g1_old=g1;g2_old=g2;
    fprintf('\nRobin-Robin iteration number %g \n',i)
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
    % save data
    sum=abs(ubc1-ubc_l{2})+abs(ubc_r{num-1}-ubc_f);
    parfor j=2:num-2
        sum=sum+abs(ubc_r{j}-ubc_l{j+1});
    end
    e(i)=norm(sum);
    criter=[abs(ubc1-ubc_l{2});abs(ubc_r{num-1}-ubc_f)];
    for j=2:num-2
        criter=[criter;abs(ubc_r{j}-ubc_l{j+1})];
    end
    if max(criter)<1.e-8
        flag=0;
    end
end
toc
%% 
x_gal=[];y_gal=x_gal;xy=[];xyp=[];p_gal=[];
parfor i=1:num
    u0=u{i};x_gal=[x_gal;u0(1:length(u0)/2)];y_gal=[y_gal;u0(1+length(u0)/2:end)];p_gal=[p_gal;p{i}];
    xy=[xy;start_data_domain{i}.xy];xyp=[xyp;start_data_domain{i}.xyp];
end
[xy,t]=unique(xy,'rows');x_gal=x_gal(t);y_gal=y_gal(t);
[xyp,t1]=unique(xyp,'rows');p_gal=p_gal(t1);
%%
u_gal=[x_gal;y_gal];
%% compute eddy length
fprintf(' Total iterations: %8.3e \n',i)
fprintf(' Eddy length for domain decomposition solutions: %8.3e \n',feng_eddy_length_compute(u_gal(1:length(u_gal)/2),xy))



