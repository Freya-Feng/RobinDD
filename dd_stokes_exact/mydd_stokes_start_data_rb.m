function start_data=mydd_stokes_start_data_rb(LR)
% generate the Galerkin system
% LR=1(for left domain)/LR=2(for right domain)
%Copyright (c) 2023 Yani Feng, Qifeng Liao
channel_domain_rb,
load channel_grid
mx=(max(xy(:,1))+min(xy(:,1)))/2;
if LR==1 % left domain
   [mv,xy] = mesh_minus(mv,xy,[mx,1,0,1]);
elseif LR==2 % right domain
   [mv,xy] = mesh_minus(mv,xy,[0,mx,0,1]);
end
x=unique(sort(xy(:,1)));
y=unique(sort(xy(:,2)));
% bound number vary with xy, but mbound doesn't change.
[mbound,bound]=find_boundary_elements(mv);
% port number
tol=10^(-10);
t=find(abs(xy(bound,1)-mx)<tol);
port=bound(t);
q_in=2;
%q_in=2;
qmethod=q_in-1;
if qmethod==1
    [ev,ee,ebound,xyp] = q1p0grid(x,y,xy,mv,bound,mbound);
    [A,B,Q,C,G,Bx,By,f,g,fn,RB] = mystokes_q1p0(xy,xyp,mv,ev,ebound);
    start_data.C=C;
    start_data.ev=ev;
    start_data.mbound=mbound;
    start_data.ebound=ebound;
end
%%
% construct start_data
start_data.xy=xy;
start_data.xyp=xyp;
start_data.x=x;
start_data.y=y;
start_data.bound=bound;
start_data.A=A;
start_data.B=B;
start_data.M=Q;
start_data.B=B;
start_data.G=G;
start_data.Bx=Bx;
start_data.By=By;
start_data.f=f;
start_data.g=g;
start_data.fn=fn;
start_data.RB=RB;
start_data.xyport=xy(port,:);
start_data.port=port;
end