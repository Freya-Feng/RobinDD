function start_data=feng_dd_cavity_start_data(LR)
% generate the Galerkin system
% LR=1(for left domain)/LR=2(for right domain)
%Copyright (c) 2023 Yani Feng, Qifeng Liao
cavity_domain
load cavity_grid.mat
mx=(max(xy(:,1))+min(xy(:,1)))/2;
if LR==1 % left domain
   [mv,xy] = mesh_minus(mv,xy,[mx,1,-1,1]);
elseif LR==2 % right domain
   [mv,xy] = mesh_minus(mv,xy,[-1,mx,-1,1]);
end
x=unique(sort(xy(:,1)));
y=unique(sort(xy(:,2)));
% bound number vary with xy, but mbound doesn't change.
[mbound,bound]=find_boundary_elements(mv);
% port number
tol=10^(-10);
t=find(abs(xy(bound,1)-mx)<tol);
port=bound(t);
% % construct the new bound by deleting the interface
% bound=setdiff(bound,port);
% % replace port edge (edge 2) by edge nan
% for i=1:length(mbound(:,1))
%     ref=mbound(i,1); ref2=mbound(i,2);
%     for k=1:length(port)
%         if(any(mv(ref,:)==port(k)) & mod(ref2,2)==0) 
%            mbound(i,2)=nan; 
%         end
%     end
% end
% % delete boudary macroelement and then formulate mbound
% t=find(mbound(:,2)>0);
% mbound=mbound(t,:);
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
% elseif qmethod==8
%     [x,y,xy,xyp] = q2p1grid(x,y,xy,mv,bound);
%     [A,BB,QQ,G,Bx,By,f,gg,fn,RB] = mystokes_q2p1(xy,xyp,mv,mbound);
%      nnp=length(gg); ppk=[1:3:nnp];
%      B=BB(ppk,:); Q=QQ(ppk,ppk); g=gg(ppk);
%     %B=BB;Q=QQ;g=gg;
%     fprintf('Reduction to Q2-P0 Stokes system done.\n')
%     start_data.ev=mv;
%     start_data.ebound=mbound;
end
%%
% construct start_data
start_data.xy=xy;
start_data.xyp=xyp;
start_data.x=x;
start_data.y=y;
%start_data.mv=mv;
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
% start_data.Mp=[];
% start_data.fn=fn;
% start_data.fr=fr;
% start_data.qmethod=1;
% start_data.nu=length(fr);
% start_data.np=0;
% start_data.dim_par=dim_par;
% start_data.Nco=Nco;
% start_data.if_plot=if_plot;
%start_data.square_size=square_size;
start_data.xyport=xy(port,:);
start_data.port=port;
% start_data.port_p=port_p;
% start_data.xyport_p=xyp(port_p,:);
% start_data.Q1pod=Q1pod;
% start_data.Q2pod=Q2pod;
end