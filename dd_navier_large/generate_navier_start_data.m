function start_data=generate_navier_start_data(LR,num,L)
% generate the Galerkin system
% LR represents the index of subdomain
% num is the number of subdomains
% L is the length of global domain
%Copyright (c) 2023 Yani Feng, Qifeng Liao
load step_grid.mat
fdx=L/num;
if LR<num+1
    if LR==1
        [mv,xy] = mesh_minus(mv,xy,[fdx,L,-1,1]);
    else
        [mv,xy] = mesh_minus(mv,xy,[LR*fdx,L,-1,1]);
        [mv,xy] = mesh_minus(mv,xy,[-1,(LR-1)*fdx,-1,1]);
    end
end
x=unique(sort(xy(:,1)));
y=unique(sort(xy(:,2)));
% bound number vary with xy, but mbound doesn't change.
[mbound,bound]=find_boundary_elements(mv);
% port number
tol=10^(-10);
if LR==1
    t=find(abs(xy(bound,1)-fdx)<tol);
    port=bound(t);
elseif LR==num
    t=find(abs(xy(bound,1)-(num-1)*fdx)<tol);
    port=bound(t);
else
    t1=find(abs(xy(bound,1)-(LR-1)*fdx)<tol);t2=find(abs(xy(bound,1)-LR*fdx)<tol);
    port=unique([bound(t1);bound(t2)]);
end
qmethod=1;
if qmethod==1
    [ev,ee,ebound,xyp] = q1p0grid(x,y,xy,mv,bound,mbound);
    [A,B,Q,C,G,Bx,By,f,g,fn,RB] = mystokes_q1p0(xy,xyp,mv,ev,ebound);
    start_data.C=C;
    start_data.ev=ev;
    start_data.mv=mv;
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
start_data.obs=obs;
start_data.bndxy=bndxy;
start_data.bnde=bnde;
% start_data.Mp=[];
% start_data.fn=fn;
% start_data.fr=fr;
% start_data.qmethod=1;
% start_data.nu=length(fr);
% start_data.np=0;
% start_data.dim_par=dim_par;
% start_data.Nco=Nco;
start_data.xyport=xy(port,:);
start_data.port=port;
end