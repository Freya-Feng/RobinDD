function start_data=myldc_stokes_start_data
channel_domain_rb,
load channel_grid
x=unique(sort(xy(:,1)));
y=unique(sort(xy(:,2)));
% bound number vary with xy, but mbound doesn't change.
[mbound,bound]=find_boundary_elements(mv);
[ev,ee,ebound,xyp] = q1p0grid(x,y,xy,mv,bound,mbound);
xi=eye(2);
global global_domaininfo

xi=eye(2);
%% matrix generating 
for i=1:2
    global_domaininfo.la=xi(i,1);
    global_domaininfo.ra=xi(i,2);
    [K{i},B,Q,fn,fr] = myfemq1_poisson(xy,ev,ebound);
end
%% matrix generating
[A,B,Q,C,G,Bx,By,f,g,fn,RB] = mystokes_q1p0(xy,xyp,mv,ev,ebound);
start_data.C=C;
start_data.ev=ev;
start_data.mbound=mbound;
start_data.ebound=ebound;

end