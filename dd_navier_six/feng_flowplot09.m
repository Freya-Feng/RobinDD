function feng_flowplot09(qmethod,sol,start_data,contourn,spc,fig,fign)
%FLOWPLOT09 plots flow data on general domain
%   flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,1,69)
%   input
%          qmethod    mixed method 
%          sol        flow solution vector
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix    
%          A          vector diffusion matrix
%          xy         velocity nodal coordinate vector  
%          xyp        pressure nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bound      boundary vertex vector
%          spc        uniform/nonuniform streamline switch 
%          fig        figure number
%
% calls function streambc.m to set boundary values
%   IFISS function: HCE; DJS; 20 September 2016.
%   Adapted from code written by M. Wu, 2009
% Copyright (c) 2009 D.J. Silvester, H.C. Elman, A. Ramage;
%
A=start_data.A;
Bx=start_data.Bx;
By=start_data.By;
x=start_data.x;
y=start_data.y;
xy=start_data.xy;
xyp=start_data.xyp;
bound=start_data.bound;
obs=start_data.obs;
bndxy=start_data.bndxy;
bnde=start_data.bnde;

if ishandle(fig), clf(fig); end
nvtx=length(xy); nu=2*nvtx;   %np=length(xyp);
Asv=A(1:nvtx,1:nvtx);
%
%% compute auxiliary quantities
u=sol(1:nu);p=sol(nu+1:end);
f=[By,-Bx]*u;
[Asv,fsv] = streambc(Asv,f,xy,bound);
phi=Asv\fsv;
%
%% plot pressure
if qmethod==2
   xx=x(1:2:end); yy=y(1:2:end);
elseif qmethod==3
   p=p(1:3:end); xx=x(1:end); yy=y(1:end);
else
   xx=x(1:end); yy=y(1:end);
end
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(xx,yy);
xysol = griddata(xyp(:,1),xyp(:,2),p,X,Y);
if size(obs,1)~=0
   [II,JJ] = findobsXY(obs,X,Y,bndxy);  xysol(II,JJ)=nan;
end
figure(fig)
colormap jet
subplot(212), mesh(X,Y,xysol),  axis('tight'), %colorbar('EastOutside')
view(10,20)
title('Pressure field','FontSize',12)
%
%% plot velocity
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);
if size(obs,1)~=0
   [II,JJ] = findobsXY(obs,X,Y,bndxy); xysol(II,JJ)=nan;
end
if spc==1
   subplot(211),contour(X,Y,xysol,contourn),axis('tight'), %colorbar('EastOutside')
   title('Streamlines: uniform','FontSize',12); 
elseif spc == 2
   cn = fix(contourn/2-1);
   ch = 21/cn;
   v=(-15:ch:6)';v=exp(v); list=sort([-v;v]);
   subplot(211),contour(X,Y,xysol,list),axis('tight'),  %colorbar('EastOutside')
   title('Streamlines: selected','FontSize',12); 
end
axis('off')
%
% plot the boundary
hold on
for i = 1:size(bnde,1)
   plot([bndxy(bnde(i,1),1), bndxy(bnde(i,2),1)],[bndxy(bnde(i,1),2),bndxy(bnde(i,2),2)],'-k')
end
hold off;
%
%% plot velocity by using the  format of unsteady_navier

figure(fign)
%set(gcf,'Position',[5,5,600,600]);
set(gcf,'Position',[10,10,1200,1200]);
[X,Y]=meshgrid(x,y);
ax = [min(x)-1 max(x)+1 min(y)-1 max(y)+1];
xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);
maxphi=max(max(xysol)); minphi=min(min(xysol));
vneg=[minphi:-minphi/6:0];
vpos=[maxphi/20:maxphi/20:19*maxphi/20];
vpospos=[79*maxphi/80: maxphi/320:maxphi];
contour(X,Y,xysol,[vneg,vpos,vpospos])
title('Streamlines: uniform');
axis equal, axx=ax(1:4); axx(2)=min(30,axx(2));
axis(axx(1:4));
stepx;
return
