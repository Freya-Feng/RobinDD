function feng_solplot_navier_stokes(sol,xy,x,y,obs,bndxy,fig)
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
if size(obs,1)~=0
   [II,JJ] = findobsXY(obs,X,Y,bndxy); xysol(II,JJ)=nan;
end
figure(fig)
mesh(X,Y,xysol),axis('square')
view(330,30)
end