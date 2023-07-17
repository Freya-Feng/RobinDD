function [nmv,nxy] = mesh_minus(mv,xy,hole_size)
%mesh_minus digs holes for given mesh 
% [nmv,nxy] = mesh_minus(mv,xy,hole_size)
% input
%  mv1 xy1 as in IFISS
% hole_size=[xmin,xmax,ymin,ymax]
%QL; 25 Sep 2012

xmin=hole_size(1); xmax=hole_size(2); ymin=hole_size(3); ymax=hole_size(4);

mel=length(mv(:,1));

tol=1d-10;

kt=[];
nxy=[];
for i=1:mel
    nodes=mv(i,:);
    coords=xy(nodes,:);
    xm=coords(9,1);
    ym=coords(9,2);
    % check if this element will be kept
    if xm<xmin | xm>xmax | ym<ymin | ym>ymax
        kt=[kt,i];
        nxy=[nxy;coords];
    end 
end
nxy=unique(nxy,'rows');

for i=1:length(kt)
    k=kt(i);
    nodes=mv(k,:);
    coords=xy(nodes,:);
    for j=1:9
        rx=coords(j,1);
        ry=coords(j,2);
        nmv(i,j)=find(abs(nxy(:,1)-rx)<tol & abs(nxy(:,2)-ry)<tol); 
    end
end


return
