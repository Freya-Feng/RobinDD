function [xigpt,yigpt]=qonemap(s,t,xy,ev)
%qonemap maps refernece a point in [-1,1]^2 to phsical elements
%[xigpt,yigpt]=qonemap(s,t,xy,ev);
%input
% s       x coordinate of a point in reference element [-1,1]^2
% t       y coordinate of a point in reference element
% xy      Q2 nodal coordinate vector 
% ev      element vertex matrix
%output
% xight   x coordinates in physical elements
% yight   y coordinates in physical elements
%   IFISS function: QL; 31 Aug 2011.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage


nel=size(ev,1);

x=xy(:,1);y=xy(:,2);
% inner loop over elements    
for ivtx = 1:4
   xl_v(:,ivtx) = x(ev(:,ivtx));
   yl_v(:,ivtx) = y(ev(:,ivtx)); 
end
xigpt=xl_v(:,1).*(s-1).*(t-1)/4-xl_v(:,2).*(s+1).*(t-1)/4 ...
        +xl_v(:,3).*(s+1).*(t+1)/4-xl_v(:,4).*(s-1).*(t+1)/4;
yigpt=yl_v(:,1).*(s-1).*(t-1)/4-yl_v(:,2).*(s+1).*(t-1)/4 ...
        +yl_v(:,3).*(s+1).*(t+1)/4-yl_v(:,4).*(s-1).*(t+1)/4;

    
return

