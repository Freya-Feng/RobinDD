function [mbound,bound]=find_boundary_elements(mv)
%find_boundary_elements construct eboundt and bound for pure Dirichlet
%boundary condition
%[mbound,bound]=find_boundary_elements(mv);
%   input
%          mv             element mapping matrix
%   output
%          mbound         element edge boundary matrix 
%          bound          boundary nodes 
%
%%%QL; 3 Aug 2012.

%% Exterior boundary 
mbound=[];
bound=[];
nel=length(mv(:,1));
for i=1:nel
    %edge1
    t=find((mv(:,1)==mv(i,1)|mv(:,2)==mv(i,1)|mv(:,3)==mv(i,1)|mv(:,4)==mv(i,1))...
          &(mv(:,1)==mv(i,2)|mv(:,2)==mv(i,2)|mv(:,3)==mv(i,2)|mv(:,4)==mv(i,2)));
    if length(t)<2, mbound=[mbound;i,1]; bound=[bound,mv(i,1),mv(i,2),mv(i,5)]; end
    %edge2
    t=find((mv(:,1)==mv(i,2)|mv(:,2)==mv(i,2)|mv(:,3)==mv(i,2)|mv(:,4)==mv(i,2))...
          &(mv(:,1)==mv(i,3)|mv(:,2)==mv(i,3)|mv(:,3)==mv(i,3)|mv(:,4)==mv(i,3)));
    if length(t)<2, mbound=[mbound;i,2]; bound=[bound,mv(i,2),mv(i,3),mv(i,6)]; end
    %edge3
    t=find((mv(:,1)==mv(i,3)|mv(:,2)==mv(i,3)|mv(:,3)==mv(i,3)|mv(:,4)==mv(i,3))...
          &(mv(:,1)==mv(i,4)|mv(:,2)==mv(i,4)|mv(:,3)==mv(i,4)|mv(:,4)==mv(i,4)));
    if length(t)<2, mbound=[mbound;i,3]; bound=[bound,mv(i,3),mv(i,4),mv(i,7)]; end
    %edge4
    t=find((mv(:,1)==mv(i,4)|mv(:,2)==mv(i,4)|mv(:,3)==mv(i,4)|mv(:,4)==mv(i,4))...
          &(mv(:,1)==mv(i,1)|mv(:,2)==mv(i,1)|mv(:,3)==mv(i,1)|mv(:,4)==mv(i,1)));
    if length(t)<2, mbound=[mbound;i,4]; bound=[bound,mv(i,4),mv(i,1),mv(i,8)]; end
end

bound=unique(bound)';
end