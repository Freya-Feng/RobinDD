function ur=changexy(ucs,xycs,xy)
nvtx=length(ucs);
tol=10^(-10);
ur=0*ucs;
for i=1:nvtx
     % find static condenstion solution wich match this point
    t=min(find(abs(xycs(:,1)-xy(i,1))<tol&abs(xycs(:,2)-xy(i,2))<tol));
    ur(i)=ucs(t);
end
end