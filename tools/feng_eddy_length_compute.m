function eddy_r=feng_eddy_length_compute(u_gal,xy)
% compute the length of the lower eddy and the upper eddy
% input 
%         u_gal  velocity solution
%         xy     nodal coordinate vector
% output
%         eddy_r [lower eddy length, the start point of upper eddy, the end point of upper eddy, upper eddy length]
%   IFISS function: YF; 1 May 2023.
% Copyright (c) 2023 Yani Feng, Qifeng Liao
nvtx=length(xy(:,1));x=unique(xy(:,1));y=unique(xy(:,2));
u0=u_gal(1:nvtx);
u_up=reshape(u0(find(xy(:,2)<1&xy(:,2)>0)),length(find(y<1&y>0)),length(x));
u_up_min=min(u_up);
x_n=x(find(x>=0));
u_low=reshape(u0(find(xy(:,2)<0&xy(:,2)>-1)),length(find(y<0&y>-1)),length(x_n));
u_low_min=min(u_low);
%% upper eddy
x_up=x(find(u_up_min<0));
up_start=min(x_up);up_end=max(x_up);
upper_l=up_end-up_start;
%% lower eddy
x_low=x_n(find(u_low_min<0));
lower_l=max(x_low);
%%
eddy_r=[lower_l,up_start,up_end,upper_l];
end