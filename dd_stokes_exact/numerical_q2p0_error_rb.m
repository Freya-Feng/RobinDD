function [uele_h1,uele_l2,pele_l2]=numerical_q2p0_error_rb(xy,ev,xst,handle_solvelocity,handle_p,handle_velocity_gradient)
% compute the exact error for q1p0 element with 10*10 points Gauss rule
% no obvious difference with using mathlab dblquad 
% input
%          xy               vertex coordinate vector  
%          ev               element mapping matrix
%          xst              q1p0 solution
%          handle_velocity  function handle of the gradient of exact veclocity
%          handle_p         function handle of the exact pressure
% output
%         uele_h1           element contribution for velocity error in H1 semi norm
%         uele_l2           element contribution for velocity error in L2 norm
%         pele_l2           element contribution for puressure error in L2 norm
% call function ana_test_q_grad
% QL; 12 Sep 2011

x=xy(:,1); y=xy(:,2);
nel=length(ev(:,1));
nv=size(xy,1);
uele_h1=zeros(nel,1);
uele_l2=zeros(nel,1);
pele_l2=zeros(nel,1);
%% set up 1D 10 Gauss points
ngpt=10;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,wt] = gausspoints_twod(oneg,onew);
nngpt=ngpt^2;     

% inner loop over elements    
        for ivtx = 1:9
          xl_v(:,ivtx) = x(ev(:,ivtx));
          yl_v(:,ivtx) = y(ev(:,ivtx)); 
          u_v(:,ivtx)=xst(ev(:,ivtx));
          v_v(:,ivtx)=xst(ev(:,ivtx)+nv);
        end
        ph=xst(1+2*nv:end);
 
    for igpt = 1:nngpt
         sigpt=s(igpt);
         tigpt=t(igpt);
         wght=wt(igpt);
% evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
         [psi,dpsidx,dpsidy] = qderiv(sigpt,tigpt,xl_v,yl_v);
	     duhdx=zeros(nel,1);
         duhdy=zeros(nel,1);
	     dvhdx=zeros(nel,1);
         dvhdy=zeros(nel,1);
         uh=zeros(nel,1);
         vh=zeros(nel,1);
         for i=1:9
             duhdx=duhdx+u_v(:,i).*dpsidx(:,i).*invjac;
             duhdy=duhdy+u_v(:,i).*dpsidy(:,i).*invjac;
             dvhdx=dvhdx+v_v(:,i).*dpsidx(:,i).*invjac;
             dvhdy=dvhdy+v_v(:,i).*dpsidy(:,i).*invjac;
             uh=uh+u_v(:,i).*psi(:,i);
             vh=vh+v_v(:,i).*psi(:,i);
         end
         [xigpt,yigpt]=qonemap(sigpt,tigpt,xy,ev);   
         [dudx,dudy,dvdx,dvdy]=handle_velocity_gradient(xigpt,yigpt);
         p=handle_p(xigpt,yigpt);
         [u,v]=handle_solvelocity(xigpt,yigpt);
         uele_h1(:)=uele_h1(:)+wght*(dudx-duhdx).^2.*jac+wght*(dudy-duhdy).^2.*jac...
                              +wght*(dvdx-dvhdx).^2.*jac+wght*(dvdy-dvhdy).^2.*jac;
         uele_l2(:)=uele_l2(:)+wght*(u-uh).^2.*jac+wght*(v-vh).^2.*jac;
         pele_l2(:)=pele_l2(:)+wght*(p-ph).^2.*jac;
% end of Gauss point loop
   end
%