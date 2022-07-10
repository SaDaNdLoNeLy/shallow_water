% Author: Karthik Velakur
% Date: July 2014
% This code was written for a numerical methods course (AMATH 745) at... 
% University of Waterloo. 
% Refer the PDF file in same archive for a mathematical background. 

% FVM for shallow water equations
clear all; close all;
% setup grid
m = 60; dx = 2/m; dy = 2/m;
x = -1-dx:dx:1+dx; y = -1-dy:dy:1+dy; [xx,yy] = meshgrid(x,y);
g = 1; c = 0.5;

h = ones(size(xx)); h(xx>=-0.5 & xx<=0.5 & yy>=-0.5 & yy<=0.5) = 2;
% U is our unknown matrix. U(:,:,1)=h,U(:,:2)=hu,U(:,:,3)=hv
U = zeros([size(h) 3]); U(:,:,1) = h;
u = zeros(size(xx)); v = u;

% vectors for circular shift
shiftp1 = circshift((1:length(x))',1);
shiftm1 = circshift((1:length(x))',-1);

mesh(x,y,U(:,:,1)), colormap gray, axis([-1 1 -1 1 0.5 2.5])
title('hit enter to continue')
xlabel x, ylabel y; zlabel h; pause;

t = 0; dt = 0; tstop = 3.0; ii = 1;
numplots = 3; tplot = [1.35;3.0];
Uplot = zeros([size(U) length(tplot)+1]); Uplot(:,:,:,1) = U;
styles = {'k:','k--','k-'};
while t < tstop
    Uold = U; uold = u; vold = v; told = t; t = t + dt;
    
    % calculate lambda = |u| + sqrt(gh) used for finding flux
    lambdau = 0.5*abs(uold+uold(:,shiftm1)) +...
        sqrt(g*0.5*(Uold(:,:,1)+Uold(:,shiftm1,1)));
    lambdav = 0.5*abs(vold+vold(shiftm1,:)) +...
        sqrt(g*0.5*(Uold(:,:,1)+Uold(shiftm1,:,1)));
    lambdamax = norm([lambdau(:); lambdav(:)],Inf);
    
    dt = c*(dx/lambdamax);
    % adjust dt to produce plots at the right time
    if (ii<=length(tplot) && tplot(ii)>=told && tplot(ii)<=t+dt)
        dt = tplot(ii)-t;
        ii = ii + 1;
    end
    
    huv = Uold(:,:,2).*Uold(:,:,3)./Uold(:,:,1);
    ghh = 0.5*g*Uold(:,:,1).^2;
    % calculate (hu,hu^2+gh^2/2,huv)
    lffu = cat(3,Uold(:,:,2),Uold(:,:,2).^2./Uold(:,:,1)+ghh,huv);
    % calcualte (hv,huv,hv^2+gh^2/2)
    lffv = cat(3,Uold(:,:,3),huv,Uold(:,:,3).^2./Uold(:,:,1)+ghh);
    % calculate fluxes
    fluxx =  0.5*(lffu+lffu(:,shiftm1,:)) - ...
        0.5*bsxfun(@times,Uold(:,shiftm1,:)-Uold,lambdau);
    fluxy =  0.5*(lffv+lffv(shiftm1,:,:)) - ...
        0.5*bsxfun(@times,Uold(shiftm1,:,:)-Uold,lambdav);
    % time step
    U = Uold - (dt/dx)*(fluxx - fluxx(:,shiftp1,:)) ...
        - (dt/dy)*(fluxy - fluxy(shiftp1,:,:));
    
    % impose boundary conditions on h
    U(1:end,end,1) =  U(1:end,end-1,1); U(1:end,1,1) =  U(1:end,2,1);
    U(end,1:end,1) =  U(end-1,1:end,1); U(1,1:end,1) =  U(2,1:end,1);
    % on hu
    U(1:end,end,2) = -U(1:end,end-1,2); U(1:end,1,2) = -U(1:end,2,2);
    U(end,1:end,2) =  U(end-1,1:end,2); U(1,1:end,2) =  U(2,1:end,2);
    % on hv
    U(1:end,end,3) =  U(1:end,end-1,3); U(1:end,1,3) =  U(1:end,2,3);
    U(end,1:end,3) = -U(end-1,1:end,3); U(1,1:end,3) = -U(2,1:end,3);
    
    % u = hu./h;            % v = hv./h;
    u = U(:,:,2)./U(:,:,1);   v = U(:,:,3)./U(:,:,1);
    
    % display a movie of height
    mesh(x,y,U(:,:,1)), colormap gray, axis([-1 1 -1 1 0 2.5])
    title(['t = ' num2str(t+dt)])
    xlabel x, ylabel y, zlabel y, pause(0.001)
    %if (ismember(t+dt,tplot))
    if (any(tplot-t-dt==0))
        Uplot(:,:,:,ii) = U;      % store U for plotting
    end
end

% % plotting
% figure; tplot = [0;tplot];
% for ii = 1:length(tplot)
%     subplot(3,3,ii), mesh(x,y,Uplot(:,:,1,ii)), axis equal
%     colormap gray, axis([-1 1 -1 1 0.5 2.5]), xlabel x, ylabel y, zlabel h
%     title(['Height, t = ', num2str(tplot(ii))])
%     subplot(3,3,3+ii), mesh(x,y,Uplot(:,:,2,ii)./Uplot(:,:,1,ii))
%     colormap gray, axis([-1 1 -1 1 -1 1]), xlabel x, ylabel y, zlabel u
%     title(['Vel. u, t = ', num2str(tplot(ii))])
%     subplot(3,3,6+ii), mesh(x,y,Uplot(:,:,3,ii)./Uplot(:,:,1,ii))
%     colormap gray, axis([-1 1 -1 1 -1 1]), xlabel x, ylabel y, zlabel v
%     title(['Vel. v, t = ', num2str(tplot(ii))])
% end
