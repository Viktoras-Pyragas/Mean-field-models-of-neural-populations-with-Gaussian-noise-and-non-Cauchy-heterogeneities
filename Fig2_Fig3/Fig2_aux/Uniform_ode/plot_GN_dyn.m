function plot_GN_dyn
%----------------------------------------------------------------------
% Here we plot the graphs from the output from truncated equations for
% Fourier components, with flat heterogeneity;
%----------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%Yz1(nq,1:4)=[t(1:end-1,1) Rv(1:end-1,1) Vv(1:end-1,1) ym(1:end-1,Ns+1)];
% Rv -> averaged spiking rate;
% Vv -> averaged potential;
% ym(:,Ns+1) -> solution of the relaxation equation
load('GM_dyn.mat','Yz1'); % time series of solution; M=300;
Yz1(:,1)=Yz1(:,1)-Yz1(1,1); % shifting the timne scale;
Yz1t=Yz1;
%-------------------------------------------------------------------------------------------
% Plotting the results:
%-------------------------------------------------------------------------------------------
figure
title('n=2, M=300, r(t)');
hold on
plot(Yz1t(:,1),Yz1t(:,2),'m-'); 
hold off

figure
title('n=2, M=300, v(t)');
hold on
plot(Yz1t(:,1),Yz1t(:,3),'m-'); 
hold off

figure
title('n=2, M=300, (r,v)');
hold on
plot(Yz1t(:,2),Yz1t(:,3),'b--'); 
hold off

end