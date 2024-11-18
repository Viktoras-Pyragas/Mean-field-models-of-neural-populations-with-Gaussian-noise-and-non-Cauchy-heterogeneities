%------------------------------------------------------------------------
%------------------------------------------------------------------------
% microscopic model; 
% interaction among neurons through delta-pulses;
% interaction inhibitory, synaptic equation involved;
% heterogeneity eta distributed as a uniform function;
% here the uniform distribution is generated deterministically;
%------------------------------------------------------------------------
close all;
clear;

format long;

N=50000; % number of neurons;
dt=0.001; % time step
Tint=60; % total time of integration;
Nt=round(Tint/dt); % total number of time steps;
J=2;  % coupling strength;

sig2=0.04;  % Gaussian noise amplitude (original)
sig=sqrt(2*sig2);
del0=0.1; % half-width (corresponds to del=d) (original)
d=del0;



muv=1; % average of parameter of heterogeneity eta; 
tau=1; % membane time
taus=1; % synaptic time;

%   tet=(linspace(0,2*pi,N))'; % initial conditions;
%   s=0; % initial conditions;
load('init_cond_micro_original_del_sig2.mat', 'tet','s'); % loading initial conditions;

ts=0.01*tau; % time window for computing the spiking rate;
tsp=10*ts; % time window for plotting the spiking rate;
tet0=pi-2*ts;
tet0p=pi-2*tsp;
eta0=(linspace(-1,1,N)).';
eta=d*eta0+muv;  % the full parameter of heterogeneity;
rm=zeros(Nt,1);
tm=linspace(0,(Nt-1)*dt,Nt);
Zm=linspace(Nt,1);
Zm=Zm.';
sk0=0;
tic;
for j=1:Nt
      sk0=sk0+1;
    if sk0==10000
        toc;
        disp([j Nt]);
        sk0=0;
        tic;
    end
   
    r=histcounts(tet,[tet0,pi])/(N*ts); % neurons spiking rate (for integration)
    rp=histcounts(tet,[tet0p,pi])/(N*tsp); % neurons spiking rate (for plotting)
    noise=sig*randn(N,1)*sqrt(dt); % generating random Gaussian noise;
    cs=cos(tet);
    sn=sin(tet);
    % numerical integration of microscopic equations:
    tet=tet+(dt*(1-cs+(1+cs).*(eta-J*s-0.5*sig^2*sn))+(1+cs).*noise)/tau;
    s=s+dt*(-s+r)/taus;
    tet=mod(tet,2*pi);
    rm(j)=rp;  % recording spiking rate value;
    Zm(j)=mean(exp(1i*tet));  % recording mean value of exponential;
end
toc;
disp(var(rm));
disp(var(abs(Zm)));

save('micro_original_del_sig2.mat','tm','rm'); % recording time series of the solution;
%-----------------------------------------------------------------------------------------
save('init_cond_micro_original_del_sig2.mat', 'tet','s'); % recording final state of the system;
% this command line should be commented if one does not wish to change the initial
% conditions; 
%-----------------------------------------------------------------------------------------

figure
subplot(4,1,1)
plot(tm,rm)
ylabel('r')
subplot(4,1,2)
plot(tm,real(Zm))
ylabel('Real(Z)')
subplot(4,1,3)
plot(tm,imag(Zm))
ylabel('Imag(Z)')
subplot(4,1,4)
plot(tm,abs(Zm))
ylabel('Abs(Z)')
xlabel('time')

figure
hold on
plot(real(Zm),imag(Zm))
%plot(real(Zt),imag(Zt),'r')
ylabel('Imag(Z)')
xlabel('Real(Z)')


