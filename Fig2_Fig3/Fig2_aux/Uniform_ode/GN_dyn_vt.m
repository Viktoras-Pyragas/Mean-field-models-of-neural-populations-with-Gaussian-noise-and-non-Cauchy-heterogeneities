%--------------------------------------------------------------------
% We integrate the ODEs of Fourier components with flat 
% heterogeneity;
% The ODEs of Fourier components is truncated; 
% We take into account the Gaussian noise (\sigma>0);
%--------------------------------------------------------------------


close all
clear all

format long;

global eta_k zeta tau taus J sig2 n M Ns njt Mt;

M=300; % number of Fourier components;
n=2; % number of poles of the flat heterogeneity distribution; 
Ns=n*M; % number of diff. eqs. related with Fourier components;
J=2.0; % global coupling strength
del=0.05; % HWHM of the heterogeneity distribution;
sig2=0.01;  % intensity of the Guassian noise;
muv=1.0; % average of eta (parameter of heterogeneity);
tau=1.0; % membrane time
taus=1.0; % time constant of spike relaxation
dt=0.01; % time step for output; 
Tpral=50; % time span of integration for warming up;
Tint=100; % time span of integration for output (to be recorded);
zeta=exp(1i*pi*(2*(1:n).'-1)/2/n); % poles of flat heterogeneity;
eta_k=del*zeta+muv; % array of parameters of heterogeneity at the poles; 
njt=(1:M).'; % indices of Fourier components;
Mt=(-1).^njt; % vector used for computation of averaged quantities;


load('GM_amt.mat','amt'); % loading initial conditions from data file;

  y0=zeros(Ns+1,1); % vector of the system state;
  % setting initial conditions from the data file:
  y0(1:Ns,1)=amt(1:Ns,1)+1i*amt(Ns+1:2*Ns,1);
  y0(Ns+1,1)=amt(2*Ns+1,1);
%     y0=zeros(1,Ns); % initial conditions for Fourier components; 
%     y0=[y0 0.0]; % initial conditions, taking into account equation of
%     relaxation;
 
odefun_vt=@MFE_vt; % functions on the r.h.s. of the ODEs;
%[t,y] = ode45(odefun,tspan,y0,options)
ode_sol=@ode45; % integrator of the ODEs;

rel_tol=1e-7; % relative tolerance;
abs_tol=1e-8; % absolute tolerance;
opts_v=odeset('RelTol',rel_tol,'AbsTol',abs_tol); % options for integrator

disp('Computing for output:');

 dTint=1.0; % small time span for integration; 
 tspan0=[0:dt:dTint]; % time span divided into dt-long steps
 szq=size(tspan0,2); % length of the time span;
 tspan=tspan0; % initial time span;

 Nt=Tint+Tpral; % number of large steps 
 
 NT=Tint; % total number of large steps for output
 nq=(1:szq-1);
 Yz1=zeros(NT*(szq-1),4); % in this array we will put the output
% Mt=(-1).^njt;
 tic;
 ng=0;
 
 for nT=1:Nt    
     ng=ng+1;
     if ng==1
     ng=0;
    toc;
    nT_v=nT % counting of computed large steps
    tic;
     end
[t,ym] = ode_sol(odefun_vt,tspan,y0,opts_v); % ode45
y0=ym(end,:); % updated initial conditions;
szy=size(ym,1);
if nT>Tpral % output is caslculated onle on the last Tint large steps;
Wr=zeros(szy,n); % auxiliary array;
for nv=1:n
    anj=ym(:,(nv-1)*M+njt.'); % Fourier components for nv-th pole;    
    Wr(:,nv)=anj*Mt; % auxiliary array;
end
if n==1 % for Cauchy heterogeneity
W=1+2*Wr; % finding the variable of the Lorentz Ansatz for the Cauchy heterogeneity;
else
W=1-2*1i*sin(pi/2/n)*zeta.'*Wr(:,1:n).'; % finding the variable of the Lorentz Ansatz;
W=W.';
end
Rv=real(W)/pi/tau; % averaged spiking rate;
Vv=-imag(W); % averaged potential;
% Here we accumulate the output:
Yz1(nq,1:4)=[t(1:end-1,1) Rv(1:end-1,1) Vv(1:end-1,1) ym(1:end-1,Ns+1)];
nq=nq+szq-1; % shifting the array of indices;
end % if nT>Tpral
tspan=tspan+dTint; % shifting the time span;

 end % for nT=1:Nt
 toc;

Yz1(1:10,:)

szw=size(ym)
save('GM_dyn.mat','Yz1','-mat'); % saving the output array in the data file;

y0=y0.';
amt=zeros(2*Ns+1,1);

amt(1:Ns,1)=real(y0(1:Ns,1));
amt(1+Ns:2*Ns,1)=imag(y0(1:Ns,1));
amt(2*Ns+1,1)=y0(Ns+1,1);

save('GM_amt.mat','amt','-mat'); % saving the final state of the system into a data file;

disp('Plotting the graphs:');

%Yz1(nq,1:4)=[t(1:end-1,1) Rv(1:end-1,1) Vv(1:end-1,1) ym(1:end-1,Ns+1)];
figure
subplot(3,1,1)
hold on
plot(Yz1(:,1),Yz1(:,2),'b-'); % dynamics of r(t);
hold off
subplot(3,1,2)
hold on
plot(Yz1(:,2),Yz1(:,3)); % phase portrait of (r(t),v(t))
hold off
subplot(3,1,3)
hold on
plot(Yz1(:,1),Yz1(:,3)); % dynamics of v(t);
hold off




