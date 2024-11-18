%--------------------------------------------------------------------
% We integrate the ODEs of Fourier components with normal 
% heterogeneity;
% The ODEs of Fourier components is truncated; 
% We take into account the Gaussian noise (\sigma>0);
%--------------------------------------------------------------------
clear all

format long;

global eta_k zeta bvk tau taus J sig2 n M Ns njt Mt;

M=300; % number of Fourier components;
n=6; % number of poles in the normal heterogeneity;
Ns=n*M; % order of ODEs without the Eq. of relaxation;
J=5.0; % global coupling strength
del=0.5; % HWHM of the normal distribution of the normal heterogeneity;
%sig2=\sigma^2; % square of the noise dispersion
sig2=0.0; % strength of the noise;
etav=1.0; % average of eta; (parameter of heterogeneity); 
tau=1.0; % membrane time
taus=1.0; % time constant of spike relaxation
dt=0.01; % time step for output; 
Tpral=50; % time span of integration for warming up;
Tint=100; % time span of integration for output (to be recorded);
M=n;
[ksi, b]=ksi_b_v(M);
zeta=ksi; % poles of normal heterogeneity;
bvk=b; % expansion coefficients for averaging over the heterogeneity distribution;
eta_k=del*zeta+etav; % array of parameters of heterogeneity at the poles; 
njt=(1:M).'; % indices of Fourier components;
Mt=(-1).^njt; % vector used for computation of averaged quantities;

%load('GM_amt.mat','amt'); % loading initial conditions from data file;

%amt_v=amt(900:910,:)
%ur=pr
  y0=zeros(Ns+1,1); % vector of the system state;
%-----------------------------------------------------------  
% setting initial conditions from the data file:
%   y0(1:Ns,1)=amt(1:Ns,1)+1i*amt(Ns+1:2*Ns,1);
%   y0(Ns+1,1)=amt(2*Ns+1,1);
%-----------------------------------------------------------  
%  y0=zeros(Ns,1); % initial conditions for Fourier components; 
%  y0=[y0; 0.0]; % initial conditions, taking into account equation of
%     relaxation;
 
odefun_vt=@MFE_vt; % functions on the r.h.s. of the ODEs;
ode_sol=@ode45; % integrator of the ODEs;

rel_tol=1e-7; % relative tolerance;
abs_tol=1e-8; % absolute tolerance;
opts_v=odeset('RelTol',rel_tol,'AbsTol',abs_tol); % options for integrator

disp('Computing for output:');

 dTint=1.0; % small time span for integration; 
 tspan0=[0:dt:dTint]; % time span divided into dt-long steps
 szq=size(tspan0,2); % length of the time span;
 tspan=tspan0; % initial time span;

 Nt=Tint+Tpral; % total number of large steps 
 
 NT=Tint; % number of large steps for output
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
    nT_v=nT  % counting of computed large steps
    tic;
     end
[t,ym] = ode_sol(odefun_vt,tspan,y0,opts_v); % ode45
y0=ym(end,:); % updated initial conditions;
szy=size(ym,1);
if nT>Tpral % output is caslculated onle on the last Tint large steps;
Wr=zeros(szy,n);
for nv=1:n
    anj=ym(:,(nv-1)*M+njt.'); % spektralines anplitudes;    
    Wr(:,nv)=anj*Mt;
end
if n==1 % for Cauchy heterogeneity
W=1+2*Wr; % finding the variable of the Lorentz Ansatz for the Cauchy heterogeneity;
else
W=1+2*bvk.'*Wr(:,1:n).'; % finding the variable of the Lorentz Ansatz;
W=W.';
end
Rv=real(W)/pi/tau; % averaged spiking rate;
Vv=-imag(W); % averaged potential;
% Here we accumulate the output:
Yz1(nq,1:4)=[t(1:end-1,1) Rv(1:end-1,1) Vv(1:end-1,1) ym(1:end-1,Ns+1)];
nq=nq+szq-1;
end % if nT>Tpral
tspan=tspan+dTint;
 end % for nT=1:Nt
 toc;

Yz1(1:10,:)

szw=size(ym)
save('GM_dyn.mat','Yz1','-mat'); % saving the output array in the data file;

%amt=[real(y0(1,1:Nj).'); imag(y0(1,1:Nj).'); y0(1,Nj+1)];
%amt=[];
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




