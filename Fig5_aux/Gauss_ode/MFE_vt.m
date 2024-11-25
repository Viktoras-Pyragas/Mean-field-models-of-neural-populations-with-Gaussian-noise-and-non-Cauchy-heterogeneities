% r.h.s of ODE for the Fourier components and realxation variable'
function dydt = MFE_vt(t,y)
global tau taus bvk eta_k J sig2 n M Ns Mt njt;
% n -> number of poles of the heterogeneity function;
% M -> number of Fourier components;
% Ns=n*M; % number of system variables (without relaxation variable);
% y(1:Ns+1,1); % vector of system state;
dydt=zeros(Ns+1,1); % array of derivatives;
% njt=(1:M).'; % array of indices
% Mt=(-1).^njt; % auxiliary array;
Wr=zeros(n,1); % auxiliary array;
% Finding the averaged spiking rate:
for nv=1:n
    anj=y((nv-1)*M+njt.',1); % Fourier amplitudes;
    Wr(nv,1)=anj.'*Mt;
end
if n==1
    W=1+2*Wr; % the case of Cauchy heterogeneity;
else
W=1+2*bvk.'*Wr; % the general case of normal heterogeneity;
end
Rv=real(W)/pi/tau; % averaged spiking rate; 
A=eta_k-J*tau*y(Ns+1,1); % auxiliary variable;
% hv=0.5*1i*(A-1.0); % auxilary variable;
% Omv=A+1.0; % auxilary variable;

    for nv=1:n   
    % Composing amplitudes vectors with pole number nv;
    % Fourier amplitudes (for pole number nv):
    ajk=y((nv-1)*M+njt.',1);     
    % Shifted arrays of Fourier amplitudes:
    ajkm1=[1.0; ajk(1:M-1,1)]; % a(m-1); % a0=1, due to normalization;    
    ajkp1=[ajk(2:M,1); 0.0]; % a(m+1); % the last is truncated; 
    ajkm2=[0; 1.0; ajk(1:M-2,1)]; % a(m-2); % the first is truncated; a0=1, due to normalization;
    ajkp2=[ajk(3:M,1); 0.0; 0.0]; % a(m+2); % the two last are truncated;

    dydtv=zeros(M,1); % array of derivatives;
    % Deterministic terms:
    dydtp=1i*ajk*(A(nv,1)+1);
    dydtv=dydtv+dydtp;
    dydtp=(ajkm1+ajkp1)*0.5*1i*(A(nv,1)-1.0);   
    dydtv=dydtv+dydtp;
    % Terms related with Gaussian noise:
    dydtp=-0.5*sig2*(3.0*njt.*ajk+(2*njt+1.0).*ajkp1+(2*njt-1).*ajkm1);
    dydtv=dydtv+dydtp;
    dydtp=-0.25*sig2*((1.0+njt).*ajkp2+(njt-1).*ajkm2);
    dydtv=dydtv+dydtp;
    dydt((nv-1)*M+njt.',1)=njt.*dydtv(:,1);
    end
    dydt(Ns+1,1)=(-y(Ns+1,1)+Rv)/taus; % equation of relaxation;
    dydt(1:Ns,1)=dydt(1:Ns,1)/tau; % equations for Fourier components;
end    














