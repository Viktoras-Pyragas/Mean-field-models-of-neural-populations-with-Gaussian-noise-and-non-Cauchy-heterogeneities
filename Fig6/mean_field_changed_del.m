% Normal-type distribution of heterogeinity eta
% Synaptic equation involved;
% Here we use the two-cumulants approach 
% with the (a1,a2) variables;
close all
clearvars
format long;

del=0.31; % half width at half maximum (HWHM);
J=5;  % coupling strength;
muv=1; % average of heterogeneity parameter eta;
tau=1; % membane time
taus=1; % synaptic time
n=6; % number of poles
sig2=0.02; % intensity of the Gaussian noise
nn=n;
[ksi, bk]=ksi_b(nn);
zeta=ksi.';
b=bk.';
P.del=del; P.J=J; P.muv=muv; P.tau=tau; P.taus=taus; P.zeta=zeta; P.b=b; P.sig2=sig2;
Tpral=3000; % warming up time;
Tint=60; % time of integration
N=1000;
x=linspace(-5, 5,N);
C=(n/pi)*sin(pi/(2*n));
g=C./(1+x.^(2*n));
gL=(1/pi)*1./(1+x.^2); % Caushy (Lorenz) distribution

figure
hold on;
plot(x,g,'-b')
plot(x,gL,'-r')


%wm0=zeros(2*n+1,1); % initial condition;
 a1=1*rand(n,1)-1i*rand(n,1); % initial condition;
 a2=a1.^2;
 wm0=[a1;a2;0];

tic;
options=odeset('RelTol',1e-7,'AbsTol',1e-12);
[~,wm] = ode45(@(t,wm) MFE(t,wm,P), [0, Tpral], wm0,options);
wm0=(wm(end,:)).';
toc;
tic;
[t,wm] = ode45(@(t,wm) MFE(t,wm,P), [0, Tint], wm0,options);
a1=wm(:,1:n);
a2=wm(:,n+1:2*n);
kap=a2-a1.*a1;
w=((1-a1)./(1+a1)+2*kap./(1+a1).^3)*b.';
toc;

r=real(w)/(pi*tau);
v=-imag(w);
%disp(wm(end,:))
figure
plot(v,r)

save('mean_field_changed_del.mat','t','r') 

figure
plot(t,r)
ylabel('$r$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')

figure
plot(t,real((2*kap./(1+a1).^3)*b.')/(pi*tau))

function dydt = MFE(~,wm,P)
% r.h.s. of the two-cumulants equations;
del=P.del; J=P.J; muv=P.muv; tau=P.tau; taus=P.taus; zeta=P.zeta; b=P.b; sig2=P.sig2;
n=(length(wm)-1)/2;
a1=wm(1:n);
a2=wm(n+1:2*n);
a3=3*a1.*a2-2*a1.^3;
a4=a1.^2.*(6*a2-5*a1.^2);
kap=a2-a1.*a1;
w=b*((1-a1)./(1+a1)+2*kap./(1+a1).^3);
r=real(w)/(tau*pi);
A=muv+del*zeta.'-J*tau*wm(2*n+1);
dydt1=1i*((A+1).*a1+(A-1).*(1+a2)/2)-sig2*(1.5*a1+0.5+1.5*a2+0.5*a3);
dydt2=1i*2*((A+1).*a2+(A-1).*(a1+a3)/2)-sig2*(6*a2+3*a1+5*a3+0.5+1.5*a4);
dydt3=(-wm(2*n+1)+r)/taus;
dydt=[dydt1/tau;dydt2/tau;dydt3];  
end 


