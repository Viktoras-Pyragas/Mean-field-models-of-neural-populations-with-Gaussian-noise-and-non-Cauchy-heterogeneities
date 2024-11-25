%--------------------------------------------------------------------
% This program performs the Newton iterations for the FPE rewritten
% in the Fourier components; 
% Given initial coordinates of the
% fixed point at the given location on the parameters plane, the program
% finds the corrected coordinated at the slighly shifted position on the
% parameters plane;
% At the start we fix a linear path on the parameters plane, consisting of
% small steps; the program performs the Newton iterations on each step, ant
% moves the values of the coordinates from step to step, until the
% whole path is covered;
% At each step we calculate the two leading eigenvalues of the Jacobi
% matrix, both their real and imaginary parts; 
% While composing the Jacobi matrix at the fixed point, this program uses
% auxiliary functions which calculate the necessary fragments of the
% Jacobi matrix;
% The linear function, involved in the Newton iterations, is also computed
% by using auxiliary functions, which evaluate the fragments of the 
% nonlinear vector-function
%--------------------------------------------------------------------close all
clear all

format long;

global eta_k zeta bvk tau taus J sig2 del n M Ns njt Mt;

M=300; % number of Fourier components;
n=2; % numper of poles of the normal distribution of heterogeneity;
Ns=n*M; % order of the ODE, without the relaxation equation;
J=2; % global coupling strength
del=5e-3; % HWHM of the normal distribution of heterogeneity;
sig2=0.05; % strength of the Gaussian noise;
muv=1.0; % average of eta
tau=1.0; % membrane time
taus=1.0; % time constant of spike relaxation
nn=n;
[ksi, b]=ksi_b(nn);
zeta=ksi % poles of normal heterogeneity;
bvk=b; % expansion coefficients for averaging over the heterogeneity distribution;
eta_k=del*zeta+muv; % array of parameters of heterogeneity at the poles;

njt=(1:M).'; % indices of Fourier components;
Mt=(-1).^njt; % vector used for computation of averaged quantities;
 
%load('GM_amt1v.mat','amt');
%load('GM_amt2v.mat','amt'); % Gv2=0.0; del=0.8; J=2.0; taus=1.0; n=2;
%load('GM_amt2_t.mat','amt'); 
%load('GM_amt_d01_Gv0_v.mat','amt'); % Gv2=0.0; del=0.1; J=2.0; taus=1.0; n=2; 
%load('GM_amt_d01_Gv001_v.mat','amt'); % Gv2=0.01; del=0.1; J=2.0; taus=1.0; n=2; 
%load('GM_amt_d001_Gv001_v.mat','amt'); % Gv2=0.01; del=0.01; J=2.0; taus=1.0; n=2; 
%load('GM_amt_d001_Gv003_v.mat','amt'); % Gv2=0.03; del=0.01; J=2.0; taus=1.0; n=2; 
%load('GM_amt_d001_Gv004_v.mat','amt'); % Gv2=0.04; del=0.01; J=2.0; taus=1.0; n=2; 
%load('GM_amt_d01_Gv005_v.mat','amt'); % Gv2=0.05; del=0.1; J=2.0; taus=1.0; n=2;
%load('GM_amt_d001_Gv005_v.mat','amt'); % Gv2=0.05; del=0.01; J=2.0; taus=1.0; n=2;
load('GM_amt_d5e-3_Gv005_v.mat','amt'); % Gv2=0.05; del=5e-3; J=2.0; taus=1.0; n=2;
% amp=[amt(1:end-1,1); amt(1:end-1,2); amt(end,1)];
% amt=amp;

% am=zeros(n*Nj+1,1);
% am(1:n*Nj,1)=amt(1:n*Nj,1)+1i*amt(n*Nj+1:2*n*Nj,1);
% am(n*Nj+1,1)=amt(2*n*Nj+1,1);
% 
% am_v=am(50:60,1)

%ur=pr
disp('Computing for output:');
Nt=3; % number of points to be computed
Nk=5; % number of Newton iterations for each point;
      Gvmn=0.05; % initial value of system parameter;
      Gvmx=0.08; % final value of system parameter;
%        Lvmn=0.01;
%        Lvmx=0.005;
ng=0; % counter
tic;
LMtot=zeros(Nt+1,5); % array of the parameters and of two corresponding 
% eigenvalues: two largest real parts and corresponding imaginary parts;
%LMtot_del=zeros(Nt+1,5);
%Lvk=Lvmn+(Lvmx-Lvmn)*(0:Nt)/Nt;
%Lvk=[Lvk 0.005];
%szL=size(Lvk,2);
%LMtot_del=zeros(szL,5);
tic;
Nkv=0;
 for nT=0:Nt 
 %for nT=1:szL 
     sig2=Gvmn+(Gvmx-Gvmn)*nT/Nt; % scanning the system parameter;
     %Lv=Lvmn+(Lvmx-Lvmn)*nT/Nt;
     %Lv=Lvk(1,nT)
     %Gv2=Gv^2;
     %del=Lv; % current system parameter;
     eta_k=del*zeta+muv; % current parameter of heterogeneity;

     ng=ng+1; % counter
     if ng==1
     ng=0;
    toc;
    Gv2_v=sig2 % showing current system parameter
%   Lv_v=Lv
%   nT_v=nT
%   Nk_v=Nkv
    nT_Nk_v=[nT Nkv]
    tic;
     end

for nk=1:Nk
%while dmx2 > 1e-4    
%nk=nk+1;
nk_nT=[nk nT]
Phin=zeros(2*M*n+1,2*M*n+1); % the total Jacobi matrix;
FN=zeros(2*M*n+1,1); % the total nonlinear function
St=amt(end,1); % variable of relaxation;
njb=(1:M); % array of indices;
for np=1:n
    amp=amt(njb,1)+1i*amt(njb+n*M,1); % complex-valued Fourier components for corrent pole;
    [LM1, Fm1] = LM1_v(amp,St,np); % finding the current fragment of the jacobi matrix, and
    % the current part of the nonlinear function;
    Bv=Bv1_v(np); % finding the current fragment of the last row of the Jacobi matrix;
    Lv1=Lv1_v(amp); % finding the current fragment of the last column of the Jacobi matrix;
    %----------------------------------------------------
    % composing current fragments of the Jacobi matrix:
    Phin(njb,njb)=real(LM1);
    Phin(njb+n*M,njb+n*M)=real(LM1);
    Phin(njb,njb+n*M)=-imag(LM1);
    Phin(njb+n*M,njb)=imag(LM1);
    Phin(2*M*n+1,njb)=2*real(Bv);
    Phin(2*M*n+1,njb+n*M)=-2*imag(Bv);
    Phin(njb,2*M*n+1)=real(Lv1);
    Phin(njb+n*M,2*M*n+1)=imag(Lv1);
    %----------------------------------------------------
    % composing the current fragments of the nonlinear function:
    FN(njb,1)=real(Fm1);
    FN(njb+n*M,1)=imag(Fm1);
    %----------------------------------------------------
    njb=njb+M; % shifting the array of indices;
end % for n=1:np
%----------------------------------------------------
% Completing the last elements of the Jacobi matrix and 
% the nonlinear function:
Phin(2*M*n+1,2*M*n+1)=-1/taus;
FN(2*n*M+1,1)=FSt_v(amt);
%----------------------------------------------------

Ln=eig(Phin);  % finding eigenvalues of the Jacobi matrix;
Lnr=real(Ln);  % finding the real parts;
Lnim=imag(Ln); % finding the imaginary parts;
[Bn,IX]=sort(Lnr,'descend'); % sorting the real parts in descending order;
Lnx=[Lnr(IX(1:end)) Lnim(IX(1:end))]; % composing eigenvalues in the sorted order;
Lnx_v=Lnx(1:2,:) % showing the two first eigenvalues possessing the largest real parts;
damt=Phin\FN; % finding the current correction for the Newton iterations;

damv=damt(1:n*M,1)+1i*damt(n*M+1:2*n*M,1);
damv=[damv; damt(2*M*n+1,1)];
dmx=abs(damv);
dmx=mean(dmx) % showing the averaged correction of all the variables;

amt=amt-damt; % correcting thefixed point coordinates (Newton iteration);
am=zeros(n*M+1,1);
am(1:n*M,1)=amt(1:n*M,1)+1i*amt(n*M+1:2*n*M,1);
am(n*M+1,1)=amt(2*n*M+1,1);
amT=abs(am);
amT=mean(amT) % showing the averaged system's variables of the fixed point;
dmx2=dmx/amT % the relation of the averaged fixed point correction 
% to the averaged fixed point;
FNt=FN(1:n*M,1)+1i*FN(n*M+1:2*n*M,1);
FNt=[FNt; FN(2*n*M+1,1)];
Fnx=abs(FNt(1:n*M+1,1));
Fnx=mean(Fnx) % averaged nonlinear function (it should converge to zero)
end
amT_v=amT % average of fixed point varibles 
dmx_v=dmx % average of the fixed point correction;
dmx2_v=dmx2 % relative correction of the averaged fixed point;
Fnx_v=Fnx % averaged nonlinear function
LMtot(nT+1,:)=[sig2 Lnx(1:2,1).' Lnx(1:2,2).']; % the two leading eigenvalues 
% with corresponding scanned parameter value (the other parameter remains fixed);
%LMtot_del(nT+1,:)=[del Lnx(1:2,1).' Lnx(1:2,2).'];
%LMtot_del(nT,:)=[del Lnx(1:2,1).' Lnx(1:2,2).'];
 end % for nT=0:Nt
toc;

  save('GM_amt2.mat','amt','-mat'); % recording the coordinates 
% of the fixed point after iterations
  save('GM_LMtot.mat','LMtot','-mat'); % puting the results into a data file
% save('GM_LMtot_del.mat','LMtot_del','-mat'); % puting the results into a data file




