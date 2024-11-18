% Computing the last element of the nonlinear function;
function FSt = FSt_v(amt)
global taus zeta n M Mt;
% amt -> all Fourier components (for all poles) and the variable St;
% n -> number of poles;
% M -> number of Fourier components;
% Ns=n*M; % the order of ODE (without equation of relaxation);
%am=amt(1:end-1,1); % harmonics of all poles
am=amt(1:n*M,1)+1i*amt(n*M+1:2*n*M,1); % complex-valued Fourier components for all poles;
St=amt(2*n*M+1,1); % relaxation variable St;
%----------------------------------------------------------------------------
% Full Jacobian (for n=1):
% LMT=[Lmr -Lmim Lv1r;
%      Lmim Lmr Lv1im;
%      2*LM2 zeros(1,M) Lv2];
%------------------------------------------------------------------
% Computing the last component of the r.h.s. function:
%------------------------------------------------------------------
AML=zeros(M,n); % we put all Fourier components into one matrix:
ntv=(1:M); % array of indices of Fourier components;
for np=1:n % scanning the poles
    AML(:,np)=am(ntv,1);
    ntv=ntv+M;
end
Rt=real(-1i*sin(pi/2/n)*zeta(:,1).'*(1+2*Mt.'*AML).')/pi/tau; % global mean firing rate;
FSt=(-St+Rt)/taus; % the last component of the r.h.s. of the ODE
%------------------------------------------------------------------
end    














