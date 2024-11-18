function Bv = Bv1_v(np)
global taus bvk Mt;
% n -> number of poles;
% M -> number of Fourier components;
% Ns=n*M; % the order of ODE (without equation of relaxation);

%njt=(1:M).'; % array of indices;
%Mt=(-1).^njt; % auxiliary vector;
Bv=bvk(np,1)*Mt.'/pi/taus; % the complex-valued row-vector
%----------------------------------------------------------------------------
% As the simplest example:
% (with n=1):
% Bv1r=real(Bv); % real part;
% Bv1im=imag(Bv); % imaginary part;
%------------------------------------------------------------------
% Full Jacobi matrix (with n=1):
% LMT=[Lmr    -Lmim    Lv1r;
%      Lmim    Lmr     Lv1im;
%      2*Bv1r -2*Bv1im Lv2];
%------------------------------------------------------------------
end    














