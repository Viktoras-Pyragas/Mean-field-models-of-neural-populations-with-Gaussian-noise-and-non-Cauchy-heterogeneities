function Lv1 = Lv1_v(amt)
global J M tau;
% n -> number of poles;
% M -> number of Fourier components;
% J -> coupling strength;
am=amt; % Fourier components corresponding to one of the poles;
Lv1=zeros(M,1); % the column-vector
%----------------------------------------------------------------------------
% Composing the fragment of the Jacobi matrix:
%----------------------------------------------------------------------------
for mp=1:M
   am0=am(mp,1);
   if mp>1
       am1=am(mp-1,1);
   else
       am1=1;
   end
   if mp<M
       ap1=am(mp+1,1);
   else
       ap1=0;
   end
   Lv1(mp,1)=-mp*1i*tau*J*(am0+0.5*am1+0.5*ap1);
end
%----------------------------------------------------------------------------
% As the simplest example:
% (with n=1):
% Lv1r=real(Lv1); % real part;
% Lv1im=imag(Lv1); % imaginary part;
%------------------------------------------------------------------
% Full Jacobi matrix (with n=1):
% LMT=[Lmr    -Lmim    Lv1r;
%      Lmim    Lmr     Lv1im;
%      2*Bv1r -2*Bv1im Lv2];
%------------------------------------------------------------------


end    














