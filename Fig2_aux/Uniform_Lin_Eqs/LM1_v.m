function [LM1, Fm] = LM1_v(amt,St,np)
global eta_k J Gv2 M Mt njt tau;
% n -> number of poles;
% np -> number of specific pole (numeration);
% M -> number of Fourier components;
% LM1 -> submatrix for one pole (np-th)
% Fm -> fragment of the nonlinear function corresponding to the np-th pole;
am=amt(:,1); % Fourier components for one pole;
% St -> synaptic activation; 
etk=eta_k(np,1); % parameter of heterogeneity for one pole;
A=etk-J*tau*St; % auxiliary variable;
%Omv=etk-J*St+1;
%hv=0.5*1i*(etk-J*St-1);

LM1=zeros(M,M); % fragment of the Jacobi matrix;
njt=(1:M).'; % array of indices for Fourier components;
Mt=(-1).^njt; % auxiliary array;
%----------------------------------------------------------------------------
% Computing a fragment of the Jacobi matrix:
%----------------------------------------------------------------------------
for mp=1:M
   LMt=zeros(1,M); % auxiliary array;
   %am0=am(mp,1);
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
   if (mp<(M-1))&&(mp>2)
   LMt(1,mp-2:mp+2)=[-0.25*(mp-1)*Gv2, 0.5*1i*(A-1)-Gv2*(mp-0.5), (1i*(A+1)-1.5*Gv2*mp), (0.5*1i*(A-1)-Gv2*(mp+0.5)),-0.25*(mp+1)*Gv2];
   end
   if mp==(M-1)
   LMt(1,mp-2:mp+1)=[-0.25*(mp-1)*Gv2, 0.5*1i*(A-1)-Gv2*(mp-0.5), (1i*(A+1)-1.5*Gv2*mp), (0.5*1i*(A-1)-Gv2*(mp+0.5))];
   end
   if mp==M
   LMt(1,mp-2:mp)=[-0.25*(mp-1)*Gv2, 0.5*1i*(A-1)-Gv2*(mp-0.5), (1i*(A+1)-1.5*Gv2*mp)];
   end
   if mp==2
   LMt(1,mp-1:mp+2)=[0.5*1i*(A-1)-Gv2*(mp-0.5), (1i*(A+1)-1.5*Gv2*mp), (0.5*1i*(A-1)-Gv2*(mp+0.5)),-0.25*(mp+1)*Gv2]; 
   end
   if mp==1
   LMt(1,mp:mp+2)=[(1i*(A+1)-1.5*Gv2*mp), (0.5*1i*(A-1)-Gv2*(mp+0.5)),-0.25*(mp+1)*Gv2];
   end
   LM1(mp,:)=mp*LMt;
end
% LM1 -> submatrix for one pole (np-th)
%------------------------------------------------------------------
% Full Jacobian (with n=1):
% LMT=[Lmr -Lmim Lv1r;
%      Lmim Lmr Lv1im;
%      2*LM2 zeros(1,M) Lv2];
%------------------------------------------------------------------
% Computing the fragment of the r.h.s. nonlinear function:
%------------------------------------------------------------------
Fm=zeros(M,1); % fragment of the r.h.s. of the ODE
for mp=1:M
%amw=[am2 am1 am0 ap1 ap2].';
am0=am(mp,1);
if mp==1
    am2=0;
    am1=1;
end
if mp==2
    am2=1;
    am1=am(mp-1,1);
end
if mp>2
    am2=am(mp-2,1);
    am1=am(mp-1,1);
end
if mp<(M-1)
    ap1=am(mp+1,1);
    ap2=am(mp+2,1);
end
if mp==(M-1)
    ap1=am(mp+1,1);
    ap2=0;
end
if mp==M
    ap1=0;
    ap2=0;
end
%Fm(mp,[mp-2:mp+2])=[-0.25*Gv2*mp*(mp-1),mp*(hv-Gv2*(mp-0.5)),mp*(1i*Omv-1.5*Gv2*mp),mp*(hv-Gv2*(mp+0.5)),-0.25*Gv2*mp*(mp+1)];
Fm(mp,1)=mp*sum([-am2*0.25*Gv2*(mp-1),am1*(0.5*1i*(A-1)-Gv2*(mp-0.5)),am0*(1i*(A+1)-1.5*Gv2*mp),ap1*(0.5*1i*(A-1)-Gv2*(mp+0.5)),-0.25*ap2*Gv2*(mp+1)]);
end
% Fr=real(Fm);
% Fim=imag(Fm);
%------------------------------------------------------------------
% r.h.s. of Eqs. :
%FN=[Fr; Fim];
%------------------------------------------------------------------

end    














