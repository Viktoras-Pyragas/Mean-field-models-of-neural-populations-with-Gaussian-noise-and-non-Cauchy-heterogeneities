% Here we compute the poles ksi_j of normal distribution of normal heterogeneity,
% and the expansion coefficients b_j
function [ksi, b]=ksi_b(nn)
%close all
%format long;
[bet, C]=betn_Cn(nn); % finding the coefficient beta, and the constant Cn;
p=zeros(1,nn);
p(1)=1;
for j=2:nn
p(j)=p(j-1)/j;
end
Mfak=1.0/p(nn); % this is n!
p=flip(p);
%p=[1./p,1];
p=[p,1];
r0 = roots(p);
z=transpose((r0(imag(r0)>0)));
if (mod(nn,2)==0) % case of even n;
    ksi=[sqrt(z),-sqrt(conj(z))]/sqrt(bet);
    b=pi*Mfak*C*[-1i*z.^(-nn-1/2),1i*(conj(z)).^(-n-1/2)]/sqrt(bet);
else % case of odd n;
zM=r0(imag(r0)==0); % additional real negative root;
ksi=[sqrt(z),-sqrt(conj(z)),sqrt(zM)]/sqrt(bet);
b=pi*Mfak*C*[-1i*z.^(-nn-1/2),1i*(conj(z)).^(-nn-1/2),-1i*zM^(-nn-1/2)]/sqrt(bet);
end

ksi=ksi.';
b=b.';
end

