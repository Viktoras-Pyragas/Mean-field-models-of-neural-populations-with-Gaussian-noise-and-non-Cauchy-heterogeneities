% Here we compute betan and normalization coefficient Cn for any n
% (even and odd);
function [bet, C]=betn_Cn(nn)
bet=bet_n(nn);
p=zeros(1,nn);
p(1)=1;
for j=2:nn
p(j)=p(j-1)/j;
end
Mfak=1.0/p(nn); % this is n!
p=flip(p);
%p=[1./p,1];
p=[p,1];
r0 = roots(p)
z=r0(imag(r0)>0);
S=sum(z.^(-nn-0.5))
if (mod(nn,2)==0)
CM=sqrt(bet)/(2*pi*Mfak*imag(S));
else
zM=r0(imag(r0)==0); % additional negative root;
CM=sqrt(bet)/(2*pi*Mfak*(imag(S)+0.5*(abs(zM))^(-nn-0.5)));
end
C=CM;
end

