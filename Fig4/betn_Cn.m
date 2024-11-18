% Skaiciujamas betM ir normavimo koef CM bet kokiam M (lyginiam ir
% nelyginiam)
%Computing coefficeints betM and CM for any n (even and odd);
function [bet, C]=betn_Cn(n)
close all
%M=1;
bet=betn(n); % finding coefficient bet;
%disp(bet)
p=zeros(1,n);
p(1)=1;
for j=2:n
p(j)=p(j-1)*j;
end
nfak=p(n); % cia yra M!
p=flip(p); % flipping the array;
p=[1./p,1]; % coefficients of polynomial;
r0 = roots(p); % finding roots of polynomial;
z=r0(imag(r0)>0); % extracting roots with positive imaginary parts;
length(z)
S=sum(z.^(-n-0.5));
if (mod(n,2)==0)
% for even n: 
Cn=sqrt(bet)/(2*pi*nfak*imag(S)); 
else
% for odd n:    
zn=r0(imag(r0)==0); % additional real-valued negative root;
Cn=sqrt(bet)/(2*pi*nfak*(imag(S)+0.5*(abs(zn))^(-n-0.5)));
end
C=Cn;
end

% This fuction computes the coefficient bet;
function betan=betn(n) 
p=zeros(1,n);
p(1)=1;
for j=2:n
p(j)=p(j-1)*j;
end
p=flip(p); % flipping the array;
p=[1./p,-1]; % coefficients of the polynomial;
r = roots(p); % computing roots of polynomial;
r=r(imag(r)==0); % extracting real-valued roots;
betan=r(r>0); % a single real-valued root of the polynomial 
end