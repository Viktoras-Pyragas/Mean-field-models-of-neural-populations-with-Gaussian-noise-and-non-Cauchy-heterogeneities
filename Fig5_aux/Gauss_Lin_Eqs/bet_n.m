%find bet
function betan=bet_n(nn)
format long;
%M=4;
p=zeros(1,nn);
p(1)=1;
for j=2:nn
p(j)=p(j-1)*j;
end
p=flip(p);
p=[1./p,-1];
%p=[p,-1];
r = roots(p)
r=r(imag(r)==0);
betan=r(r>0)

end