% This function interpolates the data from the Hopf curve
function [xq,vq]=interp_v2(Yt,Nq,method)

x=Yt(:,2); % array of del
v=Yt(:,1); % array of sig2
%Nq -> number of del points to be interpolated;
a1=x(1,1); a2=0.0; % the bounds of del interval
xq=a1+(a2-a1)*(0:Nq)/Nq; % the interval [a1,a2] (of del) of interpolation;
% method='spline'; -> method of interpolation;
vq = interp1(x,v,xq,method); % interpolating in the points xq (of del);
% exchanging the arrays;
xt=xq;
xq=vq; % sig2 points;
vq=xt; % del points;

end