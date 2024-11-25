% Plotting the dependence of leading real parts of the eigenvalues vs a
% system parameter;
% Here we also find the value of system parameter at which the real part of
% leading eigenvalue crosses zero; this is a point of Hopf bifurcation; 
function plot_GN_eigvs
format long
%----------------------------------------------------------------------
% n=2:
%-------------------------------------------------------------------------------------------
%load('GM_LMtot.mat','LMtot'); % from equations for Fourier components; M=200;
load('GM_LMtot_del.mat','LMtot_del'); % from equations for Fourier components; M=200;
%load('GM_LMtot_del_d025_Gv0.mat','LMtot_del'); % del=0.25; Gv2=0.0; J=2.0;
%load('GM_LMtot_del_Gv001.mat','LMtot_del'); % del=0.25:-0.01:0.14; Gv2=0.01; J=2.0;
%load('GM_LMtot_del_Gv002.mat','LMtot_del'); % del=0.25:-0.01:0.1; Gv2=0.02; J=2.0;
%load('GM_LMtot_del_Gv003.mat','LMtot_del'); % del=0.15:-0.01:0.1; Gv2=0.03; J=2.0;
%load('GM_LMtot_del_Gv004.mat','LMtot_del'); % del=0.15:-0.01:0.05; Gv2=0.04; J=2.0;
%load('GM_LMtot_d002.mat','LMtot'); % del=0.02; Gv2=0.06:-0.001:0.04; J=2.0;
%load('GM_LMtot_d001.mat','LMtot'); % del=0.01; Gv2=0.06:-0.001:0.04; J=2.0;
%load('GM_LMtot_d0.mat','LMtot'); % del=0.0; Gv2=0.06:-0.001:0.04; J=2.0;
LMtot=LMtot_del;
% Flip the array if the scanning of parameter was backwards:
for np=1:5
    LMtot(:,np)=flip(LMtot(:,np));
end

szL=size(LMtot,1)

Lv=LMtot(:,1); % array of the scanned parameter;
L1=[LMtot(:,2) LMtot(:,4)]; % the real and imaginary parts of the first eigenvalue;
L2=[LMtot(:,3) LMtot(:,5)]; % the real and imaginary parts of the second eigenvalue;

%-------------------------------------------------------------------------------------------
Zr=[[Lv(1,1); Lv(end,1)] zeros(2,1)]; % location of the line at zero ordinates;
figure
hold on
%plot(Gv,L1(:,1),'b*');
plot(Lv,L1(:,1),'b*'); % the real part of the leading eigenvalue;
plot(Zr(:,1),Zr(:,2),'m--'); % zero line
%plot(Gv(:,1),L2(:,1),'r*');
hold off

% Finding the parameter value at which 
% the leading real part crosses zero:
sz1=size(Lv,1)
for nw=1:sz1-1
    Lg1=L1(nw,1);
    Lg2=L1(nw+1,1);
    if Lg1*Lg2<0
        m1=nw;
        m2=nw+1;
    end
end

Pt=[Lv L1(:,1)]
nt=[m1 m2]
n1=nt(1);
n2=nt(2);
x1=Lv(n1,1)
y1=L1(n1,1)
x2=Lv(n2,1)
y2=L1(n2,1)

x3=x1-y1*(x2-x1)/(y2-y1) % linear interpolation;
xv=[x3 0]; % the value of system parameter at which the
% real part of the leading eigenvalue crosses zero;

figure
hold on
%plot(Gv,L1(:,1),'b*');
plot(Lv,L1(:,1),'b*'); % the real part of the leading eigenvalue;
plot(Zr(:,1),Zr(:,2),'m--'); % zero line;
plot(xv(1,1),xv(1,2),'r*'); % the point where eigv crosses zero;
hold off

end