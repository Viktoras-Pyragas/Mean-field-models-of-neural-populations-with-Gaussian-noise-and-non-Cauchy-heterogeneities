% Plotting the dependence of leading real parts of the eigenvalues vs a
% system parameter;
% Here we also find the value of system parameter at which the real part of
% leading eigenvalue crosses zero; this is a point of Hopf bifurcation; 
function plot_GN_eigvs
format long
%----------------------------------------------------------------------
% n=2:
%-------------------------------------------------------------------------------------------
load('GM_LMtot.mat','LMtot'); % spektralines lygtys; Nj=200;
%load('GM_LMtot_d2_Gv0-01.mat','LMtot');
%load('GM_LMtot_del.mat','LMtot_del'); % spektralines lygtys; Nj=300;
%load('GM_LMtot_del_d1.8-1.5_Gv01_v.mat','LMtot_del'); % del=0.38; Gv2=0.35; J=20; n=10;
%load('GM_LMtot_del_d1.5-1.4_Gv01_v.mat','LMtot_del'); % del=0.25; Gv2=0.0; J=2.0;
%load('GM_LMtot_del_d1.2-0.9_Gv02_v.mat','LMtot_del'); % del=0.25:-0.01:0.14; Gv2=0.01; J=2.0;
%load('GM_LMtot_d1.4_Gv01-02_v.mat','LMtot'); % del=0.25:-0.01:0.1; Gv2=0.02; J=2.0;
%load('GM_LMtot_del_Gv003.mat','LMtot_del'); % del=0.15:-0.01:0.1; Gv2=0.03; J=2.0;
%load('GM_LMtot_del_Gv004.mat','LMtot_del'); % del=0.15:-0.01:0.05; Gv2=0.04; J=2.0;
%load('GM_LMtot_d002.mat','LMtot'); % del=0.02; Gv2=0.06:-0.001:0.04; J=2.0;
%load('GM_LMtot_d001.mat','LMtot'); % del=0.01; Gv2=0.06:-0.001:0.04; J=2.0;
%load('GM_LMtot_d0.mat','LMtot'); % del=0.0; Gv2=0.06:-0.001:0.04; J=2.0;
%LMtot=LMtot_del;
%LMtot(nT+1,:)=[Gv Lnx(1:2,1).' Lnx(1:2,2).'];
%Gv=LMtot(:,1);
% Flip the array if the scanning of parameter was backwards:
% for np=1:5
%     LMtot(:,np)=flip(LMtot(:,np));
% end

szL=size(LMtot,1)

Lv=LMtot(:,1);  % array of the scanned parameter;
L1=[LMtot(:,2) LMtot(:,4)]; % the real and imaginary parts of the first eigenvalue;
L2=[LMtot(:,3) LMtot(:,5)]; % the real and imaginary parts of the second eigenvalue;

% Lvmn=1.5;
% Lvmx=1.2;
% Nt=15;
% Lv=Lvmn+(Lvmx-Lvmn)*(0:Nt).'/Nt;
% Lv=flip(Lv);
%-------------------------------------------------------------------------------------------
%Zr=[[Lv(1,1); Lv(end,1)] zeros(2,1)]; % location of the line at zero ordinates;
% figure
% hold on
% %plot(Gv,L1(:,1),'b*');
% plot(Lv,L1(:,1),'b*'); % the real part of the leading eigenvalue;
% plot(Zr(:,1),Zr(:,2),'m--'); % zero line
% %plot(Gv(:,1),L2(:,1),'r*');
% hold off

% Gvmn=0.2;
%       Gvmx=0.3;
%     Lvmn=1.04;
%     Lvmx=0.94;
% Lvt=(1.04:-0.02:0.94).';
% Lvt=flip(Lvt);
% L1t=flip(L1);
% Gvt=(0.2:0.02:0.3).';

% Finding the parameter value at which 
% the leading real part crosses zero:
sz1=size(Lv,1)
 %sz1=size(Lvt,1)
for nw=1:sz1-1
    Lg1=L1(nw,1);
    Lg2=L1(nw+1,1);
    if Lg1*Lg2<0
        m1=nw;
        m2=nw+1;
    end
end
%m1=1; m2=2;

%Pt=[Lv L1(:,1)]
nt=[m1 m2]
n1=nt(1);
n2=nt(2);
%  x1=Lv(n1,1)
%  y1=L1(n1,1)
%  x2=Lv(n2,1)
%  y2=L1(n2,1)

x1=Lv(n1,1)
y1=L1(n1,1)
x2=Lv(n2,1)
y2=L1(n2,1)



x3=x1-y1*(x2-x1)/(y2-y1) % linear interpolation;
xv=[x3 0]; % the value of system parameter at which the
% real part of the leading eigenvalue crosses zero;

% sz1=size(Gvt,1)
% for nw=1:sz1-1
%     Lg1=L1(nw,1);
%     Lg2=L1(nw+1,1);
%     if Lg1*Lg2<0
%         m1=nw;
%         m2=nw+1;
%     end
% end
%m1=1; m2=2;
% 
% %Pt=[Lv L1(:,1)]
% nt=[m1 m2]
% n1=nt(1);
% n2=nt(2);

% x1t=Gvt(n1,1)
% x2t=Gvt(n2,1)
% y1=L1(n1,1)
% y2=L1(n2,1)
% 
% x3t=x1t-y1*(x2t-x1t)/(y2-y1)
% xvt=[x3t 0];

% figure
% hold on
% %plot(Gv,L1(:,1),'b*');
% plot(Lv,L1(:,1),'b*');
% plot(Zr(:,1),Zr(:,2),'m--');
% plot(xv(1,1),xv(1,2),'r*');
% %plot(Gv(:,1),L2(:,1),'r*');
% hold off

%Lvt=flip(Lvt);
%L1t=flip(L1);
Hr=[[0.42; 0.42],[-0.1; 0.1]];
Hl=[[0.38; 0.38],[-0.1; 0.1]];

Zr=[[Lv(1,1); Lv(end,1)],[0; 0]]; % location of the line at zero ordinates;
figure
hold on
%plot(Gv,L1(:,1),'b*');
plot(Lv,L1(:,1),'b-'); % the real part of the leading eigenvalue;
plot(Zr(:,1),Zr(:,2),'m--'); % zero line;
%plot(Hr(:,1),Hr(:,2),'m--');
%plot(Hl(:,1),Hl(:,2),'m--');
plot(xv(1,1),xv(1,2),'r*'); % the point where eigv crosses zero;
%plot(Gv(:,1),L2(:,1),'r*');
%xlim([0.2 0.8]);
hold off

%Zr=[[Gvt(1,1); Gvt(end,1)] zeros(2,1)];
% figure
% hold on
% %plot(Gv,L1(:,1),'b*');
% plot(Gvt,L1(:,1),'b*');
% plot(Zr(:,1),Zr(:,2),'m--');
% plot(xvt(1,1),xvt(1,2),'r*');
% %plot(Gv(:,1),L2(:,1),'r*');
% hold off


% figure
% hold on
% % plot(Gv,L1(:,2),'b*');
% % plot(Gv(:,1),L2(:,2),'r*');
% plot(Lv,L1(:,2),'b*');
% plot(Lv(:,1),L2(:,2),'r*');
% hold off
% 
% Zr=[min(L1(:,1)); max(L1(:,1))];
% Zr=[Zr, zeros(2,1)];
% figure
% hold on
% plot(L1(:,1),L1(:,2),'b*');
% plot(L2(:,1),L2(:,2),'r*');
% plot(Zr(:,1),Zr(:,2),'m--');
% hold off


Zr_v=Zr


end