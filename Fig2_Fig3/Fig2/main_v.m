function main_v
clear all
close all
format long
%----------------------------------------------------------------------
% In this file we plot the graphs from different sources;
% -> from QIF meanfield model (two-cumulants approach);
% -> from truncated equations for Fourier components;
%----------------------------------------------------------------------
% Matcont result:
% J=2.0;
%----------------------------------------------------------------------
load('Xn1_Flat_J2.mat','X_Gd');
Xt1=X_Gd(:,1:2);
load('Xn2_Flat_J2.mat','X_Gd');
Xt2=X_Gd(:,1:2);
load('Xn6_Flat_J2.mat','X_Gd');
Xt6=X_Gd(:,1:2);
load('Xn10_Flat_J2.mat','X_Gd');
Xt10=X_Gd(:,1:2);
%----------------------------------------------------------------------
% Result from Fourier components:
% J=2.0;
%----------------------------------------------------------------------
load('H_n1_sp_J2.mat','Yt');
Yn1=Yt
load('H_n2_sp_J2.mat','Yt');
Yn2=Yt
load('H_n6_sp_J2.mat','Yt');
Yn6=Yt
load('H_n10_sp_J2.mat','Yt');
Yn10=Yt
%----------------------------------------------------------------------
% Matcont result:
% J=3.0;
%----------------------------------------------------------------------
load('Xn1_Flat_J3.mat','X_Gd');
Xt1p=X_Gd(:,1:2);
load('Xn2_Flat_J3.mat','X_Gd');
Xt2p=X_Gd(:,1:2);
load('Xn6_Flat_J3.mat','X_Gd');
Xt6p=X_Gd(:,1:2);
load('Xn10_Flat_J3.mat','X_Gd');
Xt10p=X_Gd(:,1:2);
%----------------------------------------------------------------------
% Result from Fourier components:
% J=3.0;
%----------------------------------------------------------------------
load('H_n1_sp_J3.mat','Yt');
Yn1p=Yt;
load('H_n2_sp_J3.mat','Yt');
Yn2p=Yt;
Yn2p=[Yn2p; Yn1p(end,:)];
load('H_n6_sp_J3.mat','Yt');
Yn6p=Yt;
Yn6p=[Yn6p; Yn1p(end,:)];
load('H_n10_sp_J3.mat','Yt');
Yn10p=Yt;
Yn10p=[Yn10p; Yn1p(end,:)];
%----------------------------------------------------------------------
% Matcont result:
% J=4.0;
%----------------------------------------------------------------------
load('Xn1_Flat_J4.mat','X_Gd');
Xt1w=X_Gd(:,1:2);
load('Xn2_Flat_J4.mat','X_Gd');
Xt2w=X_Gd(:,1:2);
load('Xn6_Flat_J4.mat','X_Gd');
Xt6w=X_Gd(:,1:2);
load('Xn10_Flat_J4.mat','X_Gd');
Xt10w=X_Gd(:,1:2);
%----------------------------------------------------------------------
% Result from Fourier components:
% J=4.0;
%----------------------------------------------------------------------
load('H_n1_sp_J4.mat','Yt');
Yn1w=Yt;
load('H_n2_sp_J4.mat','Yt');
Yn2w=Yt;
Yn2w=[Yn2w; Yn1w(end,:)];
load('H_n6_sp_J4.mat','Yt');
Yn6w=Yt;
Yn6w=[Yn6w; Yn1w(end,:)];
load('H_n10_sp_J4.mat','Yt');
Yn10w=Yt;
Yn10w=[Yn10w; Yn1w(end,:)];
%----------------------------------------------------------------------
% Matcont result:
% J=5.0;
%----------------------------------------------------------------------
load('Xn1_Flat_J5.mat','X_Gd');
Xt1k=X_Gd(:,1:2);
load('Xn2_Flat_J5.mat','X_Gd');
Xt2k=X_Gd(:,1:2);
load('Xn6_Flat_J5.mat','X_Gd');
Xt6k=X_Gd(:,1:2);
load('Xn10_Flat_J5.mat','X_Gd');
Xt10k=X_Gd(:,1:2);
%----------------------------------------------------------------------
% Result from Fourier components:
% J=5.0;
%----------------------------------------------------------------------
load('H_n1_sp_J5.mat','Yt');
Yn1k=Yt;
load('H_n2_sp_J5.mat','Yt');
Yn2k=Yt;
load('H_n6_sp_J5.mat','Yt');
Yn6k=Yt;
load('H_n10_sp_J5.mat','Yt');
Yn10k=Yt;
%----------------------------------------------------------------------
x1=0.4; % The coords where we write: J=(2,3,4,5)
y1=0.85;

nfk=8; % fonts in legend
mb1=2; % thickness of the gray line for mean-field solutions (MF)
mbs=1; % thickness of the colored lines for results obtaines from 
% equations for Fourier components

kc=0.6; % brightness of the gray line;
Xat=0.03; % coordinates of the text (a),(b),(c),(d); 
Yat=0.95;
%----------------------------------------------------------------------    
% circle; triangle; square (coordinates);
  Psw=[[0.04; 0.1], [0.046; 0.1], [0.04; 0.15]];
% circle -> a stable limit cycle  
% tringle, square -> stable fixed points;  
%----------------------------------------------------------------------    
nfa=10; % size of fonts in (a),(b),(c),(d);
nf=10;  % size of axis labels;
nfax=9; % size of numbers on axis;

figure
h1=subplot(2,2,1);
Yt=Yn1(1:end-1,1:2);
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq1=xq;
vq1=vq;
%---------------------------------------------
Yt=Yn2;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq2=xq;
vq2=vq;
%---------------------------------------------%---------------------------------------------
Yt=Yn6;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq6=xq;
vq6=vq;
%---------------------------------------------
Yt=Yn10;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq10=xq;
vq10=vq;
%---------------------------------------------J=2
hold on
xmn=0; xmx=0.06;
ymn=0; ymx=0.3;

% plotting the results from two-cumulants apprioach:
plot(Xt1(:,1),Xt1(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt2(:,1),Xt2(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt6(:,1),Xt6(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt10(:,1),Xt10(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);

% plotting results from the eqs. for Fourier components:
plot(xq1,vq1,'-.b','LineWidth',mbs);
plot(xq2,vq2,'--r','LineWidth',mbs);
plot(xq6,vq6,':b','LineWidth',mbs);
plot(xq10,vq10,'-r','LineWidth',mbs);
%---------------------------------------------------------------------

ms=6;
plot(Psw(1,1),Psw(2,1),'ok','MarkerSize',ms,'Color','black'); % circle
plot(Psw(1,2),Psw(2,2),'k^','MarkerSize',ms,'Color','black'); % triangle
plot(Psw(1,3),Psw(2,3),'square','MarkerSize',1.05*ms,'Color','black'); % square

xlim([xmn xmx]);
ylim([ymn ymx]);
%---------------------------------------------------------------------
text(x1,y1,'$\tilde{J}=2$','Interpreter','latex','Units','normalized','FontSize',nf);
ax = gca;
ax.XAxis.FontSize = nfax;
ax.YAxis.FontSize = nfax;
xlabel('$\tilde{\sigma}^2$','Interpreter','latex','FontSize',nf);
ylabel('$\tilde{\Delta}$','Interpreter','latex','FontSize',nf);
text(Xat,Yat,'(a)','Units','normalized','FontSize',nfa);


h2=subplot(2,2,2);
hold on
xmnp=0.0; xmxp=0.07;
ymnp=0.0; ymxp=0.45;
Yt=Yn1p(1:end-1,1:2);
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq1=xq;
vq1=vq;
%---------------------------------------------
Yt=Yn2p;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq2=xq;
vq2=vq;
%---------------------------------------------
Yt=Yn6p;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq6=xq;
vq6=vq;
%---------------------------------------------
Yt=Yn10p;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq10=xq;
vq10=vq;
%---------------------------------------------
% J=3
% plotting the results from two-cumulants apprioach:
plot(Xt1p(:,1),Xt1p(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt2p(:,1),Xt2p(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt6p(:,1),Xt6p(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt10p(:,1),Xt10p(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);

% plotting results from the eqs. for Fourier components:
plot(xq1,vq1,'-.b','LineWidth',mbs);
plot(xq2,vq2,'--r','LineWidth',mbs);
plot(xq6,vq6,':b','LineWidth',mbs);
plot(xq10,vq10,'-r','LineWidth',mbs);

xlim([xmnp xmxp]);
ylim([ymnp ymxp]);
text(x1,y1,'$\tilde{J}=3$','Interpreter','latex','Units','normalized','FontSize',nf);
ax = gca;
ax.XAxis.FontSize = nfax;
ax.YAxis.FontSize = nfax;
xlabel('$\tilde{\sigma}^2$','Interpreter','latex','FontSize',nf);
ylabel('$\tilde{\Delta}$','Interpreter','latex','FontSize',nf);
text(Xat,Yat,'(b)','Units','normalized','FontSize',nfa);

hold off
h3=subplot(2,2,3);
hold on
xmnp=0.0; xmxp=0.12;
ymnp=0.0; ymxp=0.55;

Yt=Yn1w;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq1=xq;
vq1=vq;
%---------------------------------------------
Yt=Yn2w;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq2=xq;
vq2=vq;
%---------------------------------------------
Yt=Yn6w;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq6=xq;
vq6=vq;
%---------------------------------------------
Yt=Yn10w;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq10=xq;
vq10=vq;
%---------------------------------------------
% J=4
% plotting the results from two-cumulants apprioach:
plot(Xt1w(:,1),Xt1w(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt2w(:,1),Xt2w(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt6w(:,1),Xt6w(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt10w(:,1),Xt10w(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);

% plotting results from the eqs. for Fourier components:
plot(xq1,vq1,'-.b','LineWidth',mbs);
plot(xq2,vq2,'--r','LineWidth',mbs);
plot(xq6,vq6,':b','LineWidth',mbs);
plot(xq10,vq10,'-r','LineWidth',mbs);

xlim([xmnp xmxp]);
ylim([ymnp ymxp]);
text(x1,y1,'$\tilde{J}=4$','Interpreter','latex','Units','normalized','FontSize',nf);
ax = gca;
ax.XAxis.FontSize = nfax;
ax.YAxis.FontSize = nfax;
xlabel('$\tilde{\sigma}^2$','Interpreter','latex','FontSize',nf);
ylabel('$\tilde{\Delta}$','Interpreter','latex','FontSize',nf);
text(Xat,Yat,'(c)','Units','normalized','FontSize',nfa);
xticks([0 0.05 0.1]);
xticklabels({'0','0.05','0.1'});
pos_leg_v=[0.405 0.38 0.1 0.1];
legend('MF (n=1,2,6,10)','','','',' "exact" (n=1)',' "exact" (n=2)',' "exact" (n=6)',' "exact" (n=10)','Position',pos_leg_v,'FontSize',nfk);
hold off

h4=subplot(2,2,4);
hold on
xmnp=0.0; xmxp=0.15;
ymnp=0.0; ymxp=0.68;

Yt=Yn1k;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq1=xq;
vq1=vq;
%---------------------------------------------
Yt=Yn2k;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq2=xq;
vq2=vq;
%---------------------------------------------
Yt=Yn6k;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq6=xq;
vq6=vq;
%---------------------------------------------
Yt=Yn10k;
Nq=1000;
method='spline';
[xq,vq]=interp_v2(Yt,Nq,method); % interpolating the results 
% from eqs. for Fourier components
xq10=xq;
vq10=vq;
%---------------------------------------------
% J=5
% plotting the results from two-cumulants apprioach:
plot(Xt1k(:,1),Xt1k(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt2k(:,1),Xt2k(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt6k(:,1),Xt6k(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);
plot(Xt10k(:,1),Xt10k(:,2),'-','color',kc*[1 1 1],'LineWidth',mb1);

% plotting results from the eqs. for Fourier components:
plot(xq1,vq1,'-.b','LineWidth',mbs);
plot(xq2,vq2,'--r','LineWidth',mbs);
plot(xq6,vq6,':b','LineWidth',mbs);
plot(xq10,vq10,'-r','LineWidth',mbs);

xlim([xmnp xmxp]);
ylim([ymnp ymxp]);
text(x1,y1,'$\tilde{J}=5$','Interpreter','latex','Units','normalized','FontSize',nf);
ax = gca;
ax.XAxis.FontSize = nfax;
ax.YAxis.FontSize = nfax;
yticks([0 0.2 0.4]);
yticklabels({'0','0.2','0.4'})
xlabel('$\tilde{\sigma}^2$','Interpreter','latex','FontSize',nf);
ylabel('$\tilde{\Delta}$','Interpreter','latex','FontSize',nf);
text(Xat,Yat,'(d)','Units','normalized','FontSize',nfa);
hold off

axes(h1);
hold on
Xv=[Psw(1,1) Psw(1,2)]; 
Yv=[Psw(2,1) Psw(2,2)];
myarrow_v(Xv,Yv); % plotting an arrow;
Xv = [Psw(1,1) Psw(1,3)]; 
Yv = [Psw(2,1) Psw(2,3)];
myarrow_v(Xv,Yv); % plotting an arrow;
hold off
hold off

set(gcf,'Units','centimeters')
set(gcf, 'PaperSize', [16.5,10.5]);
set(gcf,'Position',[2,2,16.5,10.5])
end