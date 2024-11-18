%----------------------------------------------------------------
% main_v.m
%-------------------------------------------------------------------------------------------
% Here we compute and plot graphs, in which are compared post-transient dynamics
% of microscopic and macroscopic models; for macroscopic model we chose
% n=6, and for microscopic - the normal heterogeneity distribution;
% Here we considered the case of the normal heterogeneity;
%-------------------------------------------------------------------------------------------
close all
format long;

%dta=-1.45; % shifting of the t-axis (for MF equations);
dta=0; % shifting of the t-axis (for MF equations);
% laiko asies stumimas (MF lygciu)
dtam=-1.2; % shifting of the t-axis (for microscopic equations);
% laiko asies stumimas (mikrosk. lygciu)
my_color=0.7*[1 1 1]; % color of the solution of MF equations;
% MF sprendinio spalva
%--------------------------------------------------
% Line widths:
lw1=0.5; % for solution of microscopic equations;
lw2=4; % for MF equations;
%--------------------------------------------------
micro='-r'; % the color of line for microscopic solutions;
%--------------------------------------------------
% (x,y) positions of texts (a), (b), (c)
Xat=1; 
Yat=0.45; 
% positions of parameters [sigma^2,Delta] values:
Xpa=8; Ypa=0.44; % (a) graph
Xpb=8; Ypb=0.44; % (b) graph
Xpc=8; Ypc=0.44; % (c) graph
%--------------------------------------------------
% x, y limits
miny=0.0;
maxy=0.5;  
minx=0;
maxx=30;  
%--------------------------------------------------
% Positions of legends:
pos_a=[0.349 0.89 0.38 0.04]; % in (a) graph
pos_b=[0.349 0.57 0.38 0.04]; % in (b) graph
pos_c=[0.349 0.25 0.38 0.04]; % in (c) graph
%--------------------------------------------------
% Loading solution for MF and microscopic equations, for
% system parameters: [sigma^2, Delta]=[0.02,0.26];
load('mean_field_original_del_sig2.mat','t','r');
load('micro_original_del_sig2.mat','tm','rm');
%-----------------------------------------------------------    
figure
% plotting the general graph:
[ha, pos]=tight_subplot_v(3,1,[0.04 0.001],[0.08 0.02],[0.13 0.05]);

Na=1; % [a] number of the graph axis;
axes(ha(Na));
hold on
plot(t+dta,r,'-','Color',my_color,'LineWidth',lw2) % solution of MF equations;
plot(tm+dtam,rm,micro,'LineWidth',lw1); % solution of microscopic equations;
hold off
xlim([0,maxx])
ylim([miny,maxy])
ylabel('$\tilde{R}$','Interpreter','latex','FontSize',10);
text(Xat,Yat,'(a)','FontSize',10);
text(Xpa,Ypa,'$\tilde{\sigma}^2=0.02 \quad\tilde{\Delta}=0.26$','Interpreter','latex','FontSize',10);
legend({'MF (n=6)','Micro (normal)'},'Position',pos_a,'NumColumns',2);
legend('boxoff');

set(ha(Na),'XTickLabel','');
set(ha(Na),'XScale','linear');
set(ha(Na),'YScale','linear');

set(ha(Na),'XTick',[0 10 20 30 40 50]);
set(ha(Na),'YTick',[0 0.25 0.5]) 
set(ha(Na),'YTickLabel',{'0','0.25','0.5','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
box on
%------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------
% Loading solution for MF and microscopic equations, for
% system parameters: [sigma^2, Delta]=[0.05,0.26];
load('mean_field_changed_sig2.mat','t','r');
load('micro_changed_sig2.mat','tm','rm');
%load('micro_changed_sig2_pr.mat','tm','rm');
delt=0.3; % shifting the solution of the MF equations;
delr=1.4; % shifting the solution of the microscopic equations;
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
Na=2; % [b] number of the graph axis;
axes(ha(Na));
hold on
plot(t,r,'-','Color',my_color,'LineWidth',lw2) % solution of MF equations;
plot(tm,rm,micro,'LineWidth',lw1) % solution of microscopic equations;
hold off

set(ha(Na),'XTickLabel','');
set(ha(Na),'XScale','linear');
set(ha(Na),'YScale','linear');
ylabel('$\tilde{R}$','Interpreter','latex','FontSize',10);
text(Xat,Yat,'(b)','FontSize',10);
text(Xpb,Ypb,'$\tilde{\sigma}^2=0.05 \quad\tilde{\Delta}=0.26$','Interpreter','latex','FontSize',10);
legend({'MF (n=6)','Micro (normal)'},'Position',pos_b,'NumColumns',2);
legend('boxoff');

set(ha(Na),'XTick',[0 10 20 30 40 50]);
set(ha(Na),'YTick',[0 0.25 0.5]) 
set(ha(Na),'YTickLabel',{'0','0.25','0.5','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
xlim([0,maxx])
ylim([miny,maxy])
 box on
%------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
Na=3; % [c] number of the graph axis;
% Loading solution for MF and microscopic equations, for
% system parameters: [sigma^2, Delta]=[0.02,0.31];
load('mean_field_changed_del.mat','t','r');
load('micro_changed_del.mat','tm','rm');
axes(ha(Na));
hold on
plot(t,r,'-','Color',my_color,'LineWidth',lw2)
plot(tm,rm,micro,'LineWidth',lw1)
hold off
ylabel('$\tilde{R}$','Interpreter','latex','FontSize',10);
xlim([0,maxx])
ylim([miny,maxy])
text(Xat,Yat,'(c)','FontSize',10);
text(Xpc,Ypc,'$\tilde{\sigma}^2=0.02 \quad\tilde{\Delta}=0.31$','Interpreter','latex','FontSize',10);
legend({'MF (n=6)','Micro (normal)'},'Position',pos_c,'NumColumns',2);
legend('boxoff');

set(ha(Na),'XTick',[0 10 20 30 40 50]);
set(ha(Na),'XTickLabel',{'0','10','20','30','40','50','Interpreter','latex','FontSize',10})
xlabel('$\tilde{t}$','Interpreter','latex','FontSize',10);
set(ha(Na),'XScale','linear');
set(ha(Na),'YScale','linear');

set(ha(Na),'YTick',[0 0.25 0.5]) 
set(ha(Na),'YTickLabel',{'0','0.25','0.5','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
box on
%------------------------------------------------------------------------------------------------
% %--------------------------------------------------------------------------------
set(gcf,'Units','centimeters')
set(gcf, 'PaperSize', [8.5 10]);
set(gcf,'Position',[2,2,8.5,10]);

%  set(gcf,'renderer','Painters');
%  print -depsc -tiff -r300 -painters Fig_nec_dyn.eps