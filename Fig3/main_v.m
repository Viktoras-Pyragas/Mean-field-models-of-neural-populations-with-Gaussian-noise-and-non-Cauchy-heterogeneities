%----------------------------------------------------------------
% main_v.m
%-------------------------------------------------------------------------------------------
% Here we compute and plot graphs, in which are compared post-transient dynamics
% of microscopic and macroscopic models; for macroscopic model we chose
% n=10, and for microscopic - the rectangular heterogeneity distribution;
% Here we considered the case of the flat heterogeneity;
%----------------------------------------------------------------
format long;
close all
dta=-3.33; % shift of the t-axis in (a) graph
my_color=0.7*[1 1 1]; % color of the MF line
%--------------------------------------------------
% Line widths:
lw1=0.5; % for solution of the microscopic equations;
lw2=4; %  for solution of the MF equations;
%--------------------------------------------------
micro='-r'; % the color of line for microscopic solutions;
%--------------------------------------------------
% (x,y) positions of texts (a), (b), (c)
Xat=1; 
Yat=0.45;
%--------------------------------------------------
xp=8;
yp=0.43;
% (a), (b), (c) positions for the values of parameters [sigma^2,Delta]:
Xpa=xp; Ypa=yp;
Xpb=xp; Ypb=yp;
Xpc=xp; Ypc=yp;
%--------------------------------------------------
% x, y limits
miny=0.0; % 
maxy=0.5; % 
minx=0;
maxx=30; % 
%--------------------------------------------------
% Loading solution for MF and microscopic equations, for
% system parameters: [sigma^2, Delta]=[0.04, 0.1];
load('mean_field_original_del_sig2.mat','t','r'); % loadig MF solution
load('micro_original_del_sig2.mat','tm','rm'); % loading microscopic equations solution;
%-----------------------------------------------------------    
figure
% plotting the general graph:
[ha, pos]=tight_subplot_v(3,1,[0.04 0.001],[0.08 0.02],[0.13 0.05]);

Na=1; % [a] number of the graph axis;
axes(ha(Na));
hold on
plot(t+dta,r,'-','Color',my_color,'LineWidth',lw2) % solution of MF equations;
plot(tm,rm,micro,'LineWidth',lw1); % solution of microscopic equations;
hold off
xlim([0,maxx])
ylim([miny,maxy])
ylabel('$\tilde{R}$','Interpreter','latex','FontSize',10);
text(Xat,Yat,'(a)','FontSize',10);
text(Xpa,Ypa,'$\tilde{\sigma}^2=0.04 \quad\tilde{\Delta}=0.1$','Interpreter','latex','FontSize',10);
legend({'MF (n=10)','Micro (uniform)'},'Location','south','NumColumns',2);
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
%--------------------------------------------------
% Loading solution for MF and microscopic equations, for
% system parameters: [sigma^2, Delta]=[0.046, 0.1];
load('mean_field_changed_sig2.mat','t','r'); % solution of MF equations;
load('micro_changed_sig2.mat','tm','rm'); % solution of microscopic equations;
% delt=0.3; % stumiame suvidurkinta sprentini
% delr=1.4; % stumiame mikroskopini sprendini
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
text(Xpb,Ypb,'$\tilde{\sigma}^2=0.046 \quad\tilde{\Delta}=0.1$','Interpreter','latex','FontSize',10);
legend({'MF (n=10)','Micro (uniform)'},'Location','south','NumColumns',2);
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
%--------------------------------------------------
Na=3; % [c] number of the graph axis;
%clear 
% Loading solution for MF and microscopic equations, for
% system parameters: [sigma^2, Delta]=[0.04, 0.15];
load('mean_field_changed_del.mat','t','r'); % loading solution of MF equations;
load('micro_changed_del.mat','tm','rm'); % loading solution of microscopic equations;
axes(ha(Na));
hold on
plot(t,r,'-','Color',my_color,'LineWidth',lw2)
plot(tm,rm,micro,'LineWidth',lw1)
hold off
ylabel('$\tilde{R}$','Interpreter','latex','FontSize',10);
xlim([0,maxx])
ylim([miny,maxy])
text(Xat,Yat,'(c)','FontSize',10);
text(Xpc,Ypc,'$\tilde{\sigma}^2=0.04 \quad\tilde{\Delta}=0.15$','Interpreter','latex','FontSize',10);
legend({'MF (n=10)','Micro (uniform)'},'Location','south','NumColumns',2);
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
%--------------------------------------------------------------------------------
set(gcf,'Units','centimeters')
set(gcf, 'PaperSize', [8.5 10]);
set(gcf,'Position',[2,2,8.5,10]);


%  set(gcf,'renderer','Painters');
%  print -depsc -tiff -r300 -vector uniform_dyn.eps