Yt=[];

Gv=0.0;
del=0.181218934142783;
Yt=[Yt; [Gv del]];

Gv=0.005;
del=0.169859384651373;
Yt=[Yt; [Gv del]];

Gv=0.01;
del=0.158155445536807; % vk2

Yt=[Yt; [Gv del]];

Gv=0.015;
del=0.146028739645944;

Yt=[Yt; [Gv del]];


Gv=0.02;
del=0.133386975791483; % vk1

Yt=[Yt; [Gv del]];

Gv=0.025;
del=0.120112556140481;

Yt=[Yt; [Gv del]];

Gv=0.03; % vk3
del=0.105916298486529;

Yt=[Yt; [Gv del]];

Gv=0.035;
del=0.090604815916409;

Yt=[Yt; [Gv del]];

Gv=0.04; % vk4
del=0.073398365317410;

Yt=[Yt; [Gv del]];

Gv=0.045;
del=0.052996543188284;

Yt=[Yt; [Gv del]];

Gv=0.0475;
del=0.040439017749046;

Yt=[Yt; [Gv del]];

Gv=0.050338763808142; % vk5
del=0.02;

Yt=[Yt; [Gv del]];


Gv=0.051077694509153;
del=0.01;

Yt=[Yt; [Gv del]];

Gv=0.051330062303728;
del=0;

Yt=[Yt; [Gv del]];

figure
hold on
plot(Yt(:,1),Yt(:,2),'b*');
xlim([0 0.06]);
ylim([0 0.2]);
hold off

save('Hopf_biff_n2_v.mat','Yt','-mat');