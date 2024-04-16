clear
close all
set(0,'DefaultFigureVisible','on')

% load('data/racrho_binding_unbinding.mat') %lf
% bindunbind=vals1;
% load('data/nosignal_racrho_yes.mat');
load('data/concdependent_racrho_nosignal.mat')


green=[0,1,0.5];
grey=[0.4,0.6,0.7];
blue=[0.2,0.2,0.9];
gg=[linspace(green(1),grey(1),20);linspace(green(2),grey(2),20);linspace(green(3),grey(3),20)];
gb=[linspace(grey(1),blue(1),20);linspace(grey(2),blue(2),20);linspace(grey(3),blue(3),20)];
ggb=flip([gg gb]');

yvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
xvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};

figure(1)
% h=heatmap(xvals,yvals,bindunbind); %lf
% h=heatmap(xvals,yvals,bindunbindyes'); %yes
h=heatmap(xvals,yvals,bindingunbindinglf);
h.Title='Proportion Leader/Follower';
% h.Title='Proportion Yes';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')

% figure(1)
% imagesc(bindunbind)
% colorbar;
% colormap(ggb)
% clim([0.2,1])
% set(gcf,'color','w')
% set(gca,'XTick',[],'YTick',[])


% load('data/racrho_binding.mat') %lf

xvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};
yvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};

figure(2)
% h=heatmap(xvals,yvals,bindingvals); %lf
% h=heatmap(xvals,yvals,bindingyes'); %yes
h=heatmap(xvals,yvals,bindinglf);
h.Title='Proportion Leader/Follower';
% h.Title='Proportion Yes';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')

% load('data/racrho_unbinding.mat') %lf

xvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
yvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};

figure(3)
% h=heatmap(xvals,yvals,unbindingvals'); %lf
% h=heatmap(xvals,yvals,unbindingyes); %yes
h=heatmap(xvals,yvals,unbindinglf);
h.Title='Proportion Leader/Follower';
% h.Title='Proportion Yes';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')