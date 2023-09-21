clear
close all;
cla

load('branched_bundled_heatmap_vals.mat')

green=[0,1,0.5];
grey=[0.4,0.6,0.7];
blue=[0.2,0.2,0.9];
gg=[linspace(green(1),grey(1),20);linspace(green(2),grey(2),20);linspace(green(3),grey(3),20)];
gb=[linspace(grey(1),blue(1),20);linspace(grey(2),blue(2),20);linspace(grey(3),blue(3),20)];
ggb=flip([gg gb]');

f=figure(1);
subplot(1,4,1)
scatter3(ka(kc==-0.9),kb(kc==-0.9),kd(kc==-0.9),1000,lf(kc==-0.9),'filled')
clim([0,max(lf)])
xlabel('ka')
ylabel('kb')
zlabel('kd')
title('kc=-0.9','Position',[0,0,2])
colormap(ggb)
alpha 1
view(-52.8368,30.7651)
axis equal
set(gca,'XTick',[],'YTick',[],'ZTick',[])
grid on

subplot(1,4,2)
scatter3(ka(kc==0),kb(kc==0),kd(kc==0),1000,lf(kc==0),'filled')
clim([0,max(lf)])
xlabel('ka')
ylabel('kb')
zlabel('kd')
title('kc=0','Position',[0,0,2])
colormap(ggb)
alpha 1
view(-52.8368,30.7651)
axis equal
set(gca,'XTick',[],'YTick',[],'ZTick',[])
grid on

subplot(1,4,3)
scatter3(ka(kc==0.9),kb(kc==0.9),kd(kc==0.9),1000,lf(kc==0.9),'filled')
clim([0,max(lf)])
xlabel('ka')
ylabel('kb')
zlabel('kd')
title('kc=0.9','Position',[0,0,2])
colormap(ggb)
alpha 1
view(-52.8368,30.7651)
axis equal
set(gca,'XTick',[],'YTick',[],'ZTick',[])
grid on

subplot(1,4,4)
cb=colorbar;
clim([0,max(lf)])
pos=get(cb,'position');
set(cb,'position',[pos(1)-0.1 pos(2)+0.2 pos(3)*2 pos(4)/2])
set(gca,'visible','off')
set(gcf,'color','w');

f.Position = [584,599,1247,559];