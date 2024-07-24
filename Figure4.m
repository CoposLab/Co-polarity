clear
close all
set(0,'DefaultFigureVisible','on')

% Figure 4B-C:
load('figure_data/figure4bc_data.mat')

outcome=Sup; %CoA produces Figure 4b; Sup produces figure 4c

labels_on=1;
title_on=1;
dot_size=2000;
line_width=2;

green=[0,1,0.5];
grey=[0.4,0.6,0.7];
blue=[0.2,0.2,0.9];
gg=[linspace(green(1),grey(1),20);linspace(green(2),grey(2),20);linspace(green(3),grey(3),20)];
gb=[linspace(grey(1),blue(1),20);linspace(grey(2),blue(2),20);linspace(grey(3),blue(3),20)];
ggb=flip([gg gb]');

f1=figure(1);
subplot(1,5,3)
    scatter3(Eb1(Ea1==Ea2),Ea1(Ea1==Ea2),Eb2(Ea1==Ea2),dot_size,outcome(Ea1==Ea2),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])
clim([0.2,1])
if labels_on==1
        xlabel('Eb1');
        ylabel('Ea1');
        zlabel('Eb2');
end
if title_on==1
        title('Ea1=Ea2')
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-12 12])
set(a,'YLim',[-12 12])
set(a,'ZLim',[-12 12])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,4)
    scatter3(Eb1(Ea1==-Ea2),Ea1(Ea1==-Ea2),Eb2(Ea1==-Ea2),dot_size,outcome(Ea1==-Ea2),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])

clim([0.2,1])
if labels_on==1
        xlabel('Eb1')
        ylabel('Ea1');
        zl=zlabel('Eb2');
end
if title_on==1
        title('Ea1=-Ea2','Position',[0,0,20])
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-12 12])
set(a,'YLim',[-12 12])
set(a,'ZLim',[-12 12])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,1)
    scatter3(Eb1(Ea1==0),Ea2(Ea1==0),Eb2(Ea1==0),dot_size,outcome(Ea1==0),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])

clim([0.2,1])
if labels_on==1
        xlabel('Eb1')
        ylabel('Ea1');
        zl=zlabel('Eb2');
end
if title_on==1
        title('Ea1=0','Position',[0,0,20])
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-12 12])
set(a,'YLim',[-12 12])
set(a,'ZLim',[-12 12])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,2)
    scatter3(Eb1(Ea2==0),Ea1(Ea2==0),Eb2(Ea2==0),dot_size,outcome(Ea2==0),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])
clim([0.2,1])
if labels_on==1
        xlabel('Eb1')
        ylabel('Ea1');
        zl=zlabel('Eb2');
end
if title_on==1
        title('Ea2=0','Position',[0,0,20])
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-12 12])
set(a,'YLim',[-12 12])
set(a,'ZLim',[-12 12])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,5)
cb=colorbar;
clim([0.2,1])
pos=get(cb,'position');
set(cb,'position',[pos(1)-0.05 pos(2)+0.2 pos(3)*2 pos(4)/2])
set(gca,'visible','off')
set(gcf,'color','w');

f1.Position = [584,599,1247,559];

if outcome==Sup
    title_text='Supracellular';
else
    title_text='Co-Aligned';
end
sgtitle(['Probability ' title_text])




% Figure 4E-F:
load('figure_data/figure4ef_data.mat')

outcome=Sup; %CoA produces Figure 4e; Sup produces figure 4f

f2=figure(2);
subplot(1,5,3)
    scatter3(kaa(kba==-kab),kab(kba==-kab),kbb(kba==-kab),dot_size,outcome(kba==-kab),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])
clim([0.2,1])
if labels_on==1
        xlabel('kaa');
        ylabel('kab');
        zlabel('kbb');
end
if title_on==1
        title('kab=-kba')
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-1.5 1.5])
set(a,'YLim',[-1.5 1.5])
set(a,'ZLim',[-1.5 1.5])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,4)
    scatter3(kaa(kba==0),kab(kba==0),kbb(kba==0),dot_size,outcome(kba==0),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])
clim([0.2,1])
if labels_on==1
        xlabel('kaa')
        ylabel('kab')
        zl=zlabel('kbb');
    % zl.Position(2) = zl.Position(2)-1;
end
if title_on==1
        title('kba=0','Position',[0,0,3])
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-1.5 1.5])
set(a,'YLim',[-1.5 1.5])
set(a,'ZLim',[-1.5 1.5])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,1)
    scatter3(kaa(kba==kab),kab(kba==kab),kbb(kba==kab),dot_size,outcome(kba==kab),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])
clim([0.2,1])
if labels_on==1
        xlabel('kaa')
        ylabel('kab')
        zl=zlabel('kbb');
end
if title_on==1
        title('kab=kba','Position',[0,0,3])
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-1.5 1.5])
set(a,'YLim',[-1.5 1.5])
set(a,'ZLim',[-1.5 1.5])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,2)
    scatter3(kaa(kab==0),kba(kab==0),kbb(kab==0),dot_size,outcome(kab==0),...
        'filled','LineWidth',line_width, 'MarkerEdgeColor',[0 0 0])
clim([0.2,1])
if labels_on==1
        xlabel('kaa')
        ylabel('kba')
        zl=zlabel('kbb');
end
if title_on==1
        title('kab=0','Position',[0,0,3])
end
colormap(ggb)
alpha 1
view(-52,30)
axis equal
a=gca;
set(a,'XLim',[-1.5 1.5])
set(a,'YLim',[-1.5 1.5])
set(a,'ZLim',[-1.5 1.5])
grid off
set(gca,'BoxStyle','full','LineWidth',line_width)
set(gca,'Ydir','reverse')
set(gcf,'color','w')

subplot(1,5,5)
cb=colorbar;
clim([0.2,1])
pos=get(cb,'position');
set(cb,'position',[pos(1)-0.05 pos(2)+0.2 pos(3)*2 pos(4)/2])
set(gca,'visible','off')
set(gcf,'color','w');

f2.Position = [584,599,1247,559];

sgtitle(['Probability ' title_text])
