set(0,'DefaultFigureVisible','on')
close all;
signal=0;


load('./vid_matfiles/collision/racup/1000RacOn1.mat')
a1=a1all(:,end-1);
a2=a2all(:,end-1);
b1=b1all(:,end-1);
b2=b2all(:,end-1);
t=2499;
L=10;
Nt=2500;
adjacent=1;

%Define colors
colorLength = 50;
white = [1,1,1];
darkyellow = [227/256,180/256,76/256];
yellow = [254/256,254/256,98/256];
pink = [211/256,95/256,183/256];
darkpink = [141/256,45/256,113/256];

whiteyellow = [linspace(white(1),yellow(1),colorLength)',linspace(white(2),yellow(2),colorLength)',linspace(white(3),yellow(3),colorLength)'];
yellowdarkyellow = [linspace(yellow(1),darkyellow(1),colorLength)',linspace(yellow(2),darkyellow(2),colorLength)',linspace(yellow(3),darkyellow(3),colorLength)'];
whitedarkyellow = [whiteyellow;yellowdarkyellow];
whitepink = [linspace(white(1),pink(1),colorLength)',linspace(white(2),pink(2),colorLength)',linspace(white(3),pink(3),colorLength)'];
pinkdarkpink = [linspace(pink(1),darkpink(1),colorLength)',linspace(pink(2),darkpink(2),colorLength)',linspace(pink(3),darkpink(3),colorLength)'];
whitedarkpink = [whitepink;pinkdarkpink];


branchedColor = whitedarkpink;
bundledColor = whitedarkyellow;
branchedColName = 'Pink';
bundledColName = 'Yellow';

% Find median for cell 1
a1New = a1;
a1New(a1New<1)=0;
if (a1New(1)~=0 && a1New(end)~=0)
    zeroInd1=find(a1New==0,1,'first');
    zeroInd2=find(a1New==0,1,'last');
    dirIndexa1=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a1New~=0,1,'first');
    ind2=find(a1New~=0,1,'last');
    dirIndexa1=ceil((ind1+ind2)/2);
end
if dirIndexa1<1
    dirIndexa1=dirIndexa1+101;
end
b1New = b1;
b1New(b1New<1)=0;
if (b1New(1)~=0 && b1New(end)~=0)
    zeroIndb1_1=find(b1New==0,1,'first');
    zeroIndb2_1=find(b1New==0,1,'last');
    dirIndexb1=ceil((zeroIndb1_1+zeroIndb2_1)/2) - 50;
else
    indb1_1=find(b1New~=0,1,'first');
    indb2_1=find(b1New~=0,1,'last');
    dirIndexb1=ceil((indb1_1+indb2_1)/2);
end
if dirIndexb1<1
    dirIndexb1=dirIndexb1+101;
end
% Find median for cell 2
a2New = a2;
a2New(a2New<1)=0;
if (a2New(1)~=0 && a2New(length(a2New))~=0)
    zeroInd1=find(a2New==0,1,'first');
    zeroInd2=find(a2New==0,1,'last');
    dirIndexa2=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a2New~=0,1,'first');
    ind2=find(a2New~=0,1,'last');
    dirIndexa2=ceil((ind1+ind2)/2);
end
if dirIndexa2<1
    dirIndexa2=dirIndexa2+101;
end
b2New = b2;
b2New(b2New<1)=0;
if (b2New(1)~=0 && b2New(end)~=0)
    zeroIndb1_2=find(b2New==0,1,'first');
    zeroIndb2_2=find(b2New==0,1,'last');
    dirIndexb2=ceil((zeroIndb1_2+zeroIndb2_2)/2) - 50;
else
    indb1_2=find(b2New~=0,1,'first');
    indb2_2=find(b2New~=0,1,'last');
    dirIndexb2=ceil((indb1_2+indb2_2)/2);
end
if dirIndexb2<1
    dirIndexb2=dirIndexb2+101;
end


% Make scatterplots
scatplot=figure(1);
subplot(1,2,1); %Cell 1
plot(Xa,a1,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
plot(Xa,b1,'-ok','color',bundledColor(end,:),'linewidth',3);
plot(s1,xC1,'-.','color',branchedColor(end,:),'linewidth',1);
plot(s1,yC1,'-.k','color',bundledColor(end,:),'linewidth',1);
% xlim([0 10]); ylim([0 2]);
%title('Time = 0');
set(gca,'fontname','times','fontsize',20); box on;
lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
lgd.NumColumns = 2;
set(gcf,'color','w');
title('Cell 1')
hold off;
%keyboard
% pause(1.0);

subplot(1,2,2); %Cell 2
plot(Xa,a2,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
plot(Xa,b2,'-ok','color',bundledColor(end,:),'linewidth',3);
plot(s2,xC2,'-.','color',branchedColor(end,:),'linewidth',1);
plot(s2,yC2,'-.k','color',bundledColor(end,:),'linewidth',1);
% xlim([0 10]); ylim([0 2]);
%title('Time = 0');
set(gca,'fontname','times','fontsize',20); box on;
lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
lgd.NumColumns = 2;
set(gcf,'color','w');
title('Cell 2')
hold off;
%keyboard
% pause(1.0);


% Define circles
gapsize=0.01;
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.85:0.01:1);
[Xcol,Ycol] = pol2cart(th,rad);
Ycol1=Ycol;
Ycol2=Ycol;
Ycol1(:,boundC1)=Ycol1(:,boundC1(1)*ones(1,length(boundC1)));
Ycol2(:,boundC2)=Ycol2(:,boundC2(1)*ones(1,length(boundC2)));
Ycol2 = Ycol2 - 2*abs(max(max(Ycol2)))-gapsize;
ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]';
ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1]';
ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2]';
ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2]';
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
[Xsm,Ysm] = pol2cart(th,rad);
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.86:0.01:0.93);
[Xmid,Ymid] = pol2cart(th,rad);

allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));


% Cell 1
figcells=figure(2);
% if max(ZBranch1)>0.5
alphaData=ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1)));
surf(Xcol,Ycol1,ZBranch1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
colormap(branchedColor)
freezeColors;
freezeColors(colorbar('Location','westoutside'));
clim([0,allmax])
shading interp
% end
hold on;
% if max(ZBund1)>0.5
alphaData=ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1)));
surf(Xcol,Ycol1,ZBund1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
colormap(bundledColor)
freezeColors;
freezeColors(jicolorbar);
clim([0,allmax])
shading interp
% end
view(2)
grid off
set(gca,'XTick',[], 'YTick', [])

% Cell 2
% if max(ZBranch2)>0.5
alphaData=ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2)));
surf(Xcol,Ycol2,ZBranch2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
colormap(branchedColor)
freezeColors;
freezeColors(colorbar('Location','westoutside'));
clim([0,allmax])
shading interp
% end
% if max(ZBund2)>0.5
alphaData=ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2)));
surf(Xcol,Ycol2,ZBund2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
colormap(bundledColor)
freezeColors;
freezeColors(jicolorbar);
clim([0,allmax])
shading interp
% end
view(2)
grid off
axis equal
set(gca,'XTick',[], 'YTick', [])
title(strcat(branchedColName, '=Branched, ', bundledColName, '=Bundled'))

% flipc2 = flip(boundC2);
% for i=1:length(boundC1)
%     plot3([Xcol(end,boundC1(i)) Xcol(end,flipc2(i))], [Ycol1(end,boundC1(i)) Ycol2(end,flipc2(i))],[allmax+1,allmax+1],'black')
% end

hold off;
box off;
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf,'color','w');


if ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
    figure(2)
    hold on;
    quiver(0,0,Xsm(dirIndexa1),Ysm(dirIndexa1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5);
    hold off;
end
if ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
    figure(2)
    hold on;
    quiver(0,-2*abs(max(max(Ycol2)))-gapsize,Xsm(dirIndexa2),Ysm(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
    hold off;
end


% Plot signal
figure(2)
[th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
[Xsig,Ysig] = pol2cart(th,rad);
if signal==1
    figure(2)
    hold on;
    scatter(Xsig(sigBound2),Ysig(sigBound2)-2*abs(max(max(Ycol2)))-gapsize,'black','.')
    hold off;
end

% Find peak for cell 1
N=length(a1);
smoothed1=zeros(1,N);
smoothAmount=20;
region=-smoothAmount:smoothAmount;
for i=1:N
    b=i-region;
    smoothed1(i) = sum(a1(b.*(b>0 & b<=N)+(b+N).*(b<=0)+(b-N).*(b>N)))/(2*smoothAmount+1);
end

% epsilon=0.01;                     % variance squared (default)
% x=linspace(0,L,N);
% a0=sqrt(pi)/(sqrt(2*pi*epsilon)); % area of normal distribution
% smoothed=zeros(length(x),N);
% for i=1:N
%     smoothed(:,i) = exp(-(x-a1(i)).^2/(2*epsilon))./(sqrt(2*pi*epsilon));
% end
% smoothed=exp(-(a1.^2)/(2*epsilon^2))./(sqrt(2*pi*epsilon));
% smoothed=smoothed/a0;
% smoothed=sum(smoothed,2);

% t = 0:.001:1-0.001;
% Fs = 1e3;
% a1New = detrend(a1,0);
% xdft = fft(a1);
% freq = 0:Fs/length(a1):Fs/2;
% xdft = xdft(1:length(a1)/2+1);
% plot(freq,abs(xdft));
% [~,I] = max(abs(xdft));

baseline=1;
peakIndices1=[];
peakIndex=0;
peakValue=0;
for i=1:N
    if smoothed1(i)>baseline
        if peakValue==0 || smoothed1(i)>peakValue
            peakIndex=i;
            peakValue=smoothed1(i);
        end
    elseif smoothed1(i)<baseline && peakIndex~=0
        peakIndices1=[peakIndices1,peakIndex];
        peakIndex=0;
        peakValue=0;
    end
end
if peakIndex~=0
    peakIndices1=[peakIndices1,peakIndex];
end

[~,maxind] = max(smoothed1(peakIndices1));
peakInd1 = peakIndices1(maxind);
if ~isempty(peakInd1)
    figure(2)
    % hold on
    % quiver(0,0,Xsm(peakInd1),Ysm(peakInd1),0,'color',[1 0 0],'LineWidth',2,'MaxHeadSize',0.5)
    % hold off;
end

% Find peak for cell 2
N=length(a2);
smoothed2=zeros(1,N);
region=-smoothAmount:smoothAmount;
for i=1:N
    b=i-region;
    smoothed2(i) = sum(a2(b.*(b>0 & b<=N)+(b+N).*(b<=0)+(b-N).*(b>N)))/(2*smoothAmount+1);
end

peakIndices2=[];
peakIndex=0;
peakValue=0;
for i=1:N
    if smoothed2(i)>baseline
        if peakValue==0 || smoothed2(i)>peakValue
            peakIndex=i;
            peakValue=smoothed2(i);
        end
    elseif smoothed2(i)<baseline && peakIndex~=0
        peakIndices2=[peakIndices2,peakIndex];
        peakIndex=0;
        peakValue=0;
    end
end
if peakIndex~=0
    peakIndices2=[peakIndices2,peakIndex];
end
[maxval,maxind] = max(smoothed2(peakIndices2));
peakInd2 = peakIndices2(maxind);
if ~isempty(peakInd2)
    figure(2)
    % hold on
    % quiver(0,-2,Xsm(peakInd2),Ysm(peakInd2),0,'color',[1 0 0],'LineWidth',2,'MaxHeadSize',0.5)
    % hold off;
end

figure(2)
ohf = findobj(gcf);
figaxes = findobj(ohf(1), 'Type', 'axes');
set(figaxes(1),'Fontsize',15)
set(figaxes(2),'Fontsize',14)
camroll(90)


figure(3)
gapsize=0.05;
clf
hold on;
th = (0:3.6:360)*pi/180;
Xvals=cos(th);
Yvals1=sin(th);
Yvals1(boundC1)=Yvals1(boundC1(1));
Yvals2=sin(th);
Yvals2(boundC2)=Yvals2(boundC2(1));
plot(Xvals,Yvals1,'black')
plot(Xvals,Yvals2-2*abs(max(Yvals2))-gapsize,'black')
YRac1=sin(posx1(1:NNx1(t),1)*2*pi/L);
YRac1(YRac1<Yvals1(boundC1(1)))=Yvals1(boundC1(1));
YRho1=sin(posy1(1:NNy1(t),1)*2*pi/L);
YRho1(YRho1<Yvals1(boundC1(1)))=Yvals1(boundC1(1));
YRac2=sin(posx2(1:NNx2(t),1)*2*pi/L)-gapsize;
YRac2(YRac2>Yvals2(boundC2(1)))=Yvals2(boundC2(1))-gapsize;
YRho2=sin(posy2(1:NNy2(t),1)*2*pi/L)-gapsize;
YRho2(YRho2>Yvals2(boundC2(1)))=Yvals2(boundC2(1))-gapsize;
scatter(cos(posx1(1:NNx1(t),1)*2*pi/L),YRac1,'MarkerEdgeColor',branchedColor(end,:),'linewidth',2)
scatter(cos(posx2(1:NNx2(t),1)*2*pi/L),YRac2-2*abs(max(Yvals2)),'MarkerEdgeColor',branchedColor(end,:),'linewidth',2)
scatter(cos(posy1(1:NNy1(t),1)*2*pi/L),YRho1,'MarkerEdgeColor',bundledColor(end,:),'linewidth',2)
scatter(cos(posy2(1:NNy2(t),1)*2*pi/L),YRho2-2*abs(max(Yvals2)),'MarkerEdgeColor',bundledColor(end,:),'linewidth',2)
hold off;
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf,'color','w')
axis([-1 1 -3 1])
axis equal

if ~isempty(dirIndexa2)
    figure(3)
    hold on;
    quiver(0,-2*abs(max(max(Ycol2)))-gapsize,Xsm(dirIndexa2),Ysm(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
    hold off;
end
if ~isempty(dirIndexa1)
    figure(3)
    hold on;
    quiver(0,0,Xsm(dirIndexa1),Ysm(dirIndexa1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5);
    hold off;
end

figure(3)
camroll(90)


% figure(3)
% subplot(1,2,1)
% hold on
% plot(Xa,a1);
% plot(Xa,smoothed1);
% hold off
% subplot(1,2,2)
% hold on
% plot(Xa,a2);
% plot(Xa,smoothed2);
% hold off


figcircle=figure(4);
clf
range=3;
% subplot(1,2,1)
plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L),'color','black','LineWidth',2)
hold on
alphaData=ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1)));
surf(Xcol,Ycol,ZBranch1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
% hold on
colormap(branchedColor)
% cb=colorbar('Location','westoutside');
freezeColors;
% freezeColors(cb);
% clim([0,allmax/2])
% cbpos=cb.Position;
% set(cb,'Position',[cbpos(1) cbpos(2) cbpos(3) cbpos(4)/2])
shading interp
alphaData=ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1)));
surf(Xcol,Ycol,ZBund1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
colormap(bundledColor)
freezeColors;
% jcb=jicolorbar;
% freezeColors(jcb);
% clim([0,allmax/2])
shading interp
view(2)
% plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L),'color','black')
% plot(2*cos(2*pi*Xa/L),2*sin(2*pi*Xa/L),'color','black')
if max(xC1)>=(range+2)
    racxvals1=(range-1)*xC1/max(xC1)+1;
    racyvals1=(range-1)*xC1/max(xC1)+1;
else
    racxvals1=(range-2)*xC1/max(xC1)+1;
    racyvals1=(range-1)*xC1/max(xC1)+1;
end
racxvals1=(racxvals1)'.*cos(2*pi*Xa/L);
racyvals1=(racyvals1)'.*sin(2*pi*Xa/L);
plot3(racxvals1,racyvals1,(allmax+1)*ones(1,length(racxvals1)),'color',...
    [branchedColor(end,:),1],'LineWidth',3)
plot3([racxvals1(end),racxvals1(1)],[racyvals1(end),racyvals1(1)],...
    [allmax+1,allmax+1],'color',[branchedColor(end,:),1],'LineWidth',3)
if max(yC1)>=(range+2)
    rhoxvals1=(range-1)*yC1/max(yC1)+1;
    rhoyvals1=(range-1)*yC1/max(yC1)+1;
else
    rhoxvals1=(range-2)*yC1/max(yC1)+1;
    rhoyvals1=(range-2)*yC1/max(yC1)+1;
end
rhoxvals1=(rhoxvals1)'.*cos(2*pi*Xa/L);
rhoyvals1=(rhoyvals1)'.*sin(2*pi*Xa/L);
plot3(rhoxvals1,rhoyvals1,(allmax+1)*ones(1,length(rhoxvals1)),'color',...
    [bundledColor(end,:),1],'LineWidth',3)
plot3([rhoxvals1(end),rhoxvals1(1)],[rhoyvals1(end),rhoyvals1(1)],...
    [allmax+1,allmax+1],'color',[bundledColor(end,:),1],'LineWidth',3)
% xlim([-3,3])
% ylim([-3,3])
% axis square
hold off

%cell 2
hold on
if adjacent==0
    plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L)-(2*range),'color','black','LineWidth',2)
else
    plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L)-2-(range-1),'color','black','LineWidth',2)
end
alphaData=ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2)));
if adjacent==0
    surf(Xcol,Ycol-(2*range),ZBranch2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
else
    surf(Xcol,Ycol-2-(range-1),ZBranch2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
end
% hold on
colormap(branchedColor)
cb=colorbar('Location','eastoutside');
freezeColors;
freezeColors(cb);
clim([0,allmax/2])
cbpos=cb.Position;
% set(cb,'Position',[cbpos(1)+2*cbpos(3),cbpos(2),cbpos(3),cbpos(4)/2])
set(cb,'Position',[0.9062    0.1097    0.0235    0.4077])
set(cb,'TickLabels',{});
cbpos=cb.Position;
shading interp
alphaData=ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2)));
if adjacent==0
    surf(Xcol,Ycol-(2*range),ZBund2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
else
    surf(Xcol,Ycol-2-(range-1),ZBund2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
end
colormap(bundledColor)
freezeColors;
jcb=jicolorbar;
freezeColors(jcb);
clim([0,allmax/2])
jcbpos=jcb.Position;
set(jcb,'Position',[cbpos(1)+cbpos(3),cbpos(2),cbpos(3),cbpos(4)])
% set(cb,'Position',[cbpos(1),cbpos(2),cbpos(3),jcbpos(4)])
shading interp
view(2)
% plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L),'color','black')
% plot(2*cos(2*pi*Xa/L),2*sin(2*pi*Xa/L),'color','black')
if max(xC2)>=(range+2)
    racxvals2=(range-1)*xC2/max(xC2)+1;
    racyvals2=(range-1)*xC2/max(xC2)+1;
else
    racxvals2=(range-2)*xC2/max(xC2)+1;
    racyvals2=(range-2)*xC2/max(xC2)+1;
end
racxvals2=(racxvals2)'.*cos(2*pi*Xa/L);
racyvals2=(racyvals2)'.*sin(2*pi*Xa/L);
if adjacent==0
    plot3(racxvals2,racyvals2-(2*range),(allmax+1)*ones(1,length(racxvals2)),'color',...
        branchedColor(end,:),'LineWidth',3)
    plot3([racxvals2(end),racxvals2(1)],[racyvals2(end),racyvals2(1)]-(2*range),...
        [allmax+1,allmax+1],'color',branchedColor(end,:),'LineWidth',3)
else
    plot3(racxvals2,racyvals2-2-(range-1),(allmax+1)*ones(1,length(racxvals2)),'color',...
        [branchedColor(end,:),0.5],'LineWidth',3)
    plot3([racxvals2(end),racxvals2(1)],[racyvals2(end),racyvals2(1)]-2-(range-1),...
        [allmax+1,allmax+1],'color',[branchedColor(end,:),0.5],'LineWidth',3)
end
if max(yC2)>=(range+2)
    rhoxvals2=(range-1)*yC2/max(yC2)+1;
    rhoyvals2=(range-1)*yC2/max(yC2)+1;
else
    rhoxvals2=(range-2)*yC2/max(yC2)+1;
    rhoyvals2=(range-2)*yC2/max(yC2)+1;
end
rhoxvals2=(rhoxvals2)'.*cos(2*pi*Xa/L);
rhoyvals2=(rhoyvals2)'.*sin(2*pi*Xa/L);
if adjacent==0
    plot3(rhoxvals2,rhoyvals2-(2*range),(allmax+1)*ones(1,length(rhoxvals2)),'color',...
        bundledColor(end,:),'LineWidth',3)
    plot3([rhoxvals2(end),rhoxvals2(1)],[rhoyvals2(end),rhoyvals2(1)]-2*range,...
        [allmax+1,allmax+1],'color',bundledColor(end,:),'LineWidth',3)
else
    plot3(rhoxvals2,rhoyvals2-2-(range-1),(allmax+1)*ones(1,length(rhoxvals1)),'color',...
        [bundledColor(end,:),0.5],'LineWidth',3)
    plot3([rhoxvals2(end),rhoxvals2(1)],[rhoyvals2(end),rhoyvals2(1)]-2-(range-1),...
        [allmax+1,allmax+1],'color',[bundledColor(end,:),0.5],'LineWidth',3)
end
if adjacent==0
    xlim([-3,3])
    ylim([-9,3])
    set(gca,'plotBoxAspectRatio',[1 2 1]);
else
    % xlim([-3,3])
    % ylim([-5,3])
    % set(gca,'plotBoxAspectRatio',[6 8 1]);
    xlim([-3,3])
    ylim([-7,3])
    set(gca,'plotBoxAspectRatio',[6 10 1]);
end
% axis square
hold off

if ~isempty(dirIndexa1)
    figure(4)
    hold on;
    quiver(0,0,Xsm(dirIndexa1),Ysm(dirIndexa1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7);
    hold off;
end
if ~isempty(dirIndexa2)
    figure(4)
    hold on;
    if adjacent==0
        quiver(0,-2*range,Xsm(dirIndexa2),Ysm(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7)
    else
        quiver(0,-2-(range-1),Xsm(dirIndexa2),Ysm(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7)
    end
    hold off;
end

grid off
set(gca,'XTick',[],'YTick',[])
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf,'color','w');
set(gcf,'Position',[209   561   682   474])
ohf = findobj(gcf);
figaxes = findobj(ohf(1), 'Type', 'axes');
set(figaxes(1),'Fontsize',15)
set(figaxes(2),'Fontsize',14)
camroll(90)




end_time=2500;
Tend=25;
figure(5)
clf
% subplot(2,1,1)
% for j=1:2:max(max([NNx1,NNy1]))
[X,Y,Z] = cylinder(0.99);
surf(Z*Tend*end_time/Nt,Y,X,ones(size(Z,1),size(Z,2),3))
hold on
surf(Z*Tend*end_time/Nt,Y-2,X,ones(size(Z,1),size(Z,2),3))
alpha 0.5
shading interp
% xlabel('time')
% ylabel('y')
% zlabel('z')
hold on
for j=1:2:200
    scatter3(linspace(0,Tend*end_time/Nt,end_time),sin(posx1(j,1:end_time)*2*pi/10),cos(posx1(j,1:end_time)*2*pi/10),...
        1,'MarkerEdgeColor',branchedColor(end,:),'MarkerFaceColor',branchedColor(end,:))
    hold on;
    scatter3(linspace(0,Tend*end_time/Nt,end_time),sin(posy1(j,1:end_time)*2*pi/10),cos(posy1(j,1:end_time)*2*pi/10),...
        1,'MarkerEdgeColor',bundledColor(end,:),'MarkerFaceColor',bundledColor(end,:))
    box on;
    % set(gca,'Color','k','fontsize',20,'fontname','times');
    % pbaspect([3 2 1]);
    % set(gcf,'color','w');
    % title('Cell 1')
end
hold off
xlabel('Time')
% set(gca,'CameraPosition',[-9.629901052021953,-7.569863519089166,5.073854444144635])

% subplot(2,1,2)
for j=1:2:200
    hold on;
    scatter3(linspace(0,Tend*end_time/Nt,end_time),sin(posx2(j,1:end_time)*2*pi/10)-2,cos(posx2(j,1:end_time)*2*pi/10),...
        1,'MarkerEdgeColor',branchedColor(end,:),'MarkerFaceColor',branchedColor(end,:))
    % hold on;
    scatter3(linspace(0,Tend*end_time/Nt,end_time),sin(posy2(j,1:end_time)*2*pi/10)-2,cos(posy2(j,1:end_time)*2*pi/10), ...
        1,'MarkerEdgeColor',bundledColor(end,:),'MarkerFaceColor',bundledColor(end,:))
    box on;
    % set(gca,'Color','k','fontsize',20,'fontname','times');
    % pbaspect([3 2 1]);
    % set(gcf,'color','w');
    % title('Cell 2')
end
hold off
set(gca,'Color','k','fontsize',20,'fontname','times');
xlim([0,Tend*end_time/Nt])
pbaspect([Tend*end_time/Nt (Tend*end_time/Nt)*2/3 (Tend*end_time/Nt)/3]);
set(gcf,'color','w');
xlabel('Time')
ylabel('y')
zlabel('z')
% set(gca,'CameraPosition',[-9.629901052021953,-7.569863519089166,5.073854444144635])



% figure(6)
% clf
% ccx = [0 0 255]/256.*ones(Nt,1);     % blue
% ccy = [255 219 88]/256.*ones(Nt,1);  % mustard yellow
% time = linspace(0,Tend,Nt);
% % for j=1:2:max(max([NNx1,NNy1]))
% for j=1:200
%     hold on;
%     scatter(linspace(0,Tend*end_time/Nt,end_time),posx1(j,1:end_time),1,branchedColor(end,:));
%     scatter(linspace(0,Tend*end_time/Nt,end_time),posy1(j,1:end_time),1,bundledColor(end,:));
%     box on;
%     set(gca,'Color','k','fontsize',20,'fontname','times');
%     pbaspect([3 1 1]);
%     set(gcf,'color','w');
%     title('Cell 1')
% end
% hold off

