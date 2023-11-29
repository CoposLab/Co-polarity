% 4 cells line

set(0,'DefaultFigureVisible','off')

filename = './simulation_results/results_line/alternateracuprhoup/1000RacOn1000RhoOn';
maxnum=100;

angletolerance=pi/4; % tolerance for yes
anglelf=pi/4; % tolerance for leader follower

countallleft=0;
countallright=0;
countyes=0;
countlf=0;
countcw=0;
countccw=0;

tic
for i=1:maxnum
    close all
    fig = openfig(strcat(filename,'Scatter_',string(i),'.fig'));
    ohf = findobj(gcf);
    figaxes = findobj(ohf(1), 'Type', 'axes');

    cell1 = findobj(figaxes(4),'Type','line');
    xBundled1 = get(cell1(3),'XData');
    yBundled1 = get(cell1(3),'YData');
    xBranched1 = get(cell1(4),'XData');
    yBranched1 = get(cell1(4),'YData');

    cell2 = findobj(figaxes(3),'Type','line');
    xBundled2 = get(cell2(3),'XData');
    yBundled2 = get(cell2(3),'YData');
    xBranched2 = get(cell2(4),'XData');
    yBranched2 = get(cell2(4),'YData');

    cell3 = findobj(figaxes(2),'Type','line');
    xBundled3 = get(cell3(3),'XData');
    yBundled3 = get(cell3(3),'YData');
    xBranched3 = get(cell3(4),'XData');
    yBranched3 = get(cell3(4),'YData');

    cell4 = findobj(figaxes(1),'Type','line');
    xBundled4 = get(cell4(3),'XData');
    yBundled4 = get(cell4(3),'YData');
    xBranched4 = get(cell4(4),'XData');
    yBranched4 = get(cell4(4),'YData');

    % figure(2)
    % subplot(1,4,1)
    % hold on
    % plot(xBundled1,yBundled1)
    % plot(xBranched1,yBranched1)
    % hold off
    % subplot(1,4,2)
    % hold on
    % plot(xBundled2,yBundled2)
    % plot(xBranched2,yBranched2)
    % hold off
    % subplot(1,4,3)
    % hold on
    % plot(xBundled3,yBundled3)
    % plot(xBranched3,yBranched3)
    % hold off
    % subplot(1,4,4)
    % hold on
    % plot(xBundled4,yBundled4)
    % plot(xBranched4,yBranched4)
    % hold off

    a1=yBranched1;
    b1=yBundled1;
    a2=yBranched2;
    b2=yBundled2;
    a3=yBranched3;
    b3=yBundled3;
    a4=yBranched4;
    b4=yBundled4;

    [th,rad] = meshgrid((0:3.6:360)*pi/180,1);

    % Find median for cell 1
    a1New = a1;
    a1New(a1New<1)=0;
    if (a1New(1)~=0 && a1New(length(a1New))~=0)
        zeroInd1=find(a1New==0,1,'first');
        zeroInd2=find(a1New==0,1,'last');
        dirIndex1=ceil((zeroInd1+zeroInd2)/2) - 50;
    else
        ind1=find(a1New~=0,1,'first');
        ind2=find(a1New~=0,1,'last');
        dirIndex1=ceil((ind1+ind2)/2);
    end
    if dirIndex1<1
        dirIndex1=dirIndex1+101;
    end

    % Find median for cell 2
    a2New = a2;
    a2New(a2New<1)=0;
    if (a2New(1)~=0 && a2New(length(a2New))~=0)
        zeroInd1=find(a2New==0,1,'first');
        zeroInd2=find(a2New==0,1,'last');
        dirIndex2=ceil((zeroInd1+zeroInd2)/2) - 50;
    else
        ind1=find(a2New~=0,1,'first');
        ind2=find(a2New~=0,1,'last');
        dirIndex2=ceil((ind1+ind2)/2);
    end
    if dirIndex2<1
        dirIndex2=dirIndex2+101;
    end

    % Find median for cell 3
    a3New = a3;
    a3New(a3New<1)=0;
    if (a3New(1)~=0 && a3New(length(a3New))~=0)
        zeroInd1=find(a3New==0,1,'first');
        zeroInd2=find(a3New==0,1,'last');
        dirIndex3=ceil((zeroInd1+zeroInd2)/2) - 50;
    else
        ind1=find(a3New~=0,1,'first');
        ind2=find(a3New~=0,1,'last');
        dirIndex3=ceil((ind1+ind2)/2);
    end
    if dirIndex3<1
        dirIndex3=dirIndex3+101;
    end

    % Find median for cell 4
    a4New = a4;
    a4New(a4New<1)=0;
    if (a4New(1)~=0 && a4New(length(a4New))~=0)
        zeroInd1=find(a4New==0,1,'first');
        zeroInd2=find(a4New==0,1,'last');
        dirIndex4=ceil((zeroInd1+zeroInd2)/2) - 50;
    else
        ind1=find(a4New~=0,1,'first');
        ind2=find(a4New~=0,1,'last');
        dirIndex4=ceil((ind1+ind2)/2);
    end
    if dirIndex4<1
        dirIndex4=dirIndex4+101;
    end

    if ~isempty(dirIndex1) && ~isempty(dirIndex2) && ~isempty(dirIndex3) && ~isempty(dirIndex4)
        medang1 = th(1,dirIndex1);
        medang2 = th(1,dirIndex2);
        medang3 = th(1,dirIndex3);
        medang4 = th(1,dirIndex4);

        if medang1>=0 && medang1<=pi && medang2>=0 && medang2<=pi ...
                && medang3>=0 && medang3<=pi && medang4>=0 && medang4<=pi
            countallleft=countallleft+1;
        elseif medang1>=pi && medang1<=2*pi && medang2>=pi && medang2<=2*pi ...
                && medang3>=pi && medang3<=2*pi && medang4>=pi && medang4<=2*pi
            countallright=countallright+1;
        end

        angdiff1=max([min(abs(medang1-medang2),2*pi-abs(medang1-medang2)),... %max angle diff between cell 1 and other cells
            min(abs(medang1-medang3),2*pi-abs(medang1-medang3)),...
            min(abs(medang1-medang4),2*pi-abs(medang1-medang4))]);
        angdiff2=max([min(abs(medang2-medang3),2*pi-abs(medang2-medang3)),... %max angle diff between c2,c3 and c2,c4
            min(abs(medang2-medang4),2*pi-abs(medang2-medang4))]);
        angdiff3=min(abs(medang3-medang4),2*pi-abs(medang3-medang4)); %angle diff between c3,c4

        if max([angdiff1,angdiff2,angdiff3])<=angletolerance
            countyes=countyes+1;
        end

        epsilon=0.01*2*pi;
        if xor(abs(medang2-3*pi/2)<anglelf, abs(medang3-pi/2)<anglelf) %cell 2 -> cell 3 or cell 2 <- cell 3
            if abs(medang1-3*pi/2)<anglelf && abs(medang2-3*pi/2)<anglelf && abs(medang3-3*pi/2)<anglelf % all pointing to cell 4
                countlf=countlf+1;
                %check CW vs CCW
                if (medang4>=0 && medang4<(pi/2)-epsilon) || (medang4>(3*pi/2)+epsilon && medang4<=2*pi)
                    countccw=countccw+1;
                elseif medang4>(pi/2)+epsilon && medang4<(3*pi/2)-epsilon
                    countcw=countcw+1;
                end
            elseif abs(medang2-pi/2)<anglelf && abs(medang3-pi/2)<anglelf && abs(medang4-pi/2)<anglelf % all pointing to cell 1
                countlf=countlf+1;
                %check CW vs CCW
                if (medang1>=0 && medang1<(pi/2)-epsilon) || (medang1>(3*pi/2)+epsilon && medang1<=2*pi)
                    countcw=countcw+1;
                elseif medang1>(pi/2)+epsilon && medang1<(3*pi/2)-epsilon
                    countccw=countccw+1;
                end
            end
        end
    end

end
toc

sprintf('All same direction: %d\nYes: %d\nLeader/Follower: %d\nCW: %d, CCW: %d',...
    countallleft+countallright, countyes, countlf, countcw, countccw)