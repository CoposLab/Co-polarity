% 1 cell

set(0,'DefaultFigureVisible','off')

filename = './singlecellresults/RhoStopsT=10/20RhoLeft_NoBind_NoUnbind_NoFB_NoDiffusion';
maxnum=100;

pol_time_10_counter=0;
pol_time_end_counter=0;
pol_both_counter=0;


tic
for i=1:maxnum
    pol_10=0;

    close all
    fig = openfig(strcat(filename,'_t10_Scatter_',string(i),'.fig'));
    ohf = findobj(gcf);
    figaxes = findobj(ohf(1), 'Type', 'axes');

    cell1 = findobj(figaxes(1),'Type','line');
    xBundled1 = get(cell1(3),'XData');
    yBundled1 = get(cell1(3),'YData');
    xBranched1 = get(cell1(4),'XData');
    yBranched1 = get(cell1(4),'YData');


    % Find median for cell 1
    a1New = yBranched1;
    a1New(a1New<1)=0;
    if (a1New(1)~=0 && a1New(end)~=0)
        zeroInda1_1=find(a1New==0,1,'first');
        zeroInda2_1=find(a1New==0,1,'last');
        dirIndexa1=ceil((zeroInda1_1+zeroInda2_1)/2) - 50;
    else
        inda1_1=find(a1New~=0,1,'first');
        inda2_1=find(a1New~=0,1,'last');
        dirIndexa1=ceil((inda1_1+inda2_1)/2);
    end
    b1New = yBundled1;
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
    
    if dirIndexa1<1
        dirIndexa1=dirIndexa1+101;
    end
    if dirIndexb1<1
        dirIndexb1=dirIndexb1+101;
    end
    

    if ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
        pol_10=1;
        pol_time_10_counter=pol_time_10_counter+1;
    end


    close all
    fig = openfig(strcat(filename,'Scatter_',string(i),'.fig'));
    ohf = findobj(gcf);
    figaxes = findobj(ohf(1), 'Type', 'axes');

    cell1 = findobj(figaxes(1),'Type','line');
    xBundled1 = get(cell1(3),'XData');
    yBundled1 = get(cell1(3),'YData');
    xBranched1 = get(cell1(4),'XData');
    yBranched1 = get(cell1(4),'YData');


    a1=yBranched1;
    b1=yBundled1;

    [th,rad] = meshgrid((0:3.6:360)*pi/180,1);

    % Find median for cell 1
    a1New = yBranched1;
    a1New(a1New<1)=0;
    if (a1New(1)~=0 && a1New(end)~=0)
        zeroInda1_1=find(a1New==0,1,'first');
        zeroInda2_1=find(a1New==0,1,'last');
        dirIndexa1=ceil((zeroInda1_1+zeroInda2_1)/2) - 50;
    else
        inda1_1=find(a1New~=0,1,'first');
        inda2_1=find(a1New~=0,1,'last');
        dirIndexa1=ceil((inda1_1+inda2_1)/2);
    end
    b1New = yBundled1;
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
    
    if dirIndexa1<1
        dirIndexa1=dirIndexa1+101;
    end
    if dirIndexb1<1
        dirIndexb1=dirIndexb1+101;
    end

    if ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
        pol_time_end_counter=pol_time_end_counter+1;
        if pol_10==1
            pol_both_counter=pol_both_counter+1;
        end
    end
    

end
toc
sprintf('Number Polarized at t=10: %d\nNumber Polarized at t=25: %d\nNumber Polarized at t=10 and t=25: %d',...
    pol_time_10_counter,pol_time_end_counter,pol_both_counter)



