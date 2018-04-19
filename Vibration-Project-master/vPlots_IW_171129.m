colors=jet(10);
C=1;

for cells=1:size(vPlots.plotData,2)
    figure(C);clf;%set(gca,'colororder', jet(10));
    colors=jet(10);
    hold on
    grid minor
    for k = 1:10
%        plot(vPlots(cells).time,vPlots(cells).ymeans(:,k),'color',colors(k,:))
        plot3(vPlots.plotData(cells).time,repmat(k,length(vPlots.plotData(cells).ymeans(:,k))),vPlots.plotData(cells).ymeans(:,k),'color',colors(k,:))
    end
    title(['#' num2str(cells) ': AH0800 C2, contact at ' num2str(vPlots.plotData(cells).dist) ' mm from base'],'FontSize',10)
    xlabel('time (ms)','FontSize',10)
    ylabel('Whisker section number','FontSize',10)
%             ylabel('mean y (mm), vibrational component','FontSize',10);
    zlabel('mean y (mm), vibrational component','FontSize',10)
    view(0,0)
    movegui east
    axis tight
    hold off
    
    figure(C+1);clf;hold on
    set(gca,'colororder', hsv(20));
    axis ij

    aidx = cell(1); %initializes aidx
    aidxp=cell(1);
    
    for i = 1:length(vPlots.plotData(cells).rawX)
        for k = 1:10
            plot(vPlots.plotData(cells).rawX{i}(vPlots.plotData(cells).aidx{i,k}),vPlots.plotData(cells).rawY{i}(vPlots.plotData(cells).aidx{i,k}),'k') %plots raw whisker data in each section aidx
            plot(vPlots.plotData(cells).polyValX{i}(vPlots.plotData(cells).aidxp{i,k}),vPlots.plotData(cells).polyValY{i}(vPlots.plotData(cells).aidxp{i,k}),'m') %plots polynomial fit to whisker data in each section aidx
        end
    end
    title(['#' num2str(cells) ': AH0800 C2, contact at ' num2str(vPlots.plotData(cells).dist) ' mm from base'],'FontSize',10)
    movegui west
    hold off
    axis tight
    
    C=C+2;
  %  print(gcf,'-dpng', [num2str(vPlots(cells).dist) '.png'])
end
figure(C+1)
imagesc(vPlots.whiskerImage.base)
colormap(gray)
title('Whisker base - follicle to 1/2 length')

figure(C+2)
imagesc(vPlots.whiskerImage.middle)
colormap(gray)
title('Whisker middle - 1/2 length to 2/3 length')

figure(C+3)
imagesc(vPlots.whiskerImage.tip)
colormap(gray)
title('Whisker tip - 2/3 length to end of whisker')