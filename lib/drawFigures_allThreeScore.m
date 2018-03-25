function [] = drawFigures_allThreeScore (sortedLoading,  IDX,  geneName,  ylbl)
    scatter(1:numel(sortedLoading),sortedLoading, 15, 'bo','filled'); %scatterplot with grouping
      box on;
    set(gca,'XTick',1:1:numel(sortedLoading),'XTickLabel',geneName(IDX),...
        'XTickLabelRotation',90,'fontname','Arial','fontsize',5.3,'TickDir', 'out');
    hy = get(gca); hy.YAxis.FontSize = 10; % set the font size of the x axis labe fontsize only
    xlabel('Proteins','fontsize',13); ylabel(ylbl,'fontsize',12);    
end