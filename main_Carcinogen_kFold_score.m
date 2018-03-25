%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform the biomarker selection for the carcinogenicity prediction
%   This study is cinducted using the yeast proteomic assay
%   Developed by: Sheikh M. Rahman
%   Date: March, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ('data/peli_carcinogen.mat'); % load carcinogenicity data
addpath('lib/');  % add necessary script library saved in the "lib" folder

%% Get the score of the biomarkers using three score Criteria
% The score methods are t-stat, MRMR-TCD, and MRMR-mrmrTCQ
% A 10-fold cross validation is used

[score, featureIdsorted] = featureRank_kfold_allThreeMethods(peli_carc, classlabel, 10);
save('results/kfold_score_three_rankingMethod_carc.mat','score','featureIdsorted');

%% plot three scores in one figure; figure 1 of the dissertation chapter
h3 = figure;
set(h3, 'PaperUnits','inches','Units','inches','Position',[2 2 9 3.5], ...
            'PaperSize',[9 3.5], 'PaperPosition',[0 0 9 3.5]);
% t-stat score plot
subplot('position',[0.055 0.23 0.27 0.75])
drawFigures_allThreeScore (score.tstatScore, featureIdsorted.tstatScore, geneName_N, 't-stat Score');
annotation('textbox','position',[0.04, 0.01, 0.27, 0.05],'String','a) Scoring Criteria: t-stat',...
       'linestyle', 'none','fontname','Arial','fontsize',12, 'fontWeight','bold',...
       'margin',0,'HorizontalAlignment','center','VerticalAlignment','bottom');

% MRMR-TCD score plot
subplot('position',[0.387 0.23 0.27 0.75])
drawFigures_allThreeScore (score.mrmrTCD, featureIdsorted.mrmrTCD, geneName_N, 'MRMR-TCD Score');
annotation('textbox','position',[0.35, 0.01, 0.33, 0.05],'String','b) Scoring Criteria: MRMR-TCD',...
       'linestyle', 'none','fontname','Arial','fontsize',12, 'fontWeight','bold',...
       'margin',0,'HorizontalAlignment','center','VerticalAlignment','bottom');
set(gca,'ylim',[-0.05, 1]);

% MRMR-TCQ score plot
subplot('position',[0.722 0.23 0.27 0.75])
drawFigures_allThreeScore (score.mrmrTCQ, featureIdsorted.mrmrTCQ, geneName_N, 'MRMR-TCQ Score');
annotation('textbox','position',[0.69, 0.01, 0.33, 0.05],'String','c) Scoring Criteria: MRMR-TCQ',...
       'linestyle', 'none','fontname','Arial','fontsize',12, 'fontWeight','bold',...
       'margin',0,'HorizontalAlignment','center','VerticalAlignment','bottom');

% save the figure as pdf file
print (h3, '-dpdf', '-r300', 'results/biomarkerScore_carc_10fold.pdf');

%% Get the prediction models for varying number of biomarekrs and calculate the accuracy

% prediction model for top five biomarkers

id = featureIdsorted.mrmrTCQ (1:5);
classModel_top5 = fitcsvm(peli_carc(:,id), classlabel, 'kfold', 10, 'kernelfunction', 'linear');
[accuracy_train(1), sensitivity_train(1), specificity_train(1), AUC_train(1)] = ...
        performanceCriteria_CV(classModel_top5, classlabel);

% prediction model for all the biomarkers
id = featureIdsorted.mrmrTCQ;
classModel_all = fitcsvm(peli_carc(:,id), classlabel, 'kfold', 10, 'kernelfunction', 'linear');
[accuracy_train(2), sensitivity_train(2), specificity_train(2), AUC_train(2)] = ...
        performanceCriteria_CV(classModel_all, classlabel);

% export the AUC and accuracy parameters (values for table 1)
writetable(table(AUC_train', accuracy_train', sensitivity_train', specificity_train',...
        'VariableNames', {'AUC_train','accuracy_train', 'sensitivity_train', 'specificity_train'},...
        'RowNames',{'Top Five','All'}),...
        'results/carc_accuracy_top5_topall.csv', 'WriteRowNames',1);

% save the models as .mat file
save ('results/predictionModels_carc.mat','classModel_top5','classModel_all','classlabel');

%% Plot the heatmap of protein PELI with the biomarker sorted based on MRMR-TCQ score (Figure 2)
% prepare data
peliGene = peli_carc;
allChemName = AllchemName_Carc;
idSorted_carc = featureIdsorted.mrmrTCQ;
id_temp = strfind(classlabel,'positive');
id_positive = ~cellfun(@isempty,id_temp);
id_negative = cellfun(@isempty,id_temp);
peliGene_ord = [peliGene(id_positive,:); peliGene(id_negative,:)];
chemName_ord = [allChemName(id_positive); allChemName(id_negative)];

[numSample, numGene] =size(peliGene);
numPos = sum(id_positive); numNeg = sum(id_negative);

% plot the figure
fh = figure;
set(fh, 'PaperUnits','inches', 'PaperSize',[7 8], 'PaperPosition',[0 0 7 8],...
    'Units','inches', 'Position',[0 0 7 8]);
axes('position',[0.34 0.1 0.645 0.88]);
imagesc(peliGene_ord(:,idSorted_carc));
d=jet(36); colormap(d([4,20:36],:));

set(gca,'xtick',1:1:numGene, 'xticklabel',geneName_N(idSorted_carc),'xticklabelrotation',90,...
    'ytick',0.5:6:120.5, 'yticklabel','', 'tickdir','out', 'fontname','arial','fontsize', 10);
set(gca, 'ygrid','on','layer', 'top','GridColor','w','gridalpha',0.7,...
    'gridlinestyle','--', 'linewidth',1);
set(gca, 'ticklength', [0.005 0.1]);

t1 = text(repmat(-0.2,numPos/6,1), 3.5:6:numPos, chemName_ord(1:6:numPos),'color','red',...
    'fontname', 'arial', 'fontsize', 12);
set(t1, 'HorizontalAlignment','right','VerticalAlignment','middle');
t2 = text(repmat(-0.2,numNeg/6,1), numPos+3.5:6:numSample, chemName_ord(numPos+1:6:end),...
    'fontname', 'arial', 'fontsize', 12);
set(t2, 'HorizontalAlignment','right','VerticalAlignment','middle');

c = colorbar('northoutside','fontname','arial','fontsize',10);
c.Label.String = 'PELI'; c.Label.FontName = 'arial'; c.Label.FontSize = 12;
c.Ticks = [1,1.5,2:1:10];
xlabel('Proteins', 'fontname', 'arial', 'fontsize', 14);
ha = axes('position',[0.0118 0.33 0.05 0.55]);
drawbrace([0 2],[0 20],0.2,'color','r');ylim([0 20]);ha.Visible = 'off';
text(-1.9, 9.2, 'Positive','rotation',90,...
    'fontname','Arial','fontsize',15, 'fontWeight','bold');
hb = axes('position',[0.0118 0.082 0.05 0.29]);
drawbrace([0 2],[0 20],0.2,'color','k');ylim([0 20]);hb.Visible = 'off';
text(-1.9, 7.5, 'Negative','rotation',90,...
    'fontname','Arial','fontsize',15, 'fontWeight','bold');

% save the figure as heatmap
print(fh, '-dpdf','-r300','results/heatmap_carcinogen_PELI_FCQ.pdf')
