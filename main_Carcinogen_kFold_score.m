%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform the gene selection from the genotoxicity 
%   Developed by: Sheikh M. Rahman
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ('../peli_carcinogen.mat');
addpath('../../lib/');

%% rank the genes based on ttest
rng(0, 'twister');
[score_ttest, idSorted_ttest] = featureRankTtest_kfold(peli_carc, classlabel,10);
drawFigures(score_ttest,  idSorted_ttest, 'ttest_kfold', geneName_N, 'PosNeg', 'T-stat Score');   
nFeature = numel(geneName_N); %krnl = 'gaussian'; 
krnl = 'linear';
[trainPerformVls, nFeatVector] = clssificationError_CV(peli_carc, classlabel, ...
            nFeature, idSorted_ttest, 'SVM', 10, 'ttest', krnl);
save('Result/featuresSelectedFromTTest_SVM_kfold.mat','score_ttest','idSorted_ttest',...
            'nFeatVector', 'trainPerformVls');
plotAccuracy_CV (trainPerformVls, nFeatVector, 'ttest_SVM_10Fold');
[trainPerformVls.AUC trainPerformVls.accuracy]
%% MRMR FCD

nFeature = numel(geneName_N); type = 'FCD';
[score_mRmR_FCD, featureIdsorted_mRmR_FCD] = maxRelevanceMinRedundancy_kfold (peli_carc, classlabel, 10, type);

drawFigures(score_mRmR_FCD,  featureIdsorted_mRmR_FCD, strcat('MRMR_', type, '_kfold'), geneName_N, 'PosNeg', strcat('MRMR',{' '},type,' Score'));   
[trainingStruct, nFeatVec] = clssificationError_CV(peli_carc, classlabel, nFeature,...
        featureIdsorted_mRmR_FCD,'SVM',10, strcat('MRMR_',type),'linear');
save('Result/featuresSelectedFromMRMR_FCD_SVM_kfold_v2.mat','score_mRmR_FCD','featureIdsorted_mRmR_FCD',...
  'nFeatVec', 'trainingStruct');
plotAccuracy_CV(trainingStruct, nFeatVec, 'MRMR_FCD_SVM_10Fold');
[trainingStruct.AUC trainingStruct.accuracy]
%AUC_mrmr_FCD = trainingStruct.AUC;

%% MRMR FCQ

nFeature = numel(geneName_N); type = 'FCQ';
rng(1);
[score_mRmR_FCQ, featureIdsorted_mRmR_FCQ] = maxRelevanceMinRedundancy_kfold (peli_carc, classlabel, 10, type);

drawFigures(score_mRmR_FCQ,  featureIdsorted_mRmR_FCQ, strcat('MRMR_', type, '_kfold'), geneName_N, 'PosNeg', strcat('MRMR',{' '},type,' Score'));   
[trainingStruct, nFeatVec] = clssificationError_CV(peli_carc, classlabel, nFeature,...
        featureIdsorted_mRmR_FCQ,'SVM',10, strcat('MRMR_',type),'linear');
save('Result/featuresSelectedFromMRMR_FCQ_SVM_kfold_v2.mat','score_mRmR_FCQ','featureIdsorted_mRmR_FCQ',...
  'nFeatVec', 'trainingStruct');
plotAccuracy_CV(trainingStruct, nFeatVec, 'MRMR_FCQ_SVM_10Fold');
[trainingStruct.AUC trainingStruct.accuracy]
%save Result/AUC.mat AUC_mrmr_FCD AUC_mrmr_FCQ AUC_ttest;
%%
% plot AUC
h4 =figure;
set(h4, 'PaperUnits','inches', 'PaperSize',[5 3.5], 'PaperPosition',[0 0 5 3.5]);
plot(AUC_ttest, 'm-o','linewidth', 1); 
hold on; box on;
plot(AUC_mrmr_FCD, 'b-v','linewidth', 1); 
plot(AUC_mrmr_FCQ, 'r-s','linewidth', 1); 

set(gca,'XTick',1:1:numel(AUC_ttest),'XTickLabel',nFeatVec,...
    'fontname', 'Arial', 'fontsize', 12, 'xlim', [0, numel(AUC_ttest)+1]);
xlabel('No. of Biomarkers','fontsize',14); ylabel('AUC','fontsize',14);
ylim([0.5 1]);
legend({'T-stat','MRMR-TCD','MRMR-TCQ'},...
    'location','SouthEast','FontSize',12,'box','off');

if ~exist('Result/', 'dir'), mkdir('Result/'); end
print(h4, '-dtiff', '-r300', strcat('Result/AUC.tiff'));
print(h4, '-dpdf', '-r300', strcat('Result/AUC.pdf'));
%% for biomarker 1- 5
id = featureIdsorted_mRmR_FCQ (1:4);
rng(0);
classModel = fitcsvm(peli_carc(:,id), classlabel, 'kfold', 10);

[accuracy_train(1), sensitivity_train(1), specificity_train(1), AUC_train(1)] = ...
        performanceCriteria_CV(classModel, classlabel);
    
 id = featureIdsorted_mRmR_FCQ (1:5);
 rng(1);
classModel = fitcsvm(peli_carc(:,id), classlabel, 'kfold', 10);

[accuracy_train(2), sensitivity_train(2), specificity_train(2), AUC_train(2)] = ...
        performanceCriteria_CV(classModel, classlabel);
    
id = featureIdsorted_mRmR_FCQ;
rng(2);
classModel = fitcsvm(peli_carc(:,id), classlabel, 'kfold', 10);
[accuracy_train(3), sensitivity_train(3), specificity_train(3), AUC_train(3)] = ...
        performanceCriteria_CV(classModel, classlabel);

  fh = figure;
  set(fh, 'paperunits', 'inches', 'papersize', [7, 3.5], 'paperposition',[0 0 7 3.5]);
subplot(1,3,1)
bar([accuracy_train']*100, 'grouped');
colormap(cool); xlim([0 4]); ylim([0 100]);
set (gca,'xticklabel',{'Top Four','Top Five','All'},'xticklabelrotation', 60,...
    'fontname','arial','fontsize', 12);
ylabel('Accuracy (%)', 'fontname','arial','fontsize', 14);

subplot(1,3,2)
bar([sensitivity_train']*100, 'grouped');
colormap(cool); xlim([0 4]); ylim([0 100]);
set (gca,'xticklabel',{'Top Four','Top Five','All'},'xticklabelrotation', 60,...
    'fontname','arial','fontsize', 12);
ylabel('Sensitivity (%)', 'fontname','arial','fontsize', 14);
xlabel('Selected Biomarker No.', 'fontname','arial','fontsize', 14);

subplot(1,3,3)
bar([specificity_train']*100, 'grouped');
colormap(cool); xlim([0 4]); ylim([0 100]);
set (gca,'xticklabel',{'Top Four','Top Five','All'},'xticklabelrotation', 60,...
    'fontname','arial','fontsize', 12);
ylabel('Specificity (%)', 'fontname','arial','fontsize', 14);

print (fh, '-dtiff','-r300','Result/compareTop_top4_topAll_withAll2.tiff');
print (fh, '-dpdf','-r300','Result/compareTop_top4_topAll_withAll2.pdf');


writetable(table(AUC_train', accuracy_train', sensitivity_train', specificity_train',... 
        'VariableNames', {'AUC_train','accuracy_train', 'sensitivity_train', 'specificity_train'},...
        'RowNames',{'Top Four','Top Five','All'}),...
        'Result/accuracy_top_top5_topall2.csv', 'WriteRowNames',1);

%%  with one file for three methods

%%  with one file for three methods

[score, featureIdsorted] = featureRank_kfold_allThreeMethods(peli_carc, classlabel, 10);
%save('kfold_score_three_rankingMethod_feb18_carc.mat','score','featureIdsorted');
%% plot three score in one figure
h3 = figure;
set(h3, 'PaperUnits','inches','Units','inches','Position',[2 2 9 3.5], ...
            'PaperSize',[9 3.5], 'PaperPosition',[0 0 9 3.5]);
subplot('position',[0.055 0.23 0.27 0.75])
drawFigures_allThreeScore (score.tstatScore, featureIdsorted.tstatScore, geneName_N, 't-stat Score');
annotation('textbox','position',[0.04, 0.01, 0.27, 0.05],'String','a) Scoring Criteria: t-stat',...
       'linestyle', 'none','fontname','Arial','fontsize',12, 'fontWeight','bold',...
       'margin',0,'HorizontalAlignment','center','VerticalAlignment','bottom');


subplot('position',[0.387 0.23 0.27 0.75])
drawFigures_allThreeScore (score.mrmrTCD, featureIdsorted.mrmrTCD, geneName_N, 'MRMR-TCD Score');
annotation('textbox','position',[0.35, 0.01, 0.33, 0.05],'String','b) Scoring Criteria: MRMR-TCD',...
       'linestyle', 'none','fontname','Arial','fontsize',12, 'fontWeight','bold',...
       'margin',0,'HorizontalAlignment','center','VerticalAlignment','bottom');
set(gca,'ylim',[-0.05, 1]);   
subplot('position',[0.722 0.23 0.27 0.75])
drawFigures_allThreeScore (score.mrmrTCQ, featureIdsorted.mrmrTCQ, geneName_N, 'MRMR-TCQ Score');
annotation('textbox','position',[0.69, 0.01, 0.33, 0.05],'String','c) Scoring Criteria: MRMR-TCQ',...
       'linestyle', 'none','fontname','Arial','fontsize',12, 'fontWeight','bold',...
       'margin',0,'HorizontalAlignment','center','VerticalAlignment','bottom');


%print (h3, '-dtiff', '-r300', 'carc_kfold_Nov_17.tiff');
print (h3, '-dpdf', '-r300', 'carc_kfold_Feb_18.pdf');

%% for biomarker 1- 5
[trainData, testData] = getTrainTestData(peli_carc, classlabel);

id = featureIdsorted.mrmrTCQ (1:5);

classModel_top5 = fitcsvm(peli_carc(:,id), classlabel, 'kfold', 10, 'kernelfunction', 'linear');

%classModel = fitcsvm(trainData.vls(:,id), trainData.lbls, 'kfold', 10);    
[accuracy_train(1), sensitivity_train(1), specificity_train(1), AUC_train(1)] = ...
        performanceCriteria_CV(classModel_top5, classlabel);
    
id = featureIdsorted.mrmrTCQ; 
classModel_all = fitcsvm(peli_carc(:,id), classlabel, 'kfold', 10, 'kernelfunction', 'linear');
[accuracy_train(2), sensitivity_train(2), specificity_train(2), AUC_train(2)] = ...
        performanceCriteria_CV(classModel_all, classlabel);

writetable(table(AUC_train', accuracy_train', sensitivity_train', specificity_train',... 
        'VariableNames', {'AUC_train','accuracy_train', 'sensitivity_train', 'specificity_train'},...
        'RowNames',{'Top Five','All'}),...
        'carc_accuracy_top5_topall_nov2017.csv', 'WriteRowNames',1);
    
    save ('svmModel.mat','classModel_top5','classModel_all','classlabel');
