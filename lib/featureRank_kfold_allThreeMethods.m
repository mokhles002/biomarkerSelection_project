function [score, featureIdsorted] = featureRank_kfold_allThreeMethods(data, classlabel,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform the gene selection from the genotoxicity
%   Developed by: Sheikh M. Rahman
%   Date: 11-14-2017
%   A fucntion that find the score using ttest, MRMR-TCD and MRMR-TCQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kfoldCVP = cvpartition(classlabel,'kfold',k);
nFeature = size(data,2);
tvls = nan(nFeature,k);
mrmrTCD = nan(nFeature,k);
mrmrTCQ = nan(nFeature,k);

for i = 1:k
    tvls(:,i) = getScore(i, 'tstat');
    mrmrTCD(:,i) = getScore(i, 'MRMR-TCD');
    mrmrTCQ(:,i) = getScore(i, 'MRMR-TCQ');
end
[score.tstatScore,featureIdsorted.tstatScore] = sort(mean(tvls,2),'descend'); % sort the features
[score.mrmrTCD,featureIdsorted.mrmrTCD] = sort(mean(mrmrTCD,2),'descend'); % sort the features
[score.mrmrTCQ,featureIdsorted.mrmrTCQ] = sort(mean(mrmrTCQ,2),'descend'); % sort the features

    function score = getScore(i, rnkgMthod)
        dataTrain = data(kfoldCVP.training(i),:);
        lblTrain = classlabel(kfoldCVP.training(i));
        id_pos = strcmp(lblTrain,'positive');
        id_neg = strcmp(lblTrain,'negative');
        dataTrainG1 = dataTrain(id_pos,:);
        dataTrainG2 = dataTrain(id_neg,:);
        
        %Relevance
        [~,~,~,stat] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal');
        tvalue = (stat.tstat).^2; % score is the t-square
        rlvnce = tvalue/nFeature;
        score = rlvnce;
        
        % Redundancy
        if ~strcmp(rnkgMthod, 'tstat')
            rho = corr(dataTrain);
            redundancy = (sum(abs(rho),1))/(nFeature^2);
            switch rnkgMthod
                case 'MRMR-TCD'
                    mRmR = rlvnce - redundancy;
                case 'MRMR-TCQ'
                    mRmR = rlvnce ./ redundancy;
            end
            score = mRmR;
        end        
    end

end