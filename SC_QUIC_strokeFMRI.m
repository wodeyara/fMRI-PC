% for lambda = 1:13
allSources = useAreas;
numSources = length(allSources);
orig_G  =double(~(eye(length(allSources)))) .* SC(allSources, allSources) > 0;
GforFit = (eye(length(orig_G)) +double(orig_G>0)); 

lambda = 10; %0.125

allSourceCov = zeros(length(allTimeSeriesFmri), length(allSources), length(allSources));
allSourceCoh =  zeros(length(allTimeSeriesFmri), length(allSources), length(allSources));
allSourcePrec =  zeros(length(allTimeSeriesFmri), length(allSources), length(allSources));
allSourceParCoh =  zeros(length(allTimeSeriesFmri), length(allSources), length(allSources));
allSourceReconCoh =  zeros(length(allTimeSeriesFmri), length(allSources), length(allSources));

for j = 1:length(allTimeSeriesFmri)
    sourceData = (allTimeSeriesFmri{j})';
    
    if ~isempty(sourceData)
        sourceData = sourceData(:, allSources);
        N = size(sourceData,1);	    

        dataCov = cov(sourceData);

        sourceCov = dataCov; % real valued
        allSourceCov(j,:,:) = sourceCov;
        sourceCoh = abs(sourceCov)./sqrt(diag(sourceCov) * diag(sourceCov)'); % getting coherence
        allSourceCoh(j,:,:) = sourceCoh;
% 
        [X,~,~, cputime] = QUIC('default', sourceCov, ... 
            allLambdas(lambda)  * max(max(triu(abs(sourceCov),1))) +  ... 
            (max(allLambdas)-allLambdas(lambda))*max(max(triu(abs(sourceCov),1)))*  ... 
                double(~(GforFit) + eye(length(GforFit))), 1e-10, 0, 300);

%         [X,~,~, cputime] = QUIC('default', sourceCov, ... 
%             allLambdas(lambda)  * max(max(triu(abs(sourceCov),1))), 1e-10, 0, 300);


        newG1 = abs(X)>0 ;
        sum(sum(triu(newG1,1)))
        sum(sum(triu(newG1.*(orig_G),1)))
        
        GforFit1 = newG1;

        P = ggmFitHtf(sourceCov+ eye(length(sourceCov)) *min(allLambdas)  * ...
            max(max(triu(abs(sourceCov),1))),GforFit1);

        %what is deviance?
        lassoMdlDev(j) = deviance(sourceCov + eye(length(sourceCov)) * min(allLambdas) ...
            * max(max(triu(abs(sourceCov),1))), P,1);

        % compare coherence reconstruction from model fit and error term
        P1 = P;
        allSourcePrec(j,:,:) = P1;
        allSourceParCoh(j,:,:) = abs(P1./ sqrt(diag(P1)*diag(P1)'));

        precReconCov_Lasso = inv(P1);
        precReconCoh_Lasso  = abs(precReconCov_Lasso)./sqrt(diag(precReconCov_Lasso) ...
            * diag(precReconCov_Lasso)');
        allSourceReconCoh(j,:,:) = precReconCoh_Lasso;

        clear sourceData dataCov dataUse sourceCov sourceData ...
            tmpCov X P1 precReconCov_Lasso precReconCoh_Lasso
    end
end
% allMdlDev(lambda,:) = lassoMdlDev;
% end