    
function [mdlDev] = est2foldDev_AGL(sourceData,allLambdas, lambda, GforFit)

    timeLength = size(sourceData,1); % assumes data is ROIs by samples
    
    sourceCovFirstHalf = cov(sourceData(1:floor(timeLength/2),:));
    sourceCovSecondHalf = cov(sourceData(1+floor(timeLength/2) :timeLength, :));

    %% FIRST HALF OF DATA

    [X,~,~] = QUIC('default', sourceCovFirstHalf, ... 
        allLambdas(lambda)  * max(max(triu(abs(sourceCovFirstHalf),1))) +  ... 
        (max(allLambdas)-allLambdas(lambda))*max(max(triu(abs(sourceCovFirstHalf),1)))*  ... 
            double(~(GforFit) + eye(length(GforFit))), 1e-10, 0, 300);
        

    GforFit1 = abs(X)>0;

    P = ggmFitHtf(sourceCovFirstHalf+ eye(length(sourceCovFirstHalf)) *min(allLambdas)  * ...
        max(max(triu(abs(sourceCovFirstHalf),1))),GforFit1);

    %what is deviance?
    mdlDev(1) = deviance(sourceCovSecondHalf + eye(length(sourceCovSecondHalf)) * min(allLambdas) ...
        * max(max(triu(abs(sourceCovSecondHalf),1))), P,1);
    
    %% SECOND HALF OF DATA 

    [X,~,~] = QUIC('default', sourceCovSecondHalf, ... 
        allLambdas(lambda)  * max(max(triu(abs(sourceCovSecondHalf),1))) +  ... 
        (max(allLambdas)-allLambdas(lambda))*max(max(triu(abs(sourceCovSecondHalf),1)))*  ... 
            double(~(GforFit) + eye(length(GforFit))), 1e-10, 0, 300);
        

    GforFit1 = abs(X)>0;

    P = ggmFitHtf(sourceCovSecondHalf+ eye(length(sourceCovSecondHalf)) *min(allLambdas)  * ...
        max(max(triu(abs(sourceCovSecondHalf),1))),GforFit1);

    %what is deviance?
    mdlDev(2) = deviance(sourceCovFirstHalf + eye(length(sourceCovFirstHalf)) * min(allLambdas) ...
        * max(max(triu(abs(sourceCovFirstHalf),1))), P,1);
