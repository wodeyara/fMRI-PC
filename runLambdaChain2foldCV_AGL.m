% I'm rewriting the SC_QUIC_strokeFMRI code to do 2-fold cross validation
% so that we can decide on the correct penalization to use. 

allSources = useAreas; % losing all subcortical areas
numSources = length(allSources);
orig_G  =double(~(eye(length(allSources)))) .* SC(allSources, allSources) > 0;
GforFit = (eye(length(orig_G)) +double(orig_G>0)); %adding a diagonal as this is necessary for QUIC line

allMdlDev  = 0;
for lambda = 1:length(lambda)
    tic
    for j = 1:length(allTimeSeriesFmri)
        sourceData = (allTimeSeriesFmri{j})';

        if ~isempty(sourceData)
            sourceData = sourceData(:, allSources);
            N = size(sourceData,1);	    
           
            mdlDev = est2foldDev_AGL(sourceData, allLambdas, lambda, GforFit);
            allMdlDev(lambda,j) = mean(mdlDev);
            
        end
    end
    toc
end
clear sourceData N orig_G GforFit allSources numSources
