function cost = costRobust(RClean,RNoise)
% Averages costs for several noise conditions.  The cost function is the
% Pearson correlation coefficient between the clean response and
% noise-corrupted response.  

cost = 0;

for i = 1:length(RClean)
    cost = cost + noiseCost(RClean{i}, RNoise{i})/length(RClean);
end

function cost = noiseCost(RClean, RNoise)
% This for a single noise condition

cost = mean((RClean - repmat(mean(RClean,2),1,size(RClean,2))).*(RNoise - repmat(mean(RNoise,2),1,size(RClean,2))),2)./ ...
       sqrt(mean((RClean - repmat(mean(RClean,2),1,size(RClean,2))).^2,2).*mean((RNoise - repmat(mean(RNoise,2),1,size(RNoise,2))).^2,2));