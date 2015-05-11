function [costSus, costKrt, cstICA, costRob] = dispCosts(RTst, RTstClean, RTstNoise, param)

disp('---')

% Display cost functions
costSus = costSF(RTst,param);
disp(['Sustained firing objective: ',num2str(mean(costSus))])
costKrt = costKurt(RTst);
disp(['Kurtosis objective: ',num2str(mean(costKrt))])
[cost1,cost2,cost3] = costICA(RTst);
cstICA = cost1;
disp(['ICA objective: ',num2str(mean(cost1))])
costRob = costRobust(RTstClean,RTstNoise);
disp(['Robustness objective: ',num2str(mean(costRob))])

figure(6); 
subplot(2,2,1); hold all;
plot(sort(costSus,'descend'))
axis tight
title('Sustained firing objective')

subplot(2,2,3); hold all;
plot(sort(cost1,'descend'))
axis tight
title('ICA objective')

subplot(2,2,4); hold all;
plot(costRob)
axis tight
title('Robustness objective')

disp('---')

figure(9); clf
plot(costSus,costRob,'x')