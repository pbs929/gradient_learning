function dispFilts(H, costs, param)

figure(7)
set(0,'defaulttextfontsize',18)
switch param.costOpt
    case 0
        plotFilters(H,5,10,[],param);
    case {-1,1,2,3,4,5}
        [costsOrd, order] = sort(costs,'descend');
        plotFilters(H(:,order),5,10,round(costsOrd),param);
end