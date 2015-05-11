function [RTst, RTstClean, RTstNoise] = strfResponses(H, MTst, MTstClean, MTstNoise)

RTst = H'*MTst;
RTstClean = cell(size(MTstClean));
RTstNoise = cell(size(MTstNoise));
for ind = 1:length(MTstClean)
    RTstClean{ind} = H'*MTstClean{ind};
    RTstNoise{ind} = H'*MTstNoise{ind};
end