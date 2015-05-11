function plotFilters(H,nRow,nCol,labels,param)

clim(1) = max(H(:)); clim(2) = min(H(:));

for i=1:min(size(H,2),nRow*nCol)

    STRF = reshape(H(:,i)', size(H,1)/param.delay, param.delay);
    
    subplot(nRow,nCol,i)
    imagesc(STRF); axis xy; 
    caxis([-max(abs(clim)) max(abs(clim))]);
    axis off
    if ~isempty(labels), title(labels(i)); end
    
end
