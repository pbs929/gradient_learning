function [A, WM, DWM] = whiten(A, nDim, plotopt)

% perform PCA
[V,D] = eig(A*A'/size(A,2));

if nargin > 2 && plotopt
    figure
    plot(flipud(diag(D)))
    xlabel('PCA #')
    ylabel('variance')
    d = diag(D);
    disp(sum(d(1:50))/sum(d))
end
    
% dimension reduction
if ~isempty(nDim)
    disp(['Keeping ',num2str(nDim),' of ',num2str(length(D)),' dimensions, ',num2str(sum(diag(D(end-nDim+1:end,end-nDim+1:end)))/sum(diag(D))*100),'% variance retained']);
    V = V(:,end:-1:end-nDim+1);
    D = D(end:-1:end-nDim+1,end:-1:end-nDim+1);
end

% whitening - If no whitening selected, take the whitening and
% dewhitening matrices to be the projections onto the principal
% components.   
% if param.whiten
%     WM = V*diag(diag(D).^-0.5);
%     DWM = (D.^0.5)*V';
% else
    WM = V;
    DWM = V';
% end

% apply to data matrix
A = WM'*A;

