function [H,R] = gradDescent(A, param)
% Train a set of filters using gradient descent.  
% Returns:
%  H - learned filters
%  A - raw data matrix
%  R - responses

% ------------------------- Pre-processing --------------------------------

disp('Performing pre-processing...')

% Whitening and dimension reduction
if param.whiten || ~isempty(param.nDim)
    [A, WM, DWM] = whiten(A, param.nDim, 1); % add argument 1 for plot of PCA variances
end

% If pre-processing only, return.
if param.costOpt==1 % PCA only
    H = WM;
    R = A;
    return
end

% ------------------------- Main algorithm --------------------------------

% initialize the filters and orthogonalize
H = randn(size(A,1),size(A,1));
H = H/sqrtm(H'*H);
R = H'*A;

% If random filters only, return. 
if param.costOpt == 0, % Random filters
    H = WM*H;
    return
end

% Set a large arbitrary initial cost
cost = 1.e14;
oldCost = cost;

% Loop through grad descent steps
for iter = 1:param.maxIter
    
    disp(['Running step ', num2str(iter)])

    % Compute new cost function
    switch param.costOpt
        
        case 3 % Kurtosis
            R_mean = repmat(mean(R,1), size(R,1), 1);
            mu2 = repmat( mean((R-R_mean).^2, 1), size(R,1), 1);
            mu3 = repmat( mean((R-R_mean).^3, 1), size(R,1), 1);
            mu4 = repmat( mean((R-R_mean).^4, 1), size(R,1), 1);
            cost = -1*mean(mu4(1,:)./(mu2(1,:).^2));
        
        case 4 % Sustained Firing
            cost = 0;
            for tau = 0:param.deltaT-1
                alpha = 1 - tau/param.deltaT;
                cost = cost - alpha*sum(sum( R.^2.*[R(:,tau+1:end).^2,zeros(size(R,1),tau)] ))/numel(R);
            end
        
    end
    
    % Compare cost functions and update if there was a decrease
    disp(['New value of cost function: ', num2str(cost)])
    disp(['Change in cost function: ', num2str(cost - oldCost)])
    if cost<oldCost
        oldCost = cost;
    else
        disp(['*** Decreasing learning rate to ',num2str(param.learnRate/2),' ***'])
        param.learnRate = param.learnRate/2;
    end
    
    % Compute gradient
    switch param.costOpt
        
        case 3 % Kurtosis
            grad = -1*A*( ((R-R_mean).^3 - mu3)./(mu2.^2) - mu4.*(R-R_mean)./(mu2.^3) )'*4/numel(A);
        
        case 4 % Sustained Firing
            grad = zeros(size(H));
            for tau = 0:param.deltaT-1
                alpha = 1 - tau/param.deltaT;
                grad = grad - 2*alpha*( A*(R.*[R(:,tau+1:end).^2,zeros(size(R,1),tau)])' + [A(:,tau+1:end),zeros(size(A,1),tau)]*([R(:,tau+1:end),zeros(size(R,1),tau)].*R.^2)' )/numel(R);
            end
            
    end
    
    % Update filters and orthogonlaize
    H = H - param.learnRate*grad;
    H = H/sqrtm(H'*H);
    R = H'*A;
    
end

% ------------------------- Post-processing -------------------------------

% De-whiten the filters
H = WM*H;
A = DWM'*A;
R = H'*A;

