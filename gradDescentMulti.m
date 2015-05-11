function [H,R,RClean,RNoise] = gradDescentMulti(M,MClean,MNoise,param)
% Train a set of filters using gradient descent.  
% This version for multi-condition training. MClean and MNoise are
% cell arrays of data matrices for each noise mix condition.
% Returns:
%  H - learned filters
%  R - responses

% ------------------------- Pre-processing --------------------------------

disp('Performing pre-processing...')

% Whitening and dimension reduction
if param.whiten || ~isempty(param.nDim)
    [M, WM, DWM] = whiten(M, param.nDim);
    % Apply to all cases
    for i = 1:length(MClean)
        MClean{i} = WM'*MClean{i};
        MNoise{i} = WM'*MNoise{i};
    end
end

% If pre-processing only, return.
if param.costOpt==0
    H = WM;
    R = M;
    RClean = MClean;
    RNoise = MNoise;
    return
end

% ------------------------- Main algorithm --------------------------------

% initialize the filters and orthogonalize
H = randn(size(M,1),size(M,1));
H = H/sqrtm(H'*H);
RClean = cell(size(MClean));
RNoise = cell(size(MNoise));
RCleanMean = cell(size(MClean));
RNoiseMean = cell(size(MNoise));
Beta = cell(size(MNoise));
for i = 1:length(MClean)
    RClean{i} = H'*MClean{i};
    RNoise{i} = H'*MNoise{i};
end

% If random filters only, return. 
if param.costOpt == 4, 
    H = WM*H;
    return
end

% Set M large arbitrary initial cost
cost = 1.e14;
oldCost = cost;

% Loop through grad descent steps
for iter = 1:param.maxIter
    
    disp(['Running step ', num2str(iter)])

    % Compute new cost function
    switch param.costOpt
        
        case 5 % Robustness
            cost = 0;
            for i = 1:length(RClean)
                RCleanMean{i} = mean(RClean{i},2);
                RNoiseMean{i} = mean(RNoise{i},2);
                Beta{i} = mean((RClean{i} - repmat(RCleanMean{i},1,size(RClean{i},2))).^2,2).*mean((RNoise{i} - repmat(RNoiseMean{i},1,size(RNoise{i},2))).^2,2);
                cost = cost - ...
                    mean( mean((RClean{i} - repmat(RCleanMean{i},1,size(RClean{i},2))).*(RNoise{i} - repmat(RNoiseMean{i},1,size(RClean{i},2))),2)./ ...
                          sqrt(Beta{i}) )...
                    /length(RClean);
            end
            
        otherwise
            error('invalid costopt parameter');
            
    end
    
    % Compare cost functions and update if there was M decrease
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

        case 5
            grad = zeros(size(H));
            % THE FOLLOWING FORMULA WORKS ONLY BECAUSE THE MEANS OF CLEAN
            % AND NOISY DATA HAVE ALREADY BEEN REMOVED.
            for i = 1:length(RClean)
                grad = grad - ...
                    (MNoise{i}*RClean{i}' + MClean{i}*RNoise{i}') /size(MClean{i},2) *diag(Beta{i}.^(-1/2)) - ...
                    cost*( MClean{i}*RClean{i}'*diag(mean(RNoise{i}.^2,2)) + ...
                           MNoise{i}*RNoise{i}'*diag(mean(RClean{i}.^2,2)) ) /size(MClean{i},2) *diag(Beta{i}.^-1) ...
                    /length(RClean);
            end
            
    end
    
    % Update filters and orthogonlaize
    H = H - param.learnRate*grad;
    H = H/sqrtm(H'*H);
    for i = 1:length(MClean)
        RClean{i} = H'*MClean{i};
        RNoise{i} = H'*MNoise{i};
    end
    
end

% ------------------------- Post-processing -------------------------------

% De-whiten the filters
H = WM*H;
M = DWM'*M;
R = H'*M;
for i = 1:length(MClean)
    MClean{i} = DWM'*MClean{i};
    MNoise{i} = DWM'*MNoise{i};
    RClean{i} = H'*MClean{i};
    RNoise{i} = H'*MNoise{i};
end





