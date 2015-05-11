function [H,R] = gradDescentDeflate(A, param)
% Train a set of filters using gradient descent.  
% Returns:
%  H - learned filters
%  A - raw data matrix
%  R - responses

% ------------------------- Pre-processing --------------------------------

disp('Performing pre-processing...')

% Whitening and dimension reduction
if param.whiten || ~isempty(param.nDim)
    [A, WM, DWM] = whiten(A, param.nDim);
end

% If pre-processing only, return.
if param.costOpt==1
    H = WM;
    R = A;
    return
end

% ------------------------- Main algorithm --------------------------------

H = zeros(size(A,1),size(A,1));

% Iterate deflation steps
for iDefl = 1:param.nDim
    
    disp(['Calculating filter #',num2str(iDefl)])
    
    % Find a basis vector
    h = deflStep(A, H, param);
    H(:,iDefl) = h;
    
end

% ------------------------- Post-processing -------------------------------

% De-whiten the filters
H = WM*H;
A = DWM'*A;
R = H'*A;


function h = deflStep(A, H, param)

% initialize the filter
h = randn(size(A,1),1);
h = h - H*H'*h;   % orthogonalize with the other vectors
h = h/norm(h);    % normalize
r = h'*A;

% If random filters only, return. 
if param.costOpt == 0, return; end

% Set a large arbitrary initial cost
cost = 1.e14;
oldCost = cost;

% Loop through grad descent steps
for iter = 1:param.maxIter
    
    disp(['running step ', num2str(iter)])

    % Compute new cost function
    switch param.costOpt
        
        case 4 % Sustained Firing
            cost = 0;
            for tau = 0:param.deltaT-1
                alpha = 1 - tau/param.deltaT;
                cost = cost - alpha*sum(sum( r.^2.*[r(:,tau+1:end).^2,zeros(size(r,1),tau)] ))/numel(r);
            end
            
        otherwise
            error('invalid costopt parameter');
        
    end
    
    % Compare cost functions and update if there was a decrease
    disp(['New value of cost function: ', num2str(cost)])
%     disp(['Change in cost function: ', num2str(cost - oldCost)])
    if cost<oldCost
        oldCost = cost;
    else
        disp(['*** Decreasing learning rate to ',num2str(param.learnRate/2),' ***'])
        param.learnRate = param.learnRate/2;
    end
    
    % Compute gradient
    switch param.costOpt
                
        case 4 % Sustained Firing
            grad = zeros(size(h));
            for tau = 0:param.deltaT-1
                alpha = 1 - tau/param.deltaT;
                grad = grad - 2*alpha*( A*(r.*[r(:,tau+1:end).^2,zeros(size(r,1),tau)])' + [A(:,tau+1:end),zeros(size(A,1),tau)]*([r(:,tau+1:end),zeros(size(r,1),tau)].*r.^2)' )/numel(r);
            end
            
    end
    
    % Update filters and orthogonlaize
    h = h - param.learnRate*grad;
    h = h - H*H'*h;
    h = h/norm(h);
    r = h'*A;
    
end


