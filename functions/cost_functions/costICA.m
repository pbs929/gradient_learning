function [cost1,cost2,cost3] = costICA(R)
% The objective function for ICA is an approximation of the negentropy,
% averaged over the output channels.  
%
% Negentropy is defined for a single random variable as 
%   J(x) = H(x_gauss) - H(x)
% where x_gauss is a gaussian variable with identical variance (and mean).
% Approximations take the form:
%   J(x) ~ c[E{G(x)} - E{G(x_gauss)}]^2
% The outputs of this function are the three cost functions in the
% Hyvarinen (1999) 'Fast and robust...'.  
%
% Not that I have subracted the gaussian form only in the first case

% G(x) = 1/a1*log(cosh(a1*x)), with a1 = 1, J(x) = [G(x) - 0.3746]^2
% This is a good general-purpose contrast function.
cost1 = (mean(log(cosh(R)),2) - 0.3746).^2;

% G(x) = -1/a2*exp(-a2*x^2/2), with z2 = 1, J(x) = [G(x) + sqrt(2)]^2
% This contrast function is good for very supergaussian distributions and 
% when robustness against outliers is important.  
% Since the expectation for the Gaussian case is just an integration of a
% Gaussian, it's easy to show that it equals -sqrt(2).  
cost2 = (-1*mean(exp(R.^2)/2,2) + sqrt(2)).^2;

% G(x) = x^4/4, J(x) = [x^4/4 - 3/4]^2 aka. (kurtosis/4)^2
% 3 is the 4th moment of the normal distribution w/ variance 1.  So this
% contrast function is equivalent to the kurtosis, is >1 for supergaussian
% distributions.
cost3 = (mean(R.^4/4,2) - 3/4).^2; 

return

% Code used to find the gaussian values (numerical integration) in case 2:
% 
% x = -5:.01:5;
% f = @(x)[log(cosh(x)).*exp(-x.^2/2)/sqrt(2*pi)*];
% y = f(x);
% trapz(x,y)



