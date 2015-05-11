function cost = costKurt(R)

R_mean = repmat(mean(R,1), size(R,1), 1);
mu2 = mean((R-R_mean).^2, 1);
mu4 = mean((R-R_mean).^4, 1);
cost = mean(mu4./(mu2.^2));
