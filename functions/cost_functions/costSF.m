function cost = costSF(R,param)

cost = 0;
for tau = 0:param.deltaT-1
    alpha = 1 - tau/param.deltaT;
    cost = cost + alpha*sum( R.^2.*[R(:,tau+1:end).^2,zeros(size(R,1),tau)], 2 )/size(R,2);
end