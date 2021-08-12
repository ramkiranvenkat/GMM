function v = gaussianEval(x,mu,sig)
% x,mu are row vectors
N = length(x);
c = 1/((2*pi)^(N/2))/(det(sig)^0.5);
v = c*exp(-0.5*(x-mu)*pinv(sig)*(x-mu)');
