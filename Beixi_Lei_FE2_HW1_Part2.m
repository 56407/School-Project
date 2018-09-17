mu =0;
sigma = 0.1;
n = 100;
x = normrnd(mu,sigma,n);

%Likelihood function
fun = @(x,mu)(1/(sqrt(2*pi)*0.1)*exp(-(x-mu)^2/(2*0.1^2)));
L = @(mu)(-sum(log(fun(x,mu))));
LL = @(theta) (-L);

%Find maximization values
%gives initial value
mu1 = 0;
[muhat,Lval] = fminsearch(LL,mu1);

%Find confidence interval and p-value
[muhat, sigmahat, muc, sigmac] = normfit(x,0.05);

%Calculate the p-value and confidence interval 
average = mean(x);
CI = norminv([0.025,0.975],muhat,0,1);
pvalue = 2*(1 - normcdf(muhat,0,sigma));

 