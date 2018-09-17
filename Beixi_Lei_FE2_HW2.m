%Define all the parameters
alpha = 1;
beta = 0.6;
f = @(x,mu) 1/sqrt(2*pi)*exp(-(((x(2:1001)-(alpha+beta*x(1:1000)) ^ 2)/2)));
mu = [alpha,beta];
x0 = [1,1];

num  = 1000; 
x = mvnrnd(mu,sigma,num);
x(1)=2;
k = rand(num+1,1);

for i = 1:num
    x(i+1) = alpha + beta*x(i) + k(i+1);
end

%Define likelihood function and max theta
objfun = @(mu) - sum(log(f(x,mu)));
options = optimset('LargeScale','off');
[mu_hat, fval, exitflag, output,grad,H] = fminunc(objfun,x0,options);
%estimate variance in MLE 
lambda_hat = (H/num)^(-1);

%p-value, H0: mu_hat = mu
z1=  sqrt(num)*(mu_hat(1)-mu(1))/sqrt(lambda_hat(1,1));
z2 = sqrt(num)*(mu_hat(2)-mu(2))/sqrt(lambda_hat(2,2));
p1 = 2*(1-normcdf(abs(z1)));
p2 = 2*(1-normcdf(abs(z2)));

%confidence interval
alpha = 0.05;
ub1 = (mu_hat(1)) + norminv(1 - alpha/2) * sqrt(lambda_hat(1,1)/num);
lb1 = (mu_hat(1)) - norminv(1 - alpha/2) * sqrt(lambda_hat(2,2)/num);
ub2 = (mu_hat(2)) + norminv(1 - alpha/2) * sqrt(lambda_hat(1,1)/num);
lb2 = (mu_hat(2)) - norminv(1 - alpha/2) * sqrt(lambda_hat(2,2)/num);

%print all the results
fprintf("The confidence interval for x1 is : " , lb1, ub1);
fprintf("The confidence interval for x2 is : " , lb2, ub2);
fprintf("The p-value for X1 is : ", p1);
fprintf("The p-value for X2 is : ", p2);
fprintf("The true mu is : ", mu(1), mu(2));
fprintf("The estimate mu is : " , mu_hat(1), mu_hat(2));