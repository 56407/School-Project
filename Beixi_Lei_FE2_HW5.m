%Initialize parameters
%[r,~,~] = xlsread('FRB_H15.csv');
r = FRBH15.VarName2;
num = length(r);

%Maximumization routine
theta0 = [0.001,0.02,0.4,0.25,0.045,0.06,0.5,0.65,0.95];
objfunc = @(theta) -likelihood(r,theta(1:2),theta(3:4),theta(5:6),theta(7),theta(8),theta(9)); 
[thetahat,fval,~,~,~,H] = fminunc(objfunc,theta0); 

%covariance matrix
lambdahat = (H/num)^(-1);
%H0: thetahat = theta
z1 = sqrt(num)*(thetahat(1) - mu(1)) ./ sqrt(lambdahat(1,1));
z2 = sqrt(num)*(thetahat(2) - mu(2)) ./ sqrt(lambdahat(2,2));
z3 = sqrt(num)*(thetahat(3) - mu(3)) ./ sqrt(lambdahat(3,3));
z4 = sqrt(num)*(thetahat(4) - mu(4)) ./ sqrt(lambdahat(4,4));
z5 = sqrt(num)*(thetahat(5) - mu(5)) ./ sqrt(lambdahat(5,5));
z6 = sqrt(num)*(thetahat(6) - mu(6)) ./ sqrt(lambdahat(6,6));
z7 = sqrt(num)*(thetahat(7) - mu(7)) ./ sqrt(lambdahat(7,7));
z8 = sqrt(num)*(thetahat(8) - mu(8)) ./ sqrt(lambdahat(8,8));
p1 = 2*(1-normcdf(abs(z1)));
p2 = 2*(1-normcdf(abs(z2)));
p3 = 2*(1-normcdf(abs(z3)));
p4 = 2*(1-normcdf(abs(z4)));
p5 = 2*(1-normcdf(abs(z5)));
p6 = 2*(1-normcdf(abs(z6)));
p7 = 2*(1-normcdf(abs(z7)));
p8 = 2*(1-normcdf(abs(z8)));

%Confidence interval
alpha = 0.05;
ub1 = (thetahat(1)) + norminv(1-alpha/2)*sqrt(lambdahat(1,1)/num);
lb1 = (thetahat(1)) - norminv(1-alpha/2)*sqrt(lambdahat(1,1)/num);
ub2 = (thetahat(2)) + norminv(1-alpha/2)*sqrt(lambdahat(2,2)/num);
lb2 = (thetahat(2)) + norminv(1-alpha/2)*sqrt(lambdahat(2,2)/num);
ub3 = (thetahat(3)) + norminv(1-alpha/2)*sqrt(lambdahat(3,3)/num);
lb3 = (thetahat(3)) - norminv(1-alpha/2)*sqrt(lambdahat(3,3)/num);
ub4 = (thetahat(4)) + norminv(1-alpha/2)*sqrt(lambdahat(4,4)/num);
lb4 = (thetahat(4)) - norminv(1-alpha/2)*sqrt(lambdahat(4,4)/num);
ub5 = (thetahat(5)) + norminv(1-alpha/2)*sqrt(lambdahat(5,5)/num);
lb5 = (thetahat(5)) - norminv(1-alpha/2)*sqrt(lambdahat(5,5)/num);
ub6 = (thetahat(6)) + norminv(1-alpha/2)*sqrt(lambdahat(6,6)/num);
lb6 = (thetahat(6)) - norminv(1-alpha/2)*sqrt(lambdahat(6,6)/num);
ub7 = (thetahat(7)) + norminv(1-alpha/2)*sqrt(lambdahat(7,7)/num);
lb7 = (thetahat(7)) - norminv(1-alpha/2)*sqrt(lambdahat(7,7)/num);
ub8 = (thetahat(8)) + norminv(1-alpha/2)*sqrt(lambdahat(8,8)/num);
lb8 = (thetahat(8)) - norminv(1-alpha/2)*sqrt(lambdahat(8,8)/num);

%Print results
fprintf("True theta = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]\n",theta(1),theta(2),theta(3),theta(4),theta(5),theta(6),theta(7),theta(8));
fprintf("Estimated theta_hat = [%.2f, %.2f, %.2f, %.2f,%.2f, %.2f, %.2f, %.2f]\n",thetahat(1),thetahat(2),thetahat(3),thetahat(4),thetahat(5),thetahat(6),thetahat(7),thetahat(8));
fprintf("Maximized loglikelihood = %.2f\n", -fval);
fprintf("p-value1 = %.8f\n", p1);
fprintf("p-value2 = %.8f\n", p2);
fprintf("p-value3 = %.8f\n", p3);
fprintf("p-value3 = %.8f\n", p4);
fprintf("p-value3 = %.8f\n", p5);
fprintf("p-value3 = %.8f\n", p6);
fprintf("p-value3 = %.8f\n", p7);
fprintf("p-value3 = %.8f\n", p8);
fprintf("Confidence interval for alpha0 at level = %.2f: [%f, %f]\n", alpha, lb1, ub1);
fprintf("Confidence interval for alpha1 at level = %.2f: [%f, %f]\n", alpha, lb2, ub2);
fprintf("Confidence interval for beta0 at level = %.2f: [%f, %f]\n", alpha, lb3, ub3);
fprintf("Confidence interval for beta1 at level = %.2f: [%f, %f]\n", alpha, lb4, ub4);
fprintf("Confidence interval for sigma0 at level = %.2f: [%f, %f]\n", alpha, lb5, ub5);
fprintf("Confidence interval for sigma1 at level = %.2f: [%f, %f]\n", alpha, lb6, ub6);
fprintf("Confidence interval for gamma0 at level = %.2f: [%f, %f]\n", alpha, lb7, ub7);
fprintf("Confidence interval for gamma1 at level = %.2f: [%f, %f]\n", alpha, lb8, ub8);


%Size Analysis
num = 14018;
numpath = 100;
SizeAnalysis = zeros(numpath,8);
P = [0.95,0.05;0.15,0.85];
statenames = ['State 1','State 2'];
mc = dtmc(P,'StateNames',statenames);
x0 = ones(1,mc.numstates)*50;
x = simulate(mc,num,'x0',x0);
alpha = [0.001,0.01];
beta = [0.30,0.20];
sigma = [0.045,0.06];
gamma = [0.50,0.99];
r = zeros(num,numpath);

for i = 1: numpath
    theta0 = [0.001, 0.01, 0.30, 0.20, 0.045, 0.06, 0.5,0.9,0.75];
    objfunc = @(theta) -likelihood(r,theta(1:2),theta(3:4),theta(5:6),theta(7),theta(8),theta(9)); 
    [thetahat,fval,~,~,~,H] = fmincon(objfunc,theta0);
    theta = [0.002,0.021,0.26,0.18,0.00038,0.049,0.47,0.89,0.75];
    lambdahat = (H/num)^(-1);
    z1 = sqrt(num)*(theta_hat(1) - theta(1)) ./ sqrt(Lambda_hat(1,1));
    z2 = sqrt(num)*(theta_hat(2) - theta(2)) ./ sqrt(Lambda_hat(2,2));
    z3 = sqrt(num)*(theta_hat(3) - theta(3)) ./ sqrt(Lambda_hat(3,3));
    z4 = sqrt(num)*(theta_hat(4) - theta(4)) ./ sqrt(Lambda_hat(4,4));
    z5 = sqrt(num)*(theta_hat(5) - theta(5)) ./ sqrt(Lambda_hat(5,5));
    z6 = sqrt(num)*(theta_hat(6) - theta(6)) ./ sqrt(Lambda_hat(6,6));
    z7 = sqrt(num)*(theta_hat(7) - theta(7)) ./ sqrt(Lambda_hat(7,7));
    z8 = sqrt(num)*(theta_hat(8) - theta(8)) ./ sqrt(Lambda_hat(8,8));
    p1 = 2 * (1 - normcdf(abs(z1)));
    p2 = 2 * (1 - normcdf(abs(z2)));
    p3 = 2 * (1 - normcdf(abs(z3)));
    p4 = 2 * (1 - normcdf(abs(z4)));
    p5 = 2 * (1 - normcdf(abs(z5)));
    p6 = 2 * (1 - normcdf(abs(z6)));
    p7 = 2 * (1 - normcdf(abs(z7)));
    p8 = 2 * (1 - normcdf(abs(z8)));
    stop = 1;
    
    SizeAnalysis(i,:) = [p1<0.05,p2<0.05,p3<0.05,p4<0.05,p5<0.05,p6<0.05,p7<0.05,p8<0.05];
end

SizeAnalysis = sum(SizeAnalysis)/numpath;
fprintf("Percentage of rejection for [alpha0,alpha1,beta0,beta1,sigma0,sigma1,gamma0,gamma1] = [[%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f]\n",SizeAnalysis(1),SizeAnalysis(2),SizeAnalysis(3),SizeAnalysis(4),SizeAnalysis(5),SizeAnalysis(6),SizeAnalysis(7),SizeAnalysis(8));

%Power Analysis 

%Get sample
num = 14018;
numpath = 100;
SizeAnalysis2 = zeros(numpath,8);
P = [0.64,0.36;0.25,0.75];
statenames = ['State 1','State 2'];
mc = dtmc(P,'StateNames',statenames);
x0 = ones(1,mc.numstates)*50;
x = simulate(mc,num,'x0',x0);
alpha = [0.002,0.02];
beta = [0.28,0.17];
sigma = [0.0045,0.03];
gamma = [0.482,0.962];
r = zeros(num,numpath);

for i = 1: numpath
    theta0 = [0.002, 0.02, 0.28, 0.17, 0.046, 0.03, 0.482, 0.68,0.90];
    objfunc = @(theta) -likelihood(r,theta(1:2),theta(3:4),theta(5:6),theta(7),theta(8),theta(9)); 
    [thetahat,fval,~,~,~,H] = fmincon(objfunc,theta0);
    theta = [0.001,0.01,0.3,0.2,-0.00045,0.00006,0.5,0.87,0.69];
    lambdahat = (H/num)^(-1);
    z1 = sqrt(num)*(theta_hat(1) - theta(1)) ./ sqrt(Lambda_hat(1,1));
    z2 = sqrt(num)*(theta_hat(2) - theta(2)) ./ sqrt(Lambda_hat(2,2));
    z3 = sqrt(num)*(theta_hat(3) - theta(3)) ./ sqrt(Lambda_hat(3,3));
    z4 = sqrt(num)*(theta_hat(4) - theta(4)) ./ sqrt(Lambda_hat(4,4));
    z5 = sqrt(num)*(theta_hat(5) - theta(5)) ./ sqrt(Lambda_hat(5,5));
    z6 = sqrt(num)*(theta_hat(6) - theta(6)) ./ sqrt(Lambda_hat(6,6));
    z7 = sqrt(num)*(theta_hat(7) - theta(7)) ./ sqrt(Lambda_hat(7,7));
    z8 = sqrt(num)*(theta_hat(8) - theta(8)) ./ sqrt(Lambda_hat(8,8));
    p1 = 2 * (1 - normcdf(abs(z1)));
    p2 = 2 * (1 - normcdf(abs(z2)));
    p3 = 2 * (1 - normcdf(abs(z3)));
    p4 = 2 * (1 - normcdf(abs(z4)));
    p5 = 2 * (1 - normcdf(abs(z5)));
    p6 = 2 * (1 - normcdf(abs(z6)));
    p7 = 2 * (1 - normcdf(abs(z7)));
    p8 = 2 * (1 - normcdf(abs(z8)));
    stop = 1;
    
    SizeAnalysis2(i,:) = [p1<0.05,p2<0.05,p3<0.05,p4<0.05,p5<0.05,p6<0.05,p7<0.05,p8<0.05];
end

SizeAnalysis2 = sum(SizeAnalysis2)/numpath;
fprintf("Percentage of rejection for [alpha0,alpha1,beta0,beta1,sigma0,sigma1,gamma0,gamma1] = [[%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f]\n",SizeAnalysis(1),SizeAnalysis(2),SizeAnalysis(3),SizeAnalysis(4),SizeAnalysis(5),SizeAnalysis(6),SizeAnalysis(7),SizeAnalysis(8));
%SizeAnalysis2;

%%Likelihood function for 2-state Markov regime switching model
function LL = likelihood(r,alpha,beta,sigma,gamma,p,q)
P00 = p;  
P10 = 1-p; 
P11 = q; 
P01 = 1-q; 
num = length(r);  
%probability matrix of S(t) given phi(t-1)
ps_matrix = zeros(num - 1, 2);    
%conditonal density function for r(t) given S(t) and phi(t-1)
pdf1 = zeros(num - 1, 2);  
for i = 1:(num - 1)
    pdf1(i,1) = 1/(sqrt(2*pi)*sigma(1)^2*r(i)^2*gamma(1))*exp(-(r(i+1)-alpha((1)+beta(1)*r(i)))^2/(2*sigma(1)^2*r(i)^2*gamma(2)));
    pdf1(i,2) = 1/(sqrt(2*pi)*sigma(2)^2*r(i)^2*gamma(2))*exp(-(r(i+1)-alpha((2)+beta(2)*r(i)))^2/(2*sigma(2)^2*r(i)^2*gamma(2)));
end
%conditional density function for r(t) given phi(t-1)
pdf2 = zeros(num - 1,1);   
PS1 = (1-p)/(2-p-q);
PS2 = (1-q)/(2-p-q);
for i = 1:(num - 1)
    ps_matrix(i,1) = p*PS1 + (1-q)*PS2;
    ps_matrix(i,2) = (1-p)*PS1 + q*PS2;
    pdf2(i) = pdf1(i,1) * ps_matrix(i,1) + pdf1(i,2)*ps_matrix(i,2);
    PS1 = (pdf1(i,1)*ps_matrix(i,1))/pdf2(i);
    PS2 = (pdf1(i,2)*ps_matrix(i,2))/pdf2(i);
end
LL = sum(log(pdf2));
end
   
%SizeAnalysis = 
%    0.0782    0.0190    0.2740    0.2870    0.3100    0.0800    0.0300  0.0012   0.0570
%SizeAnalysis2 = 
%    0.0001    0.0280    0.0430    0.0250    0.0620    0.1200    0.1100  0.0970    0.0820
