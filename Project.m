%% CIR
%% 
syms a b c 
syms h x xs
muX=b*(a-x);
sigmaX=c*sqrt(x);
CIR_Density = Density(muX, sigmaX, 4, 5);
%% Plot density function
figure(1)
cir=subs(CIR_Density,{a,b,c,h,xs},{1,1,2,1/250,1});
subplot(3,1,1)
fplot(cir,[-1,5]);
%parameters for real density function
cc=2*a/((1-exp(-a*h))*c^2);
q=2*a*b/c^2-1;
density1=cc*exp(-cc*(x+exp(-a*h)*xs))*(x/(exp(-a*h)*xs))^(q/2)*besseli(q,2*cc*sqrt(x*exp(-a*h)*xs));
v3=subs(density1,{a,b,c,h,xs},{1,1,2,1/250,1});
subplot(3,1,2)
fplot(v3,[-1,5]);
% Error term
subplot(3,1,3)
fplot(v3-cir,[-1,5]);

%% Simulate Data
alpha = 0.01;
beta = 0.01;
sigma = 0.03;
dt = 1/252;
size = 1000;
initial_X = 0.01;
X = simulate(alpha, beta, sigma, dt, size, initial_X);
%%  
Xx = X(1:size-1);
Xs = X(2:size);
names = {'alpha', 'beta', 'sigma'};
theta = [alpha, beta, sigma];
func = matlabFunction(CIR_Density);
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt, Xx, Xs)));
theta_initial = [0.2, 0.1, 0.3];
options = optimset('LargeScale','off');
[theta_hat, fval, exitflag, output, grad,hessian] = fminunc(mleProb,theta_initial,options);
n = length(theta_hat);
% %estimate the sigma
 theta_Sig = sqrt(diag(inv(hessian/n)));
% %p_value and confidence interval
 p_value = zeros(1,n);
 confidence_int = zeros(n,2);
 for i = 1:n
    z_score = abs((theta_hat(i)-theta(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n', names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));
 end

 %% Size Analysis
 Test_Times = 1000;
 reject = zeros(3,1);
 for path = 1:Test_Times
    alpha = 0.01;
    beta =  0.01;
    sigma = 0.03;
    theta = [alpha, beta, sigma];
    dt = 1/252;
    size = 1000;
    initial_X = 0.01;
    X = simulate(alpha, beta, sigma, dt, size, initial_X);
    Xx = X(1:size-1);
    Xs = X(2:size);
    names = {'alpha', 'beta', 'sigma'};
    func = matlabFunction(CIR_Density);
    mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt, Xx, Xs)));
    theta_initial = [0.02, 0, 0.05];
    options = optimset('LargeScale','off');
    [theta_hat, fval, exitflag, output, grad,hessian] = fminunc(mleProb,theta_initial,options);
    n = length(theta_hat);
    % %estimate the sigma
     theta_Sig = sqrt(diag(inv(hessian/n)));
    % %p_value and confidence interval
     p_value_power = zeros(1,n);
     for i = 1:n
        z_score = abs((theta_hat(i)-theta(i))/(theta_Sig(i)*sqrt(n)));
        p_value(i) = 2*(1-normcdf(z_score));
        if p_value(i) < 0.05
            reject(i) = reject(i) + 1;
        end
     end        
 end
 %% Power Analysis
 Test_Times = 1000;
 reject_power = zeros(3,1);
 for path = 1:Test_Times
    alpha = 0.01;
    beta =  0.01;
    sigma = 0.03;
    theta_power = [0.1, 0.2, 0.2];
    dt = 1/252;
    size = 1000;
    initial_X = 0.01;
    X = simulate(alpha, beta, sigma, dt, size, initial_X);
    Xx = X(1:size-1);
    Xs = X(2:size);
    names = {'alpha', 'beta', 'sigma'};
    func = matlabFunction(CIR_Density);
    mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt, Xx, Xs)));
    theta_initial = [0.02, 0, 0.05];
    options = optimset('LargeScale','off');
    [theta_hat, fval, exitflag, output, grad,hessian] = fminunc(mleProb,theta_initial,options);
    n = length(theta_hat);
    % %estimate the sigma
     theta_Sig = sqrt(diag(inv(hessian/n)));
    % %p_value and confidence interval
     p_value_power = zeros(1,n);
     for i = 1:n
        z_score_power = abs((theta_hat(i)-theta_power(i))/(theta_Sig(i)*sqrt(n)));
        p_value_power(i) = 2*(1-normcdf(z_score_power));
        if p_value_power(i) < 0.05
            reject_power(i) = reject_power(i) + 1;
        end
     end        
 end
 
 %% Empirical Part with Fed Funds Rate
data = csvread('FRB_H15.csv',1, 1);
X_data = data/100;
len = length(data);
Xs = X_data(1:len-1);
Xx = X_data(2:len);
%-------------Using Linear Regression to Find Parameters as Theta0
x = X_data(1:end-1); % Time series of interest rates observations
dx = diff(X_data)./sqrt(x); %dx/sqrt(x)
regressors = [1./sqrt(x) sqrt(x)];
[coefficients, intervals, residuals] = ...
    regress(dx,regressors);
%Get the parameters
beta = - coefficients(2)/dt;
alpha = - coefficients(1)/coefficients(2);
sigma = std(residuals)/sqrt(dt);
InitialParams = [alpha, beta, sigma]; % Vector of initial parameters
%----------------------------------------------------------------
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.1, 0.4];
[theta_hat, fval, exitflag, output, grad, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i)-InitialParams(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n',names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));

end
%% Subsample Test 1
X_data = data/100;
len = length(data)/2;
Xs = X_data(1:len-1);
Xx = X_data(2:len);
%-------------Using Linear Regression to Find Parameters as Theta0
x = X_data(1:len-1); % Time series of interest rates observations
dx = diff(X_data(1:len))./sqrt(x); %dx/sqrt(x)
regressors = [1./sqrt(x) sqrt(x)];
[coefficients, intervals, residuals] = ...
    regress(dx,regressors);
%Get the parameters
beta = - coefficients(2)/dt;
alpha = - coefficients(1)/coefficients(2);
sigma = std(residuals)/sqrt(dt);
InitialParams = [alpha, beta, sigma]; % Vector of initial parameters
%----------------------------------------------------------------
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.1, 0.4];
[theta_hat, fval, exitflag, output, grad, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i)-InitialParams(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n',names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));
end
%% Subsample 2
X_data = data/100;
len = length(data)/2;
Xs = X_data(len+1:end-1);
Xx = X_data(len+2:end);
%-------------Using Linear Regression to Find Parameters as Theta0
x = X_data(len+1:end-1);%Time series of interest rates observations
dx = diff(X_data(len+1:end))./sqrt(x); %dx/sqrt(x)
regressors = [1./sqrt(x) sqrt(x)];
[coefficients, intervals, residuals] = ...
    regress(dx,regressors);
%Get the parameters
beta = - coefficients(2)/dt;
alpha = - coefficients(1)/coefficients(2);
sigma = std(residuals)/sqrt(dt);
InitialParams = [alpha, beta, sigma]; % Vector of initial parameters
%----------------------------------------------------------------
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.1, 0.3];
[theta_hat, fval, exitflag, output, grad, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i)-InitialParams(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n',names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));
end
 

%% Function to Simulate Data
function data = simulate(a, b, c, h, size, initial_value)
    X = zeros(size,1);
    X(1) = initial_value;
    for i = 2:size
        X(i) =b * (a- X(i-1)) *h + c * sqrt(X(i-1)*h)* normrnd(0,1) + X(i-1);
    end
    data = X;
end


%% Transition Density Function
function TDF = Density(muX, sigmaX, K, J)
    syms a  b  c
    syms xs ys zs
    syms x  y  z
    syms h t s 
    %Change X to Y
    fX2Y=int(1/sigmaX,x);
    fY2X=subs((finverse(fX2Y)), x,y);
    %Y's Drift and Diffusion
    muY_temp=muX/sigmaX-sym('1')/sym('2')*diff(sigmaX,x);
    muY=simplify(subs(muY_temp, x, fY2X));
    sigmaY=sym('1');
    
    %Change Y to Z
    fY2Z=h^(-1/2)*(y-ys);
    
    syms Htemp Expectation 
    sym Beta
    clear Beta Htemp Expectation 
    %slides 25 (7)
    for n=1:K
         HTemp=subs(Hermite(n), z, fY2Z);
         Expectation=HTemp;
         for k=1:J 
           HTemp=muY*diff(HTemp,y,1)+sym('1')/sym('2')*sigmaY*diff(HTemp, y, 2);
           %h = t - s
           Expectation=Expectation + h^k/factorial(k)*HTemp;
         end
         Beta{n}= sym('1')/factorial(n-1) * subs(Expectation, y, ys);
    end
    
    pZ=sym('0');
    for m=1:K
      pZ=pZ+Beta{m}*Hermite(m);
    end
    
    pZ=exp(-z^2/2)/sqrt(2*pi)*pZ;
    pY=(h^(-1/2))*subs(pZ, z, fY2Z);
    pX=(sigmaX^(-1))*subs(pY, y, fX2Y);
    pX=subs(pX, ys, subs(fX2Y, x, xs)) ;
    TDF=simplify(pX);
end

