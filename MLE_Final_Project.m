%% Size Analysis 
 Test_Times = 1000;
 reject = zeors(3,1);
 for path = 1:Test_Times
    alpha = 0.01;
    beta = -0.03;
    sigma = 0.04;
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
    theta_initial = [0.2, -0.1, 0.4];
    options = optimset('LargeScale','off');
    [theta_hat, fval, exitflag, output, grad,hessian] = fminunc(mleProb,theta_initial,options);
    n = length(theta_hat);
    % %estimate the sigma
     theta_Sig = sqrt(diag(inv(hessian/n)));
    % %p_value and confidence interval
     p_value = zeros(1,n);
     for i = 1:n
        z_score = abs((theta_hat(i)-theta(i))/(theta_Sig(i)*sqrt(n)));
        p_value(i) = 2*(1-normcdf(z_score));
        if p_value(i) < 0.05
            reject(i) = reject(i) + 1;
        end
     end        
 end