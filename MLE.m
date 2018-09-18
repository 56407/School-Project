[raw,~,~] = xlsread('interest rate.xls');
Price = raw(1:end,3)/100;
Return = price2ret(Price);

theta0 = [0.012, 0.025];
objfunc = @(theta)-likelihood(theta(1),theta(2), Price);
[thetahat,fval,exitflag,output,fval1,fval2] = fminunc(objfunc,theta0);

theta0 = [0.01, 0.02];
objfunc = @(theta)-likelihood(theta(1),theta(2),1/252, Price);
A = [];
b = [];
Aeq = [];
Ceq = [];
lb = [0 0];
ub = [1000 10000];
[thetahat,fval,~,~,~,~,fval2] = fmincon(objfunc,theta0, A, b, Aeq, Ceq, lb, ub);