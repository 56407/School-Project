% Final Project Part: MLE for drift
% BEIXI LEI, ZETIAN WU, JUNYANG NIU

% Mydata as 252*68 table (returns) opened in csv file generated from R
% volatility as 68*1 table

data = table2array(Mydata);
vol = table2array(volatility);

f = @(x,v,theta) (1./(sqrt(2*pi)*v)*exp(-(x-theta).^2)/(2*v^2));

drift = zeros(1,size(data,2));
for i = 1:size(data,2)
    X = data(:,i);
    L = @(theta) -sum(log(f(X,vol(i),theta))); 
    [theta_hat, ~] = fminsearch(L,0);
    drift(i) = theta_hat;
end

drift = drift';

csvwrite('drift.csv',drift);