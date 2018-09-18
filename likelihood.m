function llfunc = likelihood(alpha, Volatility,dt, Y)
    llfunc = 0;
    Y_avg = mean(Y);
    for i = 2:length(Y)
        mu = Y(i-1)*exp(-alpha*dt)+(1-exp(-alpha*dt))*Y_avg;
        sigma = Volatility*sqrt((1-exp(-2*alpha*dt))/(2*alpha));
        llfunc = llfunc+log(normpdf(Y(i),mu,sigma));
    end
end