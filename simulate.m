function data = simulate(a, b, c, h, size, initial_value)
    X = zeros(size,1);
    X(1) = initial_value;
    for i = 2:size
        X(i) =b*(a - X(i-1))*h + X(i-1) + c * sqrt(X(i-1)*h)* normrnd(0,1);
    end
    data = X;
end