alpha = 0.1;
beta = 0.02;
sigma = 0.02;
dt = 1/252;
size = 600;
initial_X = 0.01;
X = simulate(alpha, beta, sigma, dt, size, initial_X);

