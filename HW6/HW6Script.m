clear all
close all
clc

max_dimension = 20;
trials = 10^6;
n_vector = 1:max_dimension;
n_ball_volume_estimation = zeros(1, max_dimension);
n_ball_volume_exact = pi.^(n_vector / 2) ./ gamma(n_vector / 2 + 1);

for n = 1:max_dimension
    p = 2*rand(trials,n) - 1;
    radius = sqrt(sum(p.^2, 2));
    n_ball_volume_estimation(n) = 2^n * sum(radius <= 1) / trials;
end

plot(n_vector, n_ball_volume_estimation, '+k', n_vector, n_ball_volume_exact, 'go', 'linewidth', 1.5)
title('Monte Carlo Estimation of N-Ball Volume', 'interpreter', 'latex')
xlabel('$n$', 'interpreter', 'latex')
ylabel('$V_n$', 'interpreter', 'latex')
legend('approximation', 'exact', 'interpreter', 'latex')