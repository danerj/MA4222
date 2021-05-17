%% Section 1: Problem 1
% 1a
clear all
close all
clc

n = 8;
dx = 2*pi/n;
xeval = linspace(0,2*pi-dx, 2^8);
fx = @(x) -x.^2 + 2*pi*x;
x = (0:(2*pi/n):2*pi - dx)';
f = fx(x);
w = exp(2*pi*1i/n);
F_n = zeros(n,n);
for j = 0:n-1
    for k = 0:n-1
        F_n(j+1,k+1) = w^(j*k);
    end
end

c = F_n' * f / n;
z = exp(1i*xeval);
px = real(polyval(flip(c), z));

figure(1)
subplot(1,2,1)
plot(xeval,fx(xeval), 'm', xeval, px, 'g', x, f, 'kp', 'linewidth', 1)
title('Problem 1a: Signal and Fourier Interpolant', 'interpreter', 'latex')
xlabel('$x$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(x) = -x^2 - 2 \pi x$', '$p(x)$', '$f_j$'},...
    'interpreter', 'latex', 'location', 'south')

% 1b
m = n/2;
F_n = zeros(n,n);
for j = 0:n-1
    for k = 0:n-1
        F_n(j+1,k+1) = w^(j*(k-m));
    end
end

c = F_n' * f / n;
z = exp(1i*xeval);
px = 0;
for t = -m:m-1
    px = px + c(t+m+1)*z.^t;
end
px = real(px);

subplot(1,2,2)
plot(xeval,fx(xeval), 'm', xeval, px, 'g', x, f, 'kp', 'linewidth', 1)
title('Problem 1b: Signal and Low Frequency Fourier Interpolant', 'interpreter', 'latex')
xlabel('$x$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(x) = -x^2 - 2 \pi x$', '$p(x)$', '$f_j$'},...
    'interpreter', 'latex', 'location', 'south')

%% Section 2: Problem 2
% 2a
clear all

n = 256;
dt = 2*pi / n;
teval = linspace(0,2*pi, 2^10);
ft = @(t) sin(10*t) + sin(20*t);
t = (0:dt:2*pi-dt).';
f = ft(t);
data = f + 2*randn(size(t));

figure(2)
subplot(3,1,1)
plot(teval,ft(teval), 'k', 'linewidth', 1), hold on
plot(t, data, 'r', 'linewidth', 1)
title('Problem 2a: Signal and Noise', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(t)$', 'data'},...
    'interpreter', 'latex', 'location', 'northeast')

m = n/2;
w = exp(2*pi*1i/n);
F_n = zeros(n,n);
for j = 0:n-1
    for k = 0:n-1
        F_n(j+1,k+1) = w^(j*(k-m));
    end
end

c = F_n' * data / n;
z = exp(1i*teval);
pt = 0;
for kk = -m:m-1
    pt = pt + c(kk+m+1)*z.^kk;
end
pt = real(pt);

subplot(3,1,2)
plot(teval,ft(teval), 'k', teval, pt, 'g', 'linewidth', 1)
title('Problem 2a: Signal and Low Frequency Fourier Interpolant',...
    'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(t)$', '$p(t)$'},...
    'interpreter', 'latex', 'location', 'northeast')

% 2b
pt_truncated = 0;
l = 20;
for kk = -l:l
    pt_truncated = pt_truncated + c(kk+m+1)*z.^kk;
end
pt_truncated = real(pt_truncated);

subplot(3,1,3)
plot(teval,ft(teval), 'k', teval, pt_truncated, 'g', 'linewidth', 1)
title('Problem 2b: Signal and Truncated Low Frequency Fourier Interpolant',...
    'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(t)$', '$p(t)$'},...
    'interpreter', 'latex', 'location', 'northeast')

%% Section 3: Problem 2 (Alternative Method)
% Denoised by determining the largest Fourier coefficients.

clear all

n = 256;
dt = 2*pi / n;
teval = linspace(0,2*pi, 2^10);
ft = @(t) sin(10*t) + sin(20*t);
t = (0:dt:2*pi-dt).';
f = ft(t);
data = f + 2*randn(size(t));

figure(3)
subplot(2,2,1)
plot(teval,ft(teval), 'k', 'linewidth', 1), hold on
plot(t, data, 'r', 'linewidth', .5)
xlim([0 2*pi])
ylim([-5 5])
title('Denoising Data with the FFT: Signal and Noise',...
    'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(t)$', 'data'},...
    'interpreter', 'latex', 'location', 'northeast')

c = fft(data,n);
PSD = c.*conj(c)/ n;
frequency = 2*pi /(dt*n) * (0:n);
L = 1:floor(n/2);
subplot(2,2,3)
plot(frequency(L),PSD(L), 'c', 'linewidth', 1)
title('Fourier Coefficient Sizes', 'interpreter', 'latex')
xlabel('$k$', 'interpreter', 'latex')
ylabel('$|c_k|^2 / n $', 'interpreter', 'latex')

indices = PSD > max(PSD) / 3;
PSD_clean = PSD.*indices;
c = indices.*c;
f_filtered = ifft(c);
subplot(2,2,4)
plot(frequency(L),PSD_clean(L), 'c', 'linewidth', 1)
title('Filtered Fourier Coefficient Sizes', 'interpreter', 'latex')
xlabel('$k$', 'interpreter', 'latex')
ylabel('$|c_k|^2 / n $', 'interpreter', 'latex')

subplot(2,2,2)
plot(teval, ft(teval), 'k', t, f_filtered, 'g', 'linewidth', 1)
title('Denoising Data with the FFT: Signal and Recovered Signal',...
    'interpreter', 'latex')
xlim([0 2*pi])
ylim([-5 5])
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(t)$', 'Recovered $f(t)$'},...
    'interpreter', 'latex', 'location', 'southeast')

%% Section 4: Testing my FFT algorithm by comparing to Strang problems

clear all

% Should produce y = (4,0,0,0,4,0,0,0)^T to test homemadefft.m
c = [1 0 1 0 1 0 1 0].';
y = homemadefft(c)
% Reproduce the c given above to test homemadeifft.m
c = homemadeifft(y)

% Should produce y = (4,0,0,0,-4,0,0,0)^T to test homemadefft.m
c = [0 1 0 1 0 1 0 1].';
y = homemadefft(c)
% Reproduce the c given above to test homemadeifft.m
c = homemadeifft(y)

% Should produce y = (20,-4,-4,-4)^T to test homemadefft.m
c = [2 6 6 6].';
y = homemadefft(c)
% Reproduce the c given above to test homemadeifft.m
c = homemadeifft(y)

% Should produce y = (0,-4i,0,4i)^T to test homemadefft.m
c = [0 -2 0 2].';
y = homemadefft(c)
% Reproduce the c given above to test homemadeifft.m
c = homemadeifft(y)

% Should produce 
% y = (20,-4-4iw_8,-4,-4+4iw_8^3,20,-4+4iw_8,-4,-4-4iw_8^3) 
% Since w_8 produces terms that are more complicated than the previous test
% calculations, use an additional piece of code to check the function.
c = [2 0 6 -2 6 0 6 2].';
y = homemadefft(c)
w_8 = exp(2*pi*1i/8);
y_check = [20, -4-4i*w_8, -4, -4+4i*w_8^3, 20, -4+4i*w_8, -4, -4-4i*w_8^3].'
% Reproduce the c given above to test homemadeifft.m
c = homemadeifft(y)

%% Section 5: Testing my FFT algorithm to see if it's actually 'Fast'.
% Testing my FFT against DFT (by standard matrix multiplication) and
% Matlab's built in FFT function (using the function ifft since Matlab uses
% the opposite convention for defining the FFT vs. inverse FFT - Matlab
% uses what is often referred to as 'omega_n = exp(-2ipi/n)' instead of
% using 'w_n = exp(2ipi/n)'. So to swap conventions you swap fft with ifft.

clear all

n = 2^10;
c = randi([-10,10], n, 1);
w = exp(2*pi*1i/n);
F_n = zeros(n,n);
for j = 0:n-1
    for k = 0:n-1
        F_n(j+1,k+1) = w^(j*k);
    end
end

tic
yfft = homemadefft(c);
toc

tic
ydft = F_n*c;
toc

norm(yfft-ydft)

tic
y_matlab_fft = ifft(c);
toc

n_vector = [2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10, 2^11, 2^12];
computation_times = zeros(1,length(n_vector));
for j = 1:length(n_vector)
    n = n_vector(j);
    c = randi([-10,10], n, 1);
    tic
    y = homemadefft(c);
    computation_times(j) = toc;
end

figure(4)
plot(n_vector, computation_times, 'g--', 'linewidth', 2)
title('Homemade FFT Computation Time', 'interpreter', 'latex')
xlabel('$n$', 'interpreter', 'latex')
ylabel('$t$ (seconds)', 'interpreter', 'latex')

% Conclusion: I think my algorithm is correct and does appear to have
% loglinear time complexity (it appears linear so I believe I have
% successfully avoided quadratic time complexity). 

%% Section 6: Problem 2 (Alternative Method using Homemade FFT and IFFT)
% Denoised by determining the largest Fourier coefficients.

clear all

n = 256;
dt = 2*pi / n;
teval = linspace(0,2*pi, 2^10);
ft = @(t) sin(10*t) + sin(20*t);
t = (0:dt:2*pi-dt).';
f = ft(t);
data = f + 2*randn(size(t));

figure(5)
subplot(2,2,1)
plot(teval,ft(teval), 'k', 'linewidth', 1), hold on
plot(t, data, 'r', 'linewidth', .5)
xlim([0 2*pi])
ylim([-5 5])
title('Denoising Data with the FFT: Signal and Noise',...
    'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(t)$', 'data'},...
    'interpreter', 'latex', 'location', 'northeast')

c = homemadefft(data);
PSD = c.*conj(c)/ n;
subplot(2,2,3)
plot(1:n,PSD, 'c', 'linewidth', 1)
title('Fourier Coefficient Sizes', 'interpreter', 'latex')
xlabel('$k$', 'interpreter', 'latex')
ylabel('$|c_k|^2 / n $', 'interpreter', 'latex')

indices = PSD > max(PSD) / 3;
PSD_clean = PSD.*indices;
c = indices.*c;
f_filtered = homemadeifft(c);
subplot(2,2,4)
plot(1:n,PSD_clean, 'c', 'linewidth', 1)
title('Filtered Fourier Coefficient Sizes', 'interpreter', 'latex')
xlabel('$k$', 'interpreter', 'latex')
ylabel('$|c_k|^2 / n $', 'interpreter', 'latex')

subplot(2,2,2)
plot(teval, ft(teval), 'k', t, f_filtered, 'g', 'linewidth', 1)
title('Denoising Data with the FFT: Signal and Recovered Signal',...
    'interpreter', 'latex')
xlim([0 2*pi])
ylim([-5 5])
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend({'$f(t)$', 'Recovered $f(t)$'},...
    'interpreter', 'latex', 'location', 'southeast')

%% Section 7: Denoising a Music file

clear all

[y, Fs] = audioread('Chopin.mp3');

%Using clip of the audio file.
y = y(10^6:10^6+5*10^5-1,:);

%Create a mono version of the clip.
y = mean(y,2);
n = length(y);


%To add noise, uncomment the line below and adjust noise level.
data = y + 10^-2*randn(n,1);

%To play the distorted clip:
%sound(data,Fs)

c = fft(data);
PSD = c.*conj(c)/ n;
figure(6)
plot(1:n, PSD)
title('Power Spectral Density Noisy Music', 'interpreter', 'latex')

indices = PSD > .001;
PSD_clean = PSD.*indices;
c = indices.*c;
y_filtered = ifft(c);

sound(y_filtered, Fs)