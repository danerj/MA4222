% Script file for comparing different solvers for the west0479 test
% matrix of the Harwell/Boing sparse matrix test set.

data = load('west0479');
A = data.west0479;
b = ones(size(A,1),1);

methods = {'CGS', 'BI-CG'};
nmethods = length(methods);
linestyles = {'rs-','bo-'};
domethods = 1:2;

tol = 1e-6;
maxit = 450;     

for mk = 1:nmethods;
  m = domethods(mk);
  disp(['Method: ' methods{m} ', maxit = ', num2str(maxit)]);
  
  switch m
%    case 1
%     % CGS
%     [x, flag, res, it, resvec] = cgs(A,b,tol,maxit);
%    case 2
%     % BI-CG
%     [x, flag, res, it, resvec] = bicg(A,b,tol,maxit);
      case 1
          [x, flag, res, it, resvec] = gmres(A,b,10,tol,maxit);
      case 2
          % bicgstab
          [x, flag, res, it, resvec] = bicgstab(A,b,tol,maxit);
  end
  disp(['Computed relative residual: ' num2str(res)]);
      
  h = semilogy(0:(length(resvec)-1),resvec/norm(b), linestyles{mk}, ...
             'LineWidth', 1);
  hold on;
end
hold off;
xlabel('number of iterations');
ylabel('relative residual');
legend(methods(domethods),1);

%%
clear all
close all
clc

maxit = 10^4;
tolerance = 10^-4;

for n = [5, 20, 100, 500, 1000]

A=toeplitz([2, -1, zeros(1, n-2)])*(n+1)^2;
z=(1:n)'/(n+1);
b=z.*sin(z);
tic
[x_cgs, flag_cgs, res_cgs, it_cgs, resvec_cgs] = cgs(A,b,tolerance,maxit);
toc
tic
[x_bicg, flag_bicg, res_bicg, it_bicg, resvec_bicg] = bicg(A,b,tolerance,maxit);
toc
tic
L = ichol(A);
spy(L)
[x, flag, res, it, resvec] = gmres(A,b,[],tol,maxit, L, L');
toc

end

plot(resvec_cgs)
figure(2)
plot(resvec_bicg)

