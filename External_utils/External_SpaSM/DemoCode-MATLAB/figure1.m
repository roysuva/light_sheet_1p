clear; close all; clc;
addpath('../SpaSM');
 
% Create synthetic data set
maxiter = 1000;
convergenceCriterion = 1e-6;
verbose = false;
stop = -25;
 
K = 3;
delta = 1e-2;
n = 600; p = 10:10:1500;
np = length(p);
times = zeros(2,np);
variances = zeros(2,np);
for i = 1:np
  t = linspace(0, 1, p(i));
  pc1 = max(0, (t - 0.5)> 0);
  pc2 = 0.8*exp(-(t - 0.5).^2/5e-3);
  pc3 = 0.4*exp(-(t - 0.15).^2/1e-3) + 0.4*exp(-(t - 0.85).^2/1e-3);
  X = [ones(n/3,1)*pc1 + randn(n/3,p(i)); ones(n/3,1)*pc2 + randn(n/3,p(i));...
    ones(n/3,1)*pc3 + randn(n/3,p(i))];
 
  disp(i)
   
  tstart = tic;
  [sv sd] = spca(X, [], K, delta, stop, maxiter, convergenceCriterion, verbose);
  times(1,i) = toc(tstart);
  variances(1,i) = sum(sd);
   
  tstart = tic;
  [sv sd] = spca_zouhastie(X, [], K, delta, stop, maxiter, convergenceCriterion, verbose);
  times(2,i) = toc(tstart);
  variances(2,i) = sum(sd);
   
end
 
figure(1)
semilogy(p, times, '*');
xlabel(['# variables (' num2str(n) ' observations)']);
ylabel('time');
legend('Sequential SPCA', 'Simultaneous SPCA', 'Location', 'SE');
 
figure(2)
plot(p, variances, '*');
xlabel(['# variables (' num2str(n) ' observations)']);
ylabel('variance explained');
legend('Sequential SPCA', 'Simultaneous SPCA', 'Location', 'SE');
