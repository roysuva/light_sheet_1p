%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo for SpaSM Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file includes all the examples from the paper:
% SpaSM: A MATLAB Toolbox for Sparse Statistical Modeling
% 
% It is recommended that the user runs the code chunks 
% separated by %% to inspect the associated output and
% figures for a given example. Individual examples are
% also included separate files. The MATLAB console is 
% cleared between examples and figures are closed.

% add path to the MATLAB toolbox SpaSM.
addpath('../SpaSM');

%---------------------------------------------------------------------------
%% Forward selection
%---------------------------------------------------------------------------

% Fix stream of random numbers
s1 = RandStream.create('mrg32k3a','Seed', 22);
s0 = RandStream.setGlobalStream(s1);

% Create data set
n = 100; p = 6;
correlation = 0.6;
Sigma = correlation*ones(p) + (1 - correlation)*eye(p);
mu = zeros(p,1);
X = mvnrnd(mu, Sigma, n);
% Model is lin.comb. of first three variables plus noise
y = X(:,1) + X(:,2) + X(:,3) + 2*randn(n,1);

% Preprocess data
X = normalize(X);
y = center(y);

% Run FORWARDSELECTION
[beta info] = forwardselection(X, y, 0, true, true);

% Inspect the information criteria output
info

% Find best fitting model
[bestAIC bestIdx] = min(info.AIC);
best_s = info.s(bestIdx);

% Plot results
h1 = figure(1);
plot(info.s, beta, '.-');
xlabel('s'), ylabel('\beta', 'Rotation', 0)
line([best_s best_s], [-6 14], 'LineStyle', ':', 'Color', [1 0 0]);
xlim([0 1]);
legend('1','2','3','4','5','6');
title('Magnitude of coefficients using forward selection')

%---------------------------------------------------------------------------
%% Least Angle Regression
%---------------------------------------------------------------------------

% Clear the workspace
clear; close all; clc;

% Fix stream of random numbers
s1 = RandStream.create('mrg32k3a','Seed', 22);
s0 = RandStream.setGlobalStream(s1);

% Create data set
n = 100; p = 6;
correlation = 0.6;
Sigma = correlation*ones(p) + (1 - correlation)*eye(p);
mu = zeros(p,1);
X = mvnrnd(mu, Sigma, n);
% Model is lin.comb. of first three variables plus noise
y = X(:,1) + X(:,2) + X(:,3) + 2*randn(n,1);

% Preprocess data
X = normalize(X);
y = center(y);

% Run LAR
[beta info] = lar(X, y, 0, true, true);

% Inspect the information criteria output
info

% Find best fitting model
[bestAIC bestIdx] = min(info.AIC);
best_s = info.s(bestIdx);

% Plot results
h1 = figure(1);
plot(info.s, beta, '.-');
xlabel('s'), ylabel('\beta', 'Rotation', 0)
line([best_s best_s], [-6 14], 'LineStyle', ':', 'Color', [1 0 0]);
xlim([0 1]);
legend('1','2','3','4','5','6');
title('Magnitude of coefficients using least angle regression')

%---------------------------------------------------------------------------
%% LASSO
%---------------------------------------------------------------------------

clear; close all; clc;

% Fix stream of random numbers
s1 = RandStream.create('mrg32k3a','Seed', 22);
s0 = RandStream.setGlobalStream(s1);

% Create data set
n = 100; p = 6;
correlation = 0.6;
Sigma = correlation*ones(p) + (1 - correlation)*eye(p);
mu = zeros(p,1);
X = mvnrnd(mu, Sigma, n);
% Model is lin.comb. of first three variables plus noise
y = X(:,1) + X(:,2) + X(:,3) + 2*randn(n,1);

% Preprocess data
X = normalize(X);
y = center(y);

% Run LASSO
[beta info] = lasso(X, y, 0, true, true);

% Inspect the information criteria output
info

% Find best fitting model
[bestAIC bestIdx] = min(info.AIC);
best_s = info.s(bestIdx);

% Plot results
h1 = figure(1);
plot(info.s, beta, '.-');
xlabel('s'), ylabel('\beta', 'Rotation', 0)
line([best_s best_s], [-6 14], 'LineStyle', ':', 'Color', [1 0 0]);
xlim([0 1]);
legend('1','2','3','4','5','6');
title('Magnitude of coefficients using lasso regression')

%---------------------------------------------------------------------------
%% Elastic Net
%---------------------------------------------------------------------------

clear; close all; clc;

% Fix stream of random numbers
s1 = RandStream.create('mrg32k3a','Seed', 42);
s0 = RandStream.setGlobalStream(s1);

% Create data set
n = 30; p = 40;
correlation = 0.1;
Sigma = correlation*ones(p) + (1 - correlation)*eye(p);
mu = zeros(p,1);
X = mvnrnd(mu, Sigma, n);
% Model is lin.comb. of first three variables plus noise
y = X(:,1) + X(:,2) + X(:,3) + 0.5*randn(n,1);

% Preprocess data
X = normalize(X);
y = center(y);

% Run elastic net
delta = 1e-3;
[beta info] = elasticnet(X, y, delta, 0, true, true);

% Plot results
h1 = figure(1);
plot(info.s, beta, '.-');
xlabel('s'), ylabel('\beta', 'Rotation', 0)
xlim([0,1])
title('Magnitude of coefficients using elastic net')

%---------------------------------------------------------------------------
%% Sparse Principal Component Analysis
%---------------------------------------------------------------------------

clear; close all; clc;

% Fix stream of random numbers
s1 = RandStream.create('mrg32k3a','Seed', 11);
s0 = RandStream.setGlobalStream(s1);

% Create synthetic data set
n = 1500; p = 500;
t = linspace(0, 1, p);
pc1 = max(0, (t - 0.5)> 0);
pc2 = 0.8*exp(-(t - 0.5).^2/5e-3);
pc3 = 0.4*exp(-(t - 0.15).^2/1e-3) + 0.4*exp(-(t - 0.85).^2/1e-3);
X = [ones(n/3,1)*pc1 + randn(n/3,p); ones(n/3,1)*pc2 + randn(n/3,p);...
  ones(n/3,1)*pc3 + randn(n/3,p)];

% PCA and SPCA
[U D V] = svd(X, 'econ');
d = sqrt(diag(D).^2/n);
K = 3;
delta = inf;
stop = -[250 125 100];
maxiter = 3000;
convCriterion = 1e-9;
verbose = true;

[SL SD] = spca(X, [], K, delta, stop, maxiter, convCriterion, verbose);

figure(1)
plot(t, [pc1; pc2; pc3]); axis([0 1 -1.2 1.2]);
title('Noiseless data');
legend('pc1','pc2','pc3','Location','NorthWest');
ylim([-0.2,1.1])
figure(2);
plot(t, X([1:5,501:505,1001:1005],:));  axis([0 1 -6 6]);
title('Data + noise');
figure(3);
plot(t, d(1:3)*ones(1,p).*(V(:,1:3)'));  axis([0 1 -1.2 1.2]);
title('PCA');
legend('pc1','pc2','pc3','Location','NorthWest');
ylim([-0.2,0.8])
figure(4)
plot(t, sqrt(SD)*ones(1,p).*(SL'));  axis([0 1 -1.2 1.2]);
title('SPCA');
legend('pc1','pc2','pc3','Location','NorthWest');
ylim([-0.2,0.8])

%---------------------------------------------------------------------------
%% Sparse Linear Discriminant Analysis
%---------------------------------------------------------------------------

clear; close all; clc;

% Fix stream of random numbers
s1 = RandStream.create('mrg32k3a','Seed', 50);
s0 = RandStream.setGlobalStream(s1);

p = 150; % number of variables
nc1 = 100; % number of observations per class
nc2 = 100; % number of observations per class
nc3 = 100; % number of observations per class
n = nc1+nc2+nc3; % total number of observations
m1 = 0.6*[ones(10,1); zeros(p-10,1)]; % c1 mean
m2 = 0.6*[zeros(10,1); ones(10,1); zeros(p-20,1)]; % c2 mean
m3 = 0.6*[zeros(20,1); ones(10,1); zeros(p-30,1)]; % c3 mean
S = 0.6*ones(p) + 0.4*eye(p); % covariance is 0.6

% training data
c1 = mvnrnd(m1,S,nc1); % class 1 data
c2 = mvnrnd(m2,S,nc2); % class 2 data
c3 = mvnrnd(m3,S,nc3); % class 3 data
X = [c1; c2; c3]; % training data set
Y = [[ones(nc1,1);zeros(nc2+nc3,1)] [zeros(nc1,1); ones(nc2,1); zeros(nc3,1)] [zeros(nc1+nc2,1); ones(nc3,1)]];

% test data
c1 = mvnrnd(m1,S,nc1);
c2 = mvnrnd(m2,S,nc2);
c3 = mvnrnd(m3,S,nc3);
X_test = [c1; c2; c3];

% SLDA parameters
delta = 1e-3; % l2-norm constraint
stop = -30; % request 30 non-zero variables
maxiter = 250; % maximum number of iterations
Q = 2; % request two discriminative directions
tol = 1e-6;

% normalize training and test data
[X mu d] = normalize(X);
X_test = (X_test-ones(n,1)*mu)./sqrt(ones(n,1)*d);

% run SLDA
[B theta] = slda(X, Y, delta, stop, Q, maxiter, tol, true);

% Project data onto the sparse directions
DC = X*B;
DC_test = X_test*B;

% Classification (LDA of projected data)
Yc = [ones(nc1,1); 2*ones(nc2,1); 3*ones(nc3,1)];
[class err] = classify(DC, DC, Yc, 'linear');
[class_test] = classify(DC_test, DC, Yc, 'linear');
err_test = sum(Yc ~= class_test)/length(Yc);
fprintf('SLDA result: training error is %2.1f %%, test error is %2.1f %%.\n', 100*err, 100*err_test);

[class err] = classify(X, X, Yc, 'linear');
[class_test] = classify(X_test, X, Yc, 'linear');
err_test = sum(Yc ~= class_test)/length(Yc);
fprintf('LDA result: training error is %2.1f %%, test error is %2.1f %%.\n', 100*err, 100*err_test);
 
% plot sparse discriminative directions for test data
figure;
plot(DC_test(1:nc1,1), DC_test(1:nc1,2),'ro'), hold on
plot(DC_test((nc1+1):(nc1+nc2),1), DC_test((nc1+1):(nc1+nc2),2),'ks')
plot(DC_test(((nc1+nc2)+1):(nc1+nc2+nc3),1), DC_test(((nc1+nc2)+1):(nc1+nc2+nc3),2),'bv')
xlabel('1st direction'), ylabel('2nd direction')
legend('C_1','C_2','C_3','Location','SouthEast')

%---------------------------------------------------------------------------
%% Regression on diabetic data set
%---------------------------------------------------------------------------
% This example illustrates the difference between the LAR and LASSO
% Three code chunks separated by %% are included in this example
% (1) The data is downloaded and preprocessed
% (2) An example with LAR is given
% (3) An example with the LASSO is given
% Observe the difference in the plots when a coefficient crosses zero.

clear; close all; clc;

% Download the data
filename = 'diabetes.txt';
urlwrite('https://www.stanford.edu/~hastie/Papers/LARS/diabetes.data',filename);
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
diabetes = table(dataArray{1:end-1}, 'VariableNames', {'AGE','SEX','BMI','BP','S1','S2','S3','S4','S5','S6','Y'});
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% Cleanup of downloaded data
delete diabetes.txt

predNames = diabetes.Properties.VariableNames(1:end-1);
X = diabetes{:,1:end-1};
y = diabetes{:,end};

X = normalize(X);
y = center(y);

%% Example with LAR

[beta info] = lar(X, y, 0, true, true);

% Find best fitting model
[bestAIC bestIdx] = min(info.AIC);
best_s = info.s(bestIdx);

% Plot results
lsty = {'.-r','+-g','o-b','*-c','s-m','d-y','^-k','v-r','>-g','<-b'};
h1 = figure(1);
set(h1, 'Position', [100 100 800 500])
hold on
for i = 1:size(beta,1)
    plot(info.s, beta(i,:), char(lsty(i)));
end
xlabel('s'), ylabel('\beta', 'Rotation', 0)
line([best_s best_s], [-600 600], 'LineStyle', ':', 'Color', [1 0 0]);
legend('1','2','3','4','5','6','7','8','9','10','Location','NorthWest');
title('Magnitude of coefficients using least angle regression on Diabetes data set')

%% Example with LASSO

[beta info] = lasso(X, y, 0, true, true);

% Find best fitting model
[bestAIC bestIdx] = min(info.AIC);
best_s = info.s(bestIdx);

% Plot results
lsty = {'.-r','+-g','o-b','*-c','s-m','d-y','^-k','v-r','>-g','<-b'};
h1 = figure(1);
set(h1, 'Position', [100 100 800 500])
hold on
for i = 1:size(beta,1)
    plot(info.s, beta(i,:), char(lsty(i)));
end
xlabel('s'), ylabel('\beta', 'Rotation', 0)
xlim([0 1])
line([best_s best_s], [-600 600], 'LineStyle', ':', 'Color', [1 0 0]);
legend('1','2','3','4','5','6','7','8','9','10','Location','NorthWest');
title('Magnitude of coefficients using lasso regression on Diabetes data set')

%---------------------------------------------------------------------------
%% PCA and SPCA on silhouette data
%---------------------------------------------------------------------------
% This example illustrates the difference between PCA and SPCA
% Two code chunks separated by %% are included in this example
% (1) The data is loaded and preprocessed
% (2) PCA and SPCA are run on the data and figures produced
% Observe the difference in the plots, SPCA gives more local 
% descriptors in the shapes compared to PCA

clear; close all; clc;

% Silhouette data set, load and extract variables
load './Silhouettes.mat'

% Initial look at the data
plot(Xa(:,1:65)',Xa(:,66:end)','-')
axis equal

% Center data and then plot it right
[X centX] = center(Xa);
XX = X + ones(39,1)*centX;
plot(XX(:,1:65)',XX(:,66:end)','-')
axis equal

%% PCA and SPCA

[coeff, scores, latent, ~, explained, mu] = pca(Xa);

figure;
imagesc(cov(Xa));
colormap(gray)
title('Covariance matrix of data')

% Mean shape and deviations from mean shape with PCA

% Vector of deviations
devVec = [-2,2,0]; % Deviances from mean along PCs
devCol = ['g','b','r']; % Colors for deviances from mean and mean
h1 = figure;
set(h1, 'Position', [100 100 800 500])
numComp = 3;
for i=1:numComp
    subplot(1,numComp,i)
    hold on
    for j=1:3
        % Deviation shape
        devShape = mu + (devVec(j)*sqrt(latent(i))*coeff(:,i))';
        % Reshape it for plotting
        PlShape = [devShape(66:end);devShape(1:65)]';
        DrawShape(PlShape,devCol(j));
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        axis equal
    end
    title(['PC ',num2str(i),', Deviation \pm ', num2str(devVec(2)), '\surd\lambda']);
    %legend('-2\surd\lambda','+2\surd\lambda','mean',2,'Location','NorthWest');
    hold off
end


% Same for spca, begin by centering it
[X centX] = center(Xa);
% Define the parameters for spca
K = 3;
delta = inf;
stop = -[70,60,40];
maxiter = 3000;
convergenceCriterion = 1e-9;
verbose = true;
[SL SD] = spca(X, [], K, delta, stop, maxiter, convergenceCriterion, verbose);

devVec = [-2,2,0];
devCol = ['g','b','r'];
h1 = figure;
set(h1, 'Position', [100 100 800 500])
for i=1:numComp
    subplot(1,numComp,i)
    hold on
    for j=1:3
        % Deviation shape
        devShape = mu'+devVec(j)*sqrt(SD(i))*SL(:,i);
        % Reshape it for plotting
        PlShape = [devShape(66:end)';devShape(1:65)']';
        DrawShape(PlShape,devCol(j));
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        axis equal
    end
    title(['SPC ',num2str(i),', Deviation \pm ', num2str(devVec(2)), '\surd\lambda']);
    %legend('-2\surd\lambda','+2\surd\lambda','mean',2,'Location','NorthWest');
    hold off
end

