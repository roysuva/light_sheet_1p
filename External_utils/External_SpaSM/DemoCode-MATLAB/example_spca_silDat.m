%% Example using SPCA and PCA on the silhouette data set

clear; close all; clc;
addpath('../SpaSM');

%% Silhouette data set, load and extract variables
load '../SpaSM/Silhouettes.mat'

% Initial look at the data
plot(Xa(:,1:65)',Xa(:,66:end)','-');
axis equal;

% Tru to center data and then plot it right
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
devVec = [-2,2,0];
devCol = ['g','b','r'];
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
        devShape = centX'+devVec(j)*sqrt(SD(i))*SL(:,i);
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


