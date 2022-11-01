%% Example using LAR and LASSO on the diabetes data set

clear; close all; clc;
addpath('../SpaSM');

%% Diabetes data set, load and extract variables
filename = 'diabetes.txt';
urlwrite('https://www.stanford.edu/~hastie/Papers/LARS/diabetes.data',filename);
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
diabetes = table(dataArray{1:end-1}, 'VariableNames', {'AGE','SEX','BMI','BP','S1','S2','S3','S4','S5','S6','Y'});
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% Delete the file
delete diabetes.txt

predNames = diabetes.Properties.VariableNames(1:end-1);
X = diabetes{:,1:end-1};
y = diabetes{:,end};

%% Normalize
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
h1 = figure(2);
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
