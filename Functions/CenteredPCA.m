function [eigenvectors, PCs, eigval, variance, rowidx] = CenteredPCA(data, variables, varargin)
%
%   data:       MXN doubles
%               Numeric only data with N variables in columns and M observations
%   
%   variables:  N char array
%               Char array of input variable names
%
%   varargin:   'plot'      include plots of eigenvectors and first 3 PCs

% eliminate rows with nan values in any of the fields
rowidx= any(isnan(data'));
data = data(~rowidx,:);

% standardize data and center around standardized mean
mean_d=mean(data);
std_d= std(data);

n=size(data,1);
data_mean=repmat(mean_d,n,1);
data_std= repmat(std_d,n,1);

dataszd=(data-data_mean)./data_std; % standardize data

% perform PCA
[eigenvectors,PCs,eigval]=pca(dataszd);

% contribution of each mode to total variance:
variance=eigval./sum(eigval);

if any(strcmp(varargin, 'plot'))
figure; 
plot([1:length(variance)],variance,'*-'); hold on;
text([1:length(variance)] + 0.1 ,variance + 0.025, num2str(round(variance,3)), 'HorizontalAlignment','left');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% plot the first 3 modes (eigenvectors and PCs)
% gives length of eigenvector projected onto PCX space, amount of
% covariance with that mode of variability
f = figure; f.Position(3:4) = [800 600];
for i=1:3
subplot(2,3,i)
plot([1:length(eigenvectors(:,1))],eigenvectors(:,i),'o-'); 
title(['e_',num2str(i)]);
xlabel('variables'); ylabel('Eigenvectors');
xlim([0 length(eigenvectors(:,1))]); ylim([-0.70 0.70]);
set(gca,'XTick',[1:1:length(eigenvectors)],'XTicklabel',variables)

subplot(2,3,i+3)
plot(PCs(:,i));
xlabel('observations'); 
title(['PC',num2str(i)]);
end
%
figure;
scatter(PCs(:,1), PCs(:,2),'filled'); hold on;
xlabel(strcat('PC1 (', num2str(variance(1)*100,'%.2f'),')')); ylabel(strcat('PC2 (', num2str(variance(2)*100,'%.2f'),')'));
%}
end
end