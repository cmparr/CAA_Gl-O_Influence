function [h, avg, lens] = PlotErrorEllipse(data, varargin)
%%%%%%%%
% add covariance ellipses to PCA plots. 
% adapted from https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
% adapted by Claire Parrott (cparrott@eoas.ubc.ca)
%
% INPUT: data Mx2 matrix
%        Should be PC1 and PC2 of categorization subset from input PCAs
%
% varargin:
%       'median' -  use median value as average center of ellipse
%       'color'  -  Set color of ellipse. Pass color variable after color.
%       'save'   -  Write eigenvector size/angles to 
%       'cond'   -  position n+1 is condition to use as tab name in writing
%
% OUTPUT: h    line handle for ellipse
%         avg  
%
%
%%%%%%%%

% Calculate the eigenvectors and eigenvalues of data covariance
[eigenvec, eigenval] = eig(cov(data));

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
else
    smallest_eigenval = max(eigenval(:,1));
end

% Calculate the angle between the x-axis and the largest eigenvector
phi = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(phi < 0)
    phi = phi + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);
if any(strcmp(varargin, 'median'))
avg = median(data);
end
% Get the 95% confidence interval error ellipse
chisquare_val = sqrt(icdf('Chisquare', 0.9, 2)); % root of 2 DOF interval from chi-squared test
theta_grid = linspace(0,2*pi);
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);
lens = [a b];

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
col = 'k';
if any(strcmp(varargin,'color'))
    col = varargin{find(strcmp(varargin, 'color'))+1};
end
h = plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-','HandleVisibility','off', 'Color',col, 'Linewidth', 2);
hold on;

% save 
if any(strcmp(varargin, 'save'))&any(strcmp(varargin,'cond'))
    
% output eigenvector values to look at numerics
sheetdata.eigv1 = largest_eigenvec(1);
sheetdata.eigv2 = largest_eigenvec(2);
sheetdata.angle = rad2deg(phi);
table = struct2table(sheetdata);
cd('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\Combined_Nuts_CTD\');
fname = varargin{find(strcmp(varargin, 'cond'))+1};
fname = regexprep(fname,' ','');

writetable(table, 'Ellipses-FWmetrics', 'FileType', 'spreadsheet', 'Sheet', fname); 
end
end