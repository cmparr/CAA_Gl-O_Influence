function [p, unique_pairs_name] = perMANOVA(data_all, groups_all)
% compute non-parametric manova p values for multivariate dataset
% 
% INPUT:    data_in = MxN double
%           matrix of "raw" numeric values. each row =  1 sample
%
%           groups = Mx1 categorical 
%           list of grouping variables corresponding to samples in data
%
%
% OUTPUT:   p
%           p-value of sigificance between groups tested. Last p-value is
%           for all groups together
%
%           unique_pairs_name   2 x length(p)-1
%           categorical array listing testing pairs relating to the order
%           in p. 


[u_grs, ~, gr_num] = unique(groups_all);
if length(u_grs) > 2
    unique_pairs = [[];[]]; unique_pairs_name = unique_pairs;
    n = length(u_grs);
    for i = 1:n
        for j = i:n
            if i == j
                continue
            end
            unique_pairs_name = [unique_pairs_name [u_grs(i); u_grs(j)]];
            unique_pairs = [unique_pairs [i;j]];
        end
    end
else 
    unique_pairs = groups_all;
end

p = zeros(length(unique_pairs)+1, 1);
for gg = 1:length(unique_pairs)+1;
    if gg <= length(unique_pairs)
        groups = unique_pairs(:,gg);
        idg = find(gr_num == groups(1) | gr_num == groups(2));
        data = data_all(idg,:);
        X = gr_num(idg);
    else
        data = data_all; 
        X = gr_num;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Standardize input data
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
    
    %% calculate euclidiean dissimilarity matrix
    dissim = f_dis(dataszd,'euc');
    
    %% use for manova
    
    result = f_npManova(dissim, X,10000,0);
    p(gg) = result.p;
end
%% adjust pvalues using Benjimani-Hochberg method 
p = mafdr(p);