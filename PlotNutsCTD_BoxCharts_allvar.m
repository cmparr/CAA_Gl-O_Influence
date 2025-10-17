%%%%%%%%%%%%%%%%%%%%%%%%%%
% Box and whisker plots for relevant variables and all years. Figure 3
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

%% set some figure properties
set(0, 'defaultAxesFontSize', 13);

% load data
load('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\Combined_CTD\combined_nuts_ctd_all_noflag_dz5.mat');

LRce = struct(); LRfwt = struct(); ll = 1;
years = 2019:1:2022;
% all the files
combined = combined((combined.dist < 20. & ~contains(string(combined.Transect), {'Jones Sound'; 'Harbour'; 'Shore'})), :);
combined = combined(((combined.Depth_m < 100.)& combined.dist < 30. & ~contains(string(combined.Transect), {'Jones Sound'; 'Harbour'; 'Shore'})), :);
combined.depth_bin = repmat(">100m", [height(combined) 1]); combined(combined.Depth_m <= 30.,:).depth_bin = repmat("Upper 30m", [height(combined(combined.Depth_m <= 30.,:)) 1]);  combined(combined.Depth_m > 30.&combined.Depth_m <= 100.,:).depth_bin = repmat("30-100m", [height(combined(combined.Depth_m > 30.& combined.Depth_m <= 100.,:)) 1]); 
combined.depth_bin = categorical(combined.depth_bin);
combined.CEd = combined.CE; combined(combined.Depth_m <= 30.,:).CEd = combined(combined.Depth_m <= 30.,:).CE30;  combined(combined.Depth_m > 30.&combined.Depth_m <= 100.,:).CEd = combined(combined.Depth_m > 30.& combined.Depth_m <= 100.,:).CE100; 

%
variables = {'sa', 'CT', 'turb', 'DO_Lmol', 'NO3_uM','PO4_uM', 'SiO4_uM', 'chla'};
labels = {'Absolute Salinity (g kg^{-1})'; 'Conservative Temperature (^{o}C)'; 'Log_{10}(Turbidity) (NTU)'; 'Dissolved O_2 (mol L^{-1})';
    'NO_3 (\mu M)'; 'PO_4 (\mu M)'; 'SiO4 (\mu M)';  'Chorophyll-\alpha (RFU)'};
titles = {'Salinity'; 'Conservative Temperature'; 'Turbidity';
    'Nitrate'; 'Phosphate'; 'Silicate'; 'DO'; 'Chorophyll-\alpha'};
lims = [[23 34.7]; [-2 5.3]; [0 1.5e3];  [250 460]; [0 12.5]; [0 1.25]; [0 21]; [0 53]];
ypos = [34 5 1e3 450 12 1.2 20 51];
%% ok now make some box charts
ax = figure; ax.Position = [0 0 1300 800];
tiledlayout(2,4,'TileSpacing','tight', 'Padding','compact');
for xx = 1:length(variables)
    % find n.o. observations in each grouping
    u19 = sum(combined.depth_bin == 'Upper 30m' & combined.Year == 2019 & ~isnan(combined.(variables{xx}))); u20 = sum(combined.depth_bin == 'Upper 30m' & combined.Year == 2020 & ~isnan(combined.(variables{xx}))); u21 = sum(combined.depth_bin == 'Upper 30m' & combined.Year == 2021 & ~isnan(combined.(variables{xx}))); u22 = sum(combined.depth_bin == 'Upper 30m' & combined.Year == 2022 & ~isnan(combined.(variables{xx})));
    b19 = sum(combined.depth_bin == '30-100m' & combined.Year == 2019 & ~isnan(combined.(variables{xx})));   b20 = sum(combined.depth_bin == '30-100m' & combined.Year == 2020 & ~isnan(combined.(variables{xx})));   b21 = sum(combined.depth_bin == '30-100m' & combined.Year == 2021 & ~isnan(combined.(variables{xx})));   b22 = sum(combined.depth_bin == '30-100m' & combined.Year == 2022 & ~isnan(combined.(variables{xx})));
    %   plot
    t = nexttile(xx);
    boxchart(combined.depth_bin, combined.(variables{xx}), 'GroupByColor', combined.Year, 'BoxEdgeColor', 'k', 'BoxFaceAlpha',0.9, 'MarkerStyle', '.', 'MarkerSize', 20); colororder(cmocean('haline',4)); 
    text([0.55 0.8 1.04 1.29 1.54 1.79 2.04 2.28] , repmat(ypos(xx),1,8), num2str([u19 u20 u21 u22 b19 b20 b21 b22]'), 'Fontsize', 12); 
    xline(1.5, ':','LineWidth', 1, 'HandleVisibility','off');
    if any(strcmp(variables{xx}, {'turb'}))
        t.YAxis.Scale ="log";
    end
    ylim(lims(xx,:)); 
    ylabel(labels{xx}); 
    %title(titles{xx});
    xaxis=get(gca,'XAxis');
    xaxis.Categories={'Upper 30m'; '30-100m'};
end
nexttile(2);
l =legend(num2str([2019:2022]'), 'Fontsize', 15, 'Position', [0.9012 0.58 0.0805 0.1029]);
% I usually pause here to move the legend to a better spot on the figure
rootfol = 'C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_plots\Plots-Claire\BoxCharts\';
saveas(ax,[rootfol 'Boxchart-newPCAVars-twodepths-notitle-logTu-CEd.png']);


