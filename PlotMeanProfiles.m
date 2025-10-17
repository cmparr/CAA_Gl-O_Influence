%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the open JS profiles for the supplementary figure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

%% set some figure properties
set(0, 'defaultAxesFontSize', 13);

% load data
load('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\Combined_CTD\combined_nuts_ctd_all_noflag.mat');

% 
clrs = cmocean('haline', 4)';
clrs(:,4) = brighten(clrs(:,4), -0.4);
yrs = 2019:1:2022;
cids = {'V1907';'V2007';'V2108'; 'V2208'};
% all the files
%
variables = {'SA', 'CT', 'CE', 'ODO_molL', 'chla', 'turb'};
labels = {'Absolute Salinity (g kg^{-1})'; 'Conservative Temperature (^{o}C)'; 'CR (N m^{-2})'; 'Dissolved Oxygen (mol L^{-1})'; 'Chorophyll (RFU)'; 'Turbidity (NTU)'};
titles = {'Salinity'; 'Conservative Temperature'; 'CR'; ''; 'Chorophyll-\alpha'; 'Turbidity'};
lims = [[27 35]; [-2 4]; [0 3500]; [200 450]; [0 30]; [0 6]];

%% ok now make some box charts
ax = figure; ax.Position = [0 0 1500 800];
tiledlayout(1,6,'TileSpacing','tight', 'Padding','compact');
for bkg = 2 %1:2
for yy = 1:4
    yr = yrs(yy);
    [mean_prof ctdf] = GetOutsideProf(cids{yy});
for xx = 1:length(variables)
    
    t = nexttile(xx);
    if bkg == 1
        for n = 1:length(ctdf)
            plot(ctdf(n).(variables{xx}), ctdf(n).depth, 'Color', brighten(clrs(:,yy)', 0.75), 'Linewidth', 0.95, 'HandleVisibility','off'); hold on;
        end
    elseif bkg ==2
        plot(mean_prof.(variables{xx}), mean_prof.depth, 'Color', clrs(:,yy), 'Linewidth', 2); hold on;
    end
    xlim(lims(xx,:)); 
    xlabel(labels{xx}); 
end
end
end
%
for nn = 1:6
    nexttile(nn);
    grid('on'); grid('minor');
    yline(30, 'LineStyle','--');
    ylabel('Depth (m)'); ylim([0 350]); 
    axis ij;
end

nexttile(6);
l =legend(num2str([2019:2022]'), 'Fontsize', 15, 'Location', 'southeast');
%
cd('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_plots\Plots-Claire\Mean_Profiles');
saveas(ax,['MeanProfiles_Year_upper350m.png']);

