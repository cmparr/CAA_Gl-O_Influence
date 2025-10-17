%%%%%%%%%%%%%%%%%%%%%%%
% Plot mean CTD variables by each glacier site. Figure 8.
%%%%%%%%%%%%%%%%%%%%%%%

clear all
close

set(0, 'defaultAxesFontSize', 14);

%% what is going on with the upwelling depth
load('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\Combined_CTD\combined_nuts_ctd_all_noflag_dz5.mat');

%% check the glaciers are upwelling year-over-year
combined = combined((combined.dist < 20. & ~contains(string([combined.Transect]), {'Jones Sound'; 'Harbour'; 'Shore';'Point';'OG';'Fram';'Talbot';'Terry';'to'})),:); % sites within 20km of shore 
combined = combined(combined.Depth_m < 30.,:);
combined_g = combined((contains(string([combined.Location]), {'Belcher'; 'Sverdrup'; 'Sydkap'; 'South';'Jakeman'})& ...
    contains(string([combined.Transect]), 'out')),:); 
combined = convertvars(combined_g, ['Location'], 'categorical');
% group by site name
MeanTable = groupsummary(combined_g,'Location','mean',{'glmelt1'});

%% replicating/extending the Bhatia 21 figure
%f1 = figure(1); f1.Position = [50 50 1000 800]; 
%tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
f2 = figure(2); f2.Position = [50 50 1300 600];
tiledlayout(2,4,"TileSpacing","tight","Padding","tight");
hold on;
glaciers = {'Jakeman';'Sverdrup';'Sydkap';'Belcher'};
gl_short = {'JK';'SV';'SD';'BL'};
%colors = flipud(cmocean('matter',4));
col_yrs = cmocean('haline',4);
yrs = 2019:2022;
markers = {'square';'o';'^';'pentagram'};
%indicator = {'turb';'NO3_uM';'chla'};
%labels = {'Turbidity (NTU)';'Mean Above-nutricline [NO_3] (\muM)';'Chlorophyll-\alpha (\mu g L^{-1})'};
indicator = {'NO3_uM';'PO4_uM'; 'SiO4_uM'; 'ODO_c';'turb_c';'chla_c';'CE30';'dz'};
%labels = {'Mean Above-nutricline [NO_3] (\muM)';'Mean Above-nutricline [PO_4] (\muM)';'Mean Above-nutricline [SiO_4] (\muM)';'Turbidity (NTU)';'Chlorophyll-\alpha (\mu g L^{-1})'; 'Convective Resistance (J/m^3)';'Dissolved O_2 M (mol L^{-1})'};
labels = {'[NO_3] (\muM)';'[PO_4] (\muM)';'[SiO_4] (\muM)'; 'Diss. O_2 (mol L^{-1})';...
        'Turbidity (NTU)';'Chlorophyll-\alpha (RFU)';'CR (J/m^3)'; '\Delta z (m)'};        
lims = [[-0.5 13];[0 1.2];[-0.5 20];[300 425]; ...
          [-0.5 25];[-0.5 20];[0 1000];[-25 45]]; 

r2s = []; ps= []; mean_vals = zeros(4,4); mkr_sz = mean_vals; cols_gl = zeros(4,4,3);
for ii = 1:length(indicator)
        X = []; Y = [];
        indicator{ii}
    for g = 1:4
        gl = glaciers(g);
        for yy = yrs
        glacier = combined(contains(string([combined.Transect]), gl)& contains(string([combined.Transect]), 'out')& combined.Year == yy & combined.Depth_m < 30. ,:);
        % compile mean value matrix, rows are years columns are glaciers
            mean_vals(yy-2018,g) = mean(glacier.(indicator{ii}), 'omitnan');
            std_vals(yy-2018,g) = std(glacier.(indicator{ii}), 'omitnan');
            mkr_sz(yy-2018,g) = mean(glacier.glmelt1)*20;
        end
    end
    mkr_sz(isnan(mkr_sz)) = 0; mkr_sz(mkr_sz == 0) = 0.1;
    %
    %% values figure
    figure(2); nexttile(ii);
    for yyy = 1:4
        col = col_yrs(yyy,:); %
        scatter(1:4, mean_vals(yyy,:), 130, 'filled', 'MarkerFaceColor', col, 'marker', 'o', 'HandleVisibility','off','MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5); hold on;
        errorbar(1:4, mean_vals(yyy,:), std_vals(yyy,:), 'Linestyle', 'none', 'Color', col,'HandleVisibility','off'); hold on;
    end
    %scatter(reshape(mkr_sz/20,[],1), reshape(mean_vals, 1,[]), 130, reshape(cols_gl,[],3), 'filled','MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8, 'HandleVisibility','off'); hold on;
   % scatter(repmat(1:4, 1,4), reshape(mean_vals, 1,[]), reshape(mkr_sz,[],1), reshape(repmat(col_yrs,4,1),[],3), 'filled','MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8, 'HandleVisibility','off'); hold on;
end
 
% legend
for yy = yrs
    figure(2);
    nexttile(length(indicator));
    scatter(-10, yy/100, 250, 'filled', 'MarkerFaceColor', col_yrs(yy-2018,:), 'MarkerEdgeColor','k'); hold on;
end

figure(2);
for ii = 1:length(indicator)
    nexttile(ii);
    ax = gca;
    ax.XTick = 1:4; 
    ax.XTickLabel = gl_short;
    xlim([0.5 4.5]); %xlim([0 210]); 
    %xlabel('Mean Daily Glacier Melt (kg m^3)'); % grid(); 
    %xlabel('Glacier Terminus Depth (m)'); % grid(); 
    yline(0, 'color', 'k', 'HandleVisibility','off');
    grid();
    ylabel(labels{ii}); ylim(lims(ii,:));
end
nexttile(1);
%yline(0.1788, '--','HandleVisibility','off'); %legend('Mean riverine [NO3] (\muM)');
nexttile(length(indicator));
legend({'2019';'2020';'2021';'2022'});
%legend(glaciers);

