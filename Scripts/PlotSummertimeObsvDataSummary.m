%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot summary of data by site classification. Figure 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear

%
load('../fucntions/Colormaps-ModelStory.mat');
set(0, 'defaultAxesFontSize', 15);


%% combineds dataset
load('../Data/combined_ctd_metadata.mat');
load('../Data/combined_nuts_ctd_all.mat');
combined.mmdd= [month(combined.Date_yyyy_mm_dd_), day(combined.Date_yyyy_mm_dd_)];
combined = combined(combined.dist < 20.,:);
idx = find(~strcmp(combined.Transect, 'Sydkap Glacier - across') & [combined.turb] > 0 & ~contains([combined.Location], {'Talbot'; 'Jones Sound'; 'Harbour';'Terry';'OG';'to'}));
combined = combined(idx,:);
%combineds.FWsourcec = categorical(combineds.FWSource); combineds.Coastlinec = categorical(combineds.Coastline); % turn these categorizations to categoricals to histogram them
vars = {'Depth_m'; 'FWSource'; 'Coastline'; 'dist'; 'Year'; 'mmdd'};
edges = [[0 10 200];[0 0 0];[0 0 0]; [0 0.5 30]; [0 10 120]; [0 0 0]];
labels = {'Sampling Depth (m)'; 'Freshwater Source'; 'Coastline'; 'Station Distance from Coast (km)'; 'Year'; 'Month/day'};

%% CTD + nuts 
f = figure; f.Position = [40 40 1200 450];
tiledlayout(1,3, 'Padding','compact', 'TileSpacing','tight');
combineds = combined(combined.Depth_m < 30.,:); % surface sample summary (s)
combinedd = combined(combined.Depth_m > 30. & combined.Depth_m < 100.,:); % deeper sample summary (d)
counts = [height(combineds(combineds.Year==2019 & combineds.combinedchar=='Glacier-Open',:)),height(combineds(combineds.Year==2019 & combineds.combinedchar=='Glacier-Fjord',:)),height(combineds(combineds.Year==2019 & combineds.combinedchar=='River-Fjord',:)),height(combineds(combineds.Year==2019 & combineds.combinedchar=='River-Open',:));
        height(combineds(combineds.Year==2020 & combineds.combinedchar=='Glacier-Open',:)),height(combineds(combineds.Year==2020 & combineds.combinedchar=='Glacier-Fjord',:)),height(combineds(combineds.Year==2020 & combineds.combinedchar=='River-Fjord',:)),height(combineds(combineds.Year==2020 & combineds.combinedchar=='River-Open',:));
        height(combineds(combineds.Year==2021 & combineds.combinedchar=='Glacier-Open',:)),height(combineds(combineds.Year==2021 & combineds.combinedchar=='Glacier-Fjord',:)),height(combineds(combineds.Year==2021 & combineds.combinedchar=='River-Fjord',:)),height(combineds(combineds.Year==2021 & combineds.combinedchar=='River-Open',:));
        height(combineds(combineds.Year==2022 & combineds.combinedchar=='Glacier-Open',:)),height(combineds(combineds.Year==2022 & combineds.combinedchar=='Glacier-Fjord',:)),height(combineds(combineds.Year==2022 & combineds.combinedchar=='River-Fjord',:)),height(combineds(combineds.Year==2022 & combineds.combinedchar=='River-Open',:))];
countd = [height(combinedd(combinedd.Year==2019 & combinedd.combinedchar=='Glacier-Open',:)),height(combinedd(combinedd.Year==2019 & combinedd.combinedchar=='Glacier-Fjord',:)),height(combinedd(combinedd.Year==2019 & combinedd.combinedchar=='River-Fjord',:)),height(combinedd(combinedd.Year==2019 & combinedd.combinedchar=='River-Open',:));
        height(combinedd(combinedd.Year==2020 & combinedd.combinedchar=='Glacier-Open',:)),height(combinedd(combinedd.Year==2020 & combinedd.combinedchar=='Glacier-Fjord',:)),height(combinedd(combinedd.Year==2020 & combinedd.combinedchar=='River-Fjord',:)),height(combinedd(combinedd.Year==2020 & combinedd.combinedchar=='River-Open',:));
        height(combinedd(combinedd.Year==2021 & combinedd.combinedchar=='Glacier-Open',:)),height(combinedd(combinedd.Year==2021 & combinedd.combinedchar=='Glacier-Fjord',:)),height(combinedd(combinedd.Year==2021 & combinedd.combinedchar=='River-Fjord',:)),height(combinedd(combinedd.Year==2021 & combinedd.combinedchar=='River-Open',:));
        height(combinedd(combinedd.Year==2022 & combinedd.combinedchar=='Glacier-Open',:)),height(combinedd(combinedd.Year==2022 & combinedd.combinedchar=='Glacier-Fjord',:)),height(combinedd(combinedd.Year==2022 & combinedd.combinedchar=='River-Fjord',:)),height(combinedd(combinedd.Year==2022 & combinedd.combinedchar=='River-Open',:))];
nexttile;
b=bar(2019:2022, counts, 'stacked', 'FaceColor','flat');
b(1).FaceColor = budget(1,:); b(2).FaceColor = budget(2,:); b(3).FaceColor = budget(3,:); b(4).FaceColor = budget(5,:);
ylabel('Upper 30 m depth Bottle Sample Count'); 
xlim([2018.5 2022.5]); ylim([0 100]);
%title('<30 m Depth');
%t = title('$(a)$', 'Interpreter','latex','Units','normalized', 'HorizontalAlignment','left'); t.Position(1) = 0; 
%saveas(gca, [rootfol 'BottleSampleSummary-BelowNutricline-ClassifStacked.png']);

nexttile;
b=bar(2019:2022, countd, 'stacked', 'FaceColor','flat');
b(1).FaceColor = budget(1,:); b(2).FaceColor = budget(2,:); b(3).FaceColor = budget(3,:); b(4).FaceColor = budget(5,:);
l = legend({'Glacierized-Open';'Glacierized-Fjorded';'Riverine-Fjorded';'Riverine-Open'}, 'Location', 'northeast');
ylabel('30-100 m depth Bottle Sample Count'); 
xlim([2018.5 2022.5]); ylim([0 100]);
%title('30-100m Depth');
%t = title('$(b)$', 'Interpreter','latex','Units','normalized', 'HorizontalAlignment','left'); t.Position(1) = 0; 
%sgtitle('Sample Summary For Bottle Samples')
sgtitle(' ');
%saveas(gca, [rootfol 'BottleSampleSummary-BelowNutricline-ClassifStacked.png']);
%% CTD 
%f = figure; f.Position = [40 40 800 600];
idx = find(~strcmp(metad.transect, 'Sydkap Glacier - across') & ~contains([metad.transect], {'Talbot'; 'Jones Sound'; 'Harbour';'Terry';'OG';'to'}));
metad = metad(idx,:);
counts = [height(metad(metad.year==2019 & metad.combchar=='Glacier-Open',:)),height(metad(metad.year==2019 & metad.combchar=='Glacier-Fjord',:)),height(metad(metad.year==2019 & metad.combchar=='River-Fjord',:)),height(metad(metad.year==2019 & metad.combchar=='River-Open',:));
          height(metad(metad.year==2020 & metad.combchar=='Glacier-Open',:)),height(metad(metad.year==2020 & metad.combchar=='Glacier-Fjord',:)),height(metad(metad.year==2020 & metad.combchar=='River-Fjord',:)),height(metad(metad.year==2020 & metad.combchar=='River-Open',:));
          height(metad(metad.year==2021 & metad.combchar=='Glacier-Open',:)),height(metad(metad.year==2021 & metad.combchar=='Glacier-Fjord',:)),height(metad(metad.year==2021 & metad.combchar=='River-Fjord',:)),height(metad(metad.year==2021 & metad.combchar=='River-Open',:));
          height(metad(metad.year==2022 & metad.combchar=='Glacier-Open',:)),height(metad(metad.year==2022 & metad.combchar=='Glacier-Fjord',:)),height(metad(metad.year==2022 & metad.combchar=='River-Fjord',:)),height(metad(metad.year==2022 & metad.combchar=='River-Open',:))];
nexttile(3);
b=bar(2019:2022, counts, 'stacked', 'FaceColor','flat');
b(1).FaceColor = budget(1,:); b(2).FaceColor = budget(2,:); b(3).FaceColor = budget(3,:); b(4).FaceColor = budget(5,:);
ylabel('CTD Cast Count'); 
xlim([2018.5 2022.5]); %ylim([0 80]);
%title('CTD Sample Summary');
%l = legend({'Glacierized-Open';'Glacierized-Fjorded';'Riverine-Fjorded';'Riverine-Open'}, 'Location', 'northwest');
%t = title('$(c)$', 'Interpreter','latex','Units','normalized', 'HorizontalAlignment','left'); t.Position(1) = 0; 

saveas(gca, ['../figures/SampleSummary-combined-ClassifStacked.png']);
%%