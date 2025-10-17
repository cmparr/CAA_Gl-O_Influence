function ax = PlotTS(data, limits, clr, LS, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot TS for given set of SA, pt values with same line color 
%
% INPUTS
%       data    struct (MxN)       
%           N profiles to plot for T/S. Must have CT and SA as fields to
%           plot
%
%       limits  double (2X2)
%           Limits of S and T to set axis. Pass empty array to free scale
%           plot. [Smin Smax; Tmin Tmax];
%
%       clr     double (3X1)
%           Color of profile to plot
%
%       LS      string 
%           Linestyle for plot. Pass empty for default (solid line)
%
%       VARARGIN
%           'thick'         Plot TS with thick line (Linewidth = 3)   
%           
%           'watermass'     Plot box outlining intermediate layer T/S
%
%
% OUTPUT
%       ax      NX1 double
%           Figure handle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Salinity contours %
    Si=20:0.1:36;
    thetai=-2:0.1:12;
    
    for j=1:length(thetai)
        for ii=1:length(Si)
            densi(j,ii)=gsw_rho(Si(ii),thetai(j),0);
        end
    end
    hold on
    v = linspace(1000, 1030, 31);
    [c,h]=contour(Si,thetai,densi,':k', 'HandleVisibility','off', 'LevelList',v);
    %[1025 1025.5 1026 1026.5 1027 1027.5 1028];
    clabel(c,h, v,'LabelSpacing',1000,'FontSize',10, 'color', [0.6 0.6 0.6])

%% loop through dataset
if isempty(LS)      % set linestyle if none selected
    LS = '-';
end

for m=1:length(data)
%% checks to make sure the data is in the correct formats
    if length(data(m).SA)~=length(data(m).pt)
        disp('Inputs are different sizes')
        continue
    elseif isempty(data(m).SA)||isempty(data(m).pt)
        continue
    end
    if any(strcmp(varargin, 'thick'))
    ax(m) = plot(movmean(data(m).SA,4, 'omitnan'), movmean(data(m).pt,4, 'omitnan'), 'Color',  clr, 'LineWidth',3, 'LineStyle',LS);
    else
    ax(m) = plot(movmean(data(m).SA,4, 'omitnan'), movmean(data(m).pt,4, 'omitnan'), 'Color',  clr, 'LineWidth',1, 'LineStyle',LS);
    end
    hold on
end

if ~isempty(limits)
xlim(limits(1,:)); 
ylim(limits(2,:));
end
ylabel('Potential Tempterature (C)', 'FontSize',12); xlabel('Absolute Salinity (g/kg)', 'FontSize',12);

%% plot Watermass ref points
if any(strcmp(varargin, 'watermass'))
xmin = 32.5 ; xmax = 33.8 ; ymin = -1.8 ; ymax = -0.5 ;
plot(34.53, 0.18, 'LineStyle','none', 'Marker','s', 'MarkerSize',7, 'MarkerFaceColor',[0.6350 0.0780 0.1840],...
        'MarkerEdgeColor', 'k','HandleVisibility','off'); % Atlantic Warm
plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin], 'LineStyle','--','color', [0.6350 0.0780 0.1840],...
 'HandleVisibility','off'); % Variable Arctic WM region
end
hold all
%%
end
