function conf = PlotSection(data, lims_in, vars, units, yval, maxy, xval, label, varargin)
%% 
% Plot section plots for given variables    
%
% conf = PlotSection(data, lims_in, vars, units, yval, maxy, xval, label, varargin)
%  
% INPUT: 
%       data    struct 
%               structure containing data to create one section
%
%       lims_in    MX3 double 
%               limits on each variable, in form upper, lower, interval
%               
%       vars    Mx1 cell 
%               variables to plot 
%
%       units   Mx1 cell 
%               corresponding units for plotted variables
%
%       yval    str 
%               which depth variable (depth or pressure) to plot with
%               respect to
%
%       xval    double 
%               Upper xlimit value. Use empty for dynamic plotting
%
%       maxy    double
%               max depth to set axes to
%
%       label    str
%               
%   VARARGIN
%       'log_chla'      sets log axis for chla
%       'log_turb'      sets log axis for turb
%       'log_n2'        sets log axis for n2
%       'obsv'          add triangles on surface to mark stations
%       'bottle'        add red triangles on surface to mark bottle stations
%       'isopycnals'    adds isopycnals to all plots
%       'lines'         plot with contour lines
%
%%
if isempty(data)
    return
end
%
log_vals = {};
if any(strcmp(varargin, 'log_chla'))
    log_vals = {log_vals; 'chla'};
elseif any(strcmp(varargin, 'log_turb'))
    log_vals = {log_vals; 'turb'};
elseif any(strcmp(varargin, 'log_n2'))
    log_vals = {log_vals; 'N2c'};
end

maxd=zeros(length(data),1);
len = [];

contour_lims = struct();
    
%% Set the x axis fill values
[distances, sorted_idcs] = find_distance(data); 
%data = data(sorted_idcs);

%% Set the y-axis fill values (depth axis)
%Determining max depth of all stations, and using as default depth
for a = 1:length(data)
len=[len length(data(a).(yval))];
maxd(a) = max(data(a).(yval));
end
[maxlen, maxpind] =max(len);
%ys=(1:maxlen)';
ys = data(maxpind).(yval);
ynans = find(isnan(ys)|isinf(ys));
if any(ynans > 1 & ynans < len(maxpind))
    indmid = ynans(find(ynans > 1 & ynans < len(maxpind)));
    ys(1) = 0.5;
    for ii = 1:length(indmid)
        ind = indmid(ii);
    ys(ind) = mean([ys(ind+1), ys(ind-1)]);
    end
end
for m = 1:length(vars)
    lims = lims_in(m,:);
%% Fill variable field
zs=[];
zs=NaN*ones(maxlen,length(data));
for i=1:length(data)
    xs = distances;
    % set station names
    stations = [data.station];
    %{
    % treat N-1 sized N2 data
    if strcmp(vars{m}, 'N2')
        arr = data(i).(vars{m});
        arr(end+1) = arr(end);
        data(i).(vars{m}) = arr;
        clear arr
    end
    %}
    %linear interp to fill in nans
    ind = find(isnan(data(i).(vars{m})));
    [~, idm] = max(data(i).(yval));

    %data(i).(vars{m})(ind(1)+1:ind(2)-1) = fillmissing(data(i).(vars{m})(ind(1)+1:ind(2)-1), 'linear', 'SamplePoints', ys(ind(1)+1:ind(2)-1));
    %data(i).(vars{m}) = fillmissing(data(i).(vars{m})(1:end), 'linear', 'SamplePoints', ys(2:length(data(i).(vars{m}))));
    if any(contains(log_vals, vars(m)))
         if strcmp('N2',vars(m))
            zs(1:length(data(i).(yval))-1,i) = log10(abs(data(i).(vars{m})));
         else
            zs(1:length(data(i).(yval)),i) = log10(abs(data(i).(vars{m})));
         end
    else
        if strcmp('N2', vars(m))   % N2 is one dimension shorter than other arrays
            zs(1:length(data(i).(yval))-1 ...
                ,i) = data(i).(vars{m});

        elseif strcmp('turb', vars(m))    % deal with negative values of turbidity (maybe shouldn't see this)
            data(i).(vars{m})(data(i).(vars{m}) < 0) = 0;
            zs(1:length(data(i).(yval)),i) = data(i).(vars{m});
        else
        
            zs(1:length(data(i).(yval)),i) = data(i).(vars{m});
        end
    end
    % fill missing values
    if isempty(ind)
    zs(2:idm, i) = fillmissing(zs(2:idm, i), 'linear', 'SamplePoints', ys(2:idm));
    else
    zs(ind(1)+1:idm, i) = fillmissing(zs(ind(1)+1:idm, i), 'linear', 'SamplePoints', ys(ind(1)+1:idm));
    end
end

% remove empty vars
%
indnan = find(sum(isnan(zs),1)==length(zs));
zs(:,indnan)     = [];
xs(indnan,:)     = [];
stations(indnan) = [];
%
if size(zs, 2)<3
    continue
end

% Contour plot
if length(vars)>1
nexttile;
end
[nX,nD,nT]=fill_tri2(xs,ys,zs);

% set contour fill values, 0.1 steps b/w min and max vals
if isempty(lims)
    if any(contains(log_vals,vars(m)))
    contour_lims(m).f = log10(abs(min(vertcat(data.(vars{m}))))):lims(3):log10(abs(max(vertcat(data.(vars{m})))));
    else
        contour_lims(m).f = min(vertcat(data.(vars{m}))):lims(3):max(vertcat(data.(vars{m})));
    end
else
    contour_lims(m).f = lims(1):lims(3):lims(2);
end
if any(strcmp(varargin, 'lines'))
    conf =contourf(nX,nD, real(nT), contour_lims(m).f, 'LineWidth', 1);
else
conf =contourf(nX,nD, real(nT), contour_lims(m).f, 'lineStyle', 'none');
end
hold on
axis ij

%% Add isopycnals
if any(strcmp(varargin, 'isopycnals'))
    % Fill density field
    [maxdval, maxdind] = max(maxd);
    rho_val = 'rho';
    rhofield = NaN*ones(length(data(maxdind).(rho_val)),length(data));
    for ii=1:length(data)
        rhofield(1:length(data(ii).(rho_val)), ii) = data(ii).(rho_val);
    end
    [X, Y, Z] = fill_tri2(distances, data(maxdind).depth, rhofield);
    levels = 1020:1:1029;
    [con,hcon] = contour(X, Y, Z, levels,'Color', 'Black', 'Linestyle', '--'); hold on;
    if isempty(maxy)
    clabel(con,hcon,'Color','black', 'Fontsize', 6, 'LabelSpacing', 550)
    else
    clabel(con,hcon,'Color','black', 'Fontsize', 8, 'LabelSpacing', 300)
    end
end

%% Mark CTD profile locations
if any(strcmp(varargin, 'obsv'))

plot([xs(:)],0, '^k', 'markerfacecolor', 'k', 'markersize', 5); hold on;
%text([xs(:)], ones(size(xs))*-1, [data.station], 'fontsize', 8, 'HorizontalAlignment','center');
%plot(0,0, 'sk', 'markerfacecolor', 'k', 'markersize', 3);
end
if any(strcmp(varargin, 'bottle'))
b_idcs = [data.bottlecast] == 'Y';
b_stations = xs(b_idcs);

% 1 for full column, 2 for upper 50 m
plot(b_stations, zeros(size(b_stations)), 'Marker', '^', 'MarkerFaceColor','r', 'MarkerEdgeColor','black', 'LineStyle','none','markersize',7); hold on;

end
%text([xs(:)],[xs(:)]*0+[-10], textstations, 'horizontalalign', 'center');
%{
if m == 1
text([xs(1)],[xs(1)]*0+(-20), stations(1), 'horizontalalign', 'center', 'Fontsize', 8);
%text([xs],[xs]*0+[-20], textstations, 'horizontalalign', 'center', 'Fontsize', 8);
text([xs(end)],[xs(end)]*0+(-20), stations(end), 'horizontalalign', 'center', 'Fontsize', 8);
end
end
%
%

%% Fill in lines for the points of observation
    for i=1:length(section)
        [mxd, ind] = max(section(i).depth);
        maxd(i) = mxd; idm(i) = ind;
    end
    [maxdval, maxdind] = max(maxd);
    obsvfill = ones([length(dists) round(maxdval/5)])*nan;
    for j = 1:length(dists)
        fill1 = 1: 50: maxd(j);
        obsvfill(j, 1:length(fill1)) = fill1;
    end
   %
        nexttile(m);
        plot(dists, obsvfill(:,:), 'Marker', 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor','black', ...
           'LineStyle','none'); hold all;
   %}
% set x limit
if ~isempty(xval)
    xlim([0 xval])
else
xlim([0 xs(end)]); %% add xlim so plots start at 0km from ref location
end
% set y limit
if isempty(maxy)
    ylim([0 max(maxd)+5]);
else
    ylim([0 maxy]);
end

colormap(parula)
%{
if m ~= length(vars) | m ~= length(vars) -1 
    set(gca, 'XTickLabel', '');
elseif m == length(vars) | m == length(vars) - 1
    xlabel('Distance from Terminus (km)');
end
%}
if mod(m,2) == 0
    set(gca, 'YTickLabel', '');
else 
    if strcmp(yval, 'depth')
    ylabel('Depth (m)');
    elseif strcmp(yval, 'press')
    ylabel('Pressure (dbar)');
end
end

%
%set(gca, 'XTickLabel', '');
if any(strcmp(varargin, 'cbar'))
cbar = colorbar; cbar.Label.String = strcat(label(m), ' ', units(m));
end
clim([min(contour_lims(m).f) max(contour_lims(m).f)]);
%}

clear indnan zs 
end
end

