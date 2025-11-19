function [FigureHandles] = mrgrad_show_gradients(RG,varargin)
%--------------------------------------------------------------------------
% mrgrad_show_gradients - Visualize mrGrad gradient trajectories
%
% INPUT:
%   RG              mrGrad output struct (one or more ROIs)
%
% Optional Name-Value Arguments:
%   'axis_names'            Axes to plot (figure columns; default: all)
%   'parameter_names'       MRI parameters to plot (figure rows; default: all)
%   'group_labels'          Group label array, one per subject
%   'show_individuals'      Binary array (1×Ngroups), show subject trajectories
%   'show_gradient_numbers' Logical, display gradient indices on curves
%   'error_name'            'STD' (default) or 'SEM'
%   'equal_ylims'           Logical, equal Y range across gradients for
%                           per parameter.
%
% OUTPUT:
%   FigureHandles.figure
%   FigureHandles.group_plot_handles
%   FigureHandles.individuals_plot_handles
%   FigureHandles.legend
%
% (C) Elior Drori, Mezer Lab, HUJI, 2021
%--------------------------------------------------------------------------

% Input Validation
if ~isa(RG,'struct')
    error('input invalid');
end
rg_names = fieldnames(RG)';
nrg = length(rg_names);

% -------------------------------------------------------------------------
% Validate all ROIs share the same structure
reqFields = {'axis_names','n_segments','parameter_names','subject_ids'};
for f = reqFields
    vals = cellfun(@(r) RG.(r).(f{1}), rg_names, 'un', 0);
    if ~all(cellfun(@(x) isequal(x, vals{1}), vals))
        error('All ROI structs must share identical %s values.', f{1});
    end
end

axisList   = RG.(rg_names{1}).axis_names;
paramList  = RG.(rg_names{1}).parameter_names;
subjectIDs = RG.(rg_names{1}).subject_ids;
nSubjects  = numel(subjectIDs);
% -------------------------------------------------------------------------
% Parse Optional Arguments
[~, parameter_names, varargin] = argParse(varargin, 'parameter_names');
if isempty(parameter_names); parameter_names = paramList; end
parameter_names = cellstr(parameter_names);

[~, axis_names, varargin] = argParse(varargin, 'axis_names');
if isempty(axis_names); axis_names = axisList; end

[~, group_labels, varargin] = argParse(varargin, 'group_labels');
if isempty(group_labels); group_labels = ones(nSubjects,1); end

group_labels = categorical(group_labels);
uniqueGroups = categories(group_labels);
nGroups = numel(uniqueGroups);
groupSize = countcats(group_labels);

[~, error_name, varargin] = argParse(varargin, 'error_name');
if isempty(error_name) || ~ismember(lower(error_name),{'std','sem'})
    error_name = 'std';
else
    error_name = lower(error_name);
end

[~, show_gradient_numbers, varargin] = argParse(varargin, 'show_gradient_numbers');
if isempty(show_gradient_numbers); show_gradient_numbers = false; end

[~, show_individuals, varargin] = argParse(varargin, 'show_individuals');
if isempty(show_individuals); show_individuals = zeros(1, nGroups); end
if isscalar(show_individuals)
    show_individuals = repmat(show_individuals, 1, nGroups);
end

[~, equal_ylims, varargin] = argParse(varargin, 'equal_ylims');
if isempty(equal_ylims); equal_ylims = true; end

%--------------------------------------------------------------------------
% Prepare Data
nAxes = length(axis_names);
nParam = length(parameter_names);

gradientData = structfun(@(x) x.Results,RG);
dataFields = fieldnames(gradientData);

% keep only specified axes and parameters
validFields = contains(dataFields, axis_names) & ...
    contains(dataFields,parameter_names);
gradientData = rmfield(gradientData,dataFields(~validFields));

useBoundedLine = exist('boundedline','file');

roiColors = lines(nrg);
groupColors = lines(nGroups);

markerStyles = {'-.','-o', '-+', '-*', '-x', '-s', '-d', '-^', '-v', '->', '-<', '-p', '-h'};


% Plot Initialization
fig = figure(WindowState="maximized");
pause(.2);
plotHandles = cell(nParam, nAxes);
individualPlotHandles = cell(nParam, nAxes);

if nAxes  > 1
    t = tiledlayout(nParam, nAxes, 'TileSpacing', 'compact');
else
    t = tiledlayout(nAxes, nParam, 'TileSpacing', 'compact');    
end
for p = 1:nParam
    param = parameter_names{p};
    for a = 1:nAxes
        axname = axis_names{a};
        fieldName = axname +"_"+ param;
        ax=nexttile;
        ax.Tag = fieldName;
        xlabel(axname);
        ylabel(param);
        axis square
        grid on
        hold on

        if nParam>1
            subplotIdx = (p-1)*nAxes + a;

            % Bold y label (parameter name) in first tile of each row
            if mod(subplotIdx-1, nAxes) == 0
                ax.YLabel.FontWeight = 'bold';
            end
            % Bold x label (axis name) in bottom row
            if p==nParam
                ax.XLabel.FontWeight = 'bold';
            end
        end

        plotHandles{p,a} = [];
        individualPlotHandles{p,a} = [];
        for rr = 1:nrg
            xVals = double((gradientData(rr).(fieldName).x));
            YVals = double(gradientData(rr).(fieldName).y{:,:});

            for g = 1:nGroups
                groupMask = group_labels == uniqueGroups{g};
                y = YVals(groupMask,:);
                mu = mean(y,1,"omitnan");
                sigma = std(y,[],1,"omitnan");
                sem = sigma / sqrt(size(y,1));

                switch error_name
                    case 'std'
                        err = double(sigma);
                    case 'sem'
                        err = double(sem);
                end


                if nrg > 1
                    cmap = roiColors(rr,:);
                    markerShape = markerStyles{g};
                else
                    cmap = groupColors(g,:);
                    markerShape = '-.';
                end

                % Actual group size for this parameter (omit nans)
                groupSize_actual = nnz(~any(isnan(y),2));

                if groupSize_actual==1 % only 1 subject in group
                    h = plot(xVals,y, markerShape, Color=cmap);
                else

                    if useBoundedLine
                        h = boundedline(xVals,mu,err,markerShape,'alpha','cmap',cmap,'transparency',0.3);
                    else
                        h = errorbar(xVals,mu,err,markerShape,Color=cmap);
                    end
                end
                plotHandles{p,a} = [plotHandles{p,a},h];

                % plot single subjects
                if show_individuals(g) && size(y,1) > 1
                     h2 = plot(xVals,y,'-');
                     arrayfun(@(x) set(x,'Color',cmap),h2);
                     
                     if show_gradient_numbers
                         arrayfun(@(i) ...
                             text(xVals(end),y(i,end), num2str(i),Color=cmap), ...
                             1:size(y,1));

                         [~, maxIdx] = max(y,[],2);
                         arrayfun(@(i) ...
                             text(xVals(maxIdx(i)),y(i,maxIdx(i)), num2str(i),Color=cmap), ...
                             1:size(y,1));
                     end
                individualPlotHandles{p,a} = [individualPlotHandles{p,a}, h2];
                end
            end
        end
        xlim(minmax(xVals))
        xticks(xVals);
    end

    if equal_ylims
        all_fig_axes = findall(fig,'type','axes');
        param_fig_axes = all_fig_axes(1:nAxes);
        y_lims = arrayfun(@(x) get(x,'YLim'), param_fig_axes,'un',0);
        y_lims = cat(1,y_lims{:});
        y_lims = [min(y_lims(:,1)), max(y_lims(:,2))];
    
        arrayfun(@(x) set(x,'YLim',y_lims), param_fig_axes);
    end
end


% Legend (All plots)
Axes = findall(t.Children, 'type', 'Axes');
Axes = flip(Axes);
for jj = 1:numel(Axes)
    ax = Axes(jj);
    y = gradientData(1).(ax.Tag).y{:,:};
    group_names = string(uniqueGroups(:));
    groupSize_actual = arrayfun(@(g) nnz(~any(isnan(y(group_labels==g,:)),2)), group_names);

    lgd_grps = compose("%s (n=%d)",group_names,groupSize_actual)';
    lgd_rois = strrep(string(rg_names(:)),'_',' ');

    lgd_entries = compose("%s: %s",lgd_rois, lgd_grps)';
%     if length(lgd_rois) > 1
%         lgd_entries = compose("%s: %s (n=%d)",lgd_rois',group_names,groupSize_actual')';
%     else
%         lgd_entries = compose("%s: %s (n=%d)",lgd_rois,group_names,groupSize_actual);
%     end
    lgd_entries = lgd_entries(:);
    legend(ax, plotHandles{jj},lgd_entries);
end



% % Legend (first plot)
% lgd_entries = strrep(string(rg_names),'_',' ');
% if nGroups>1
%     group_names = unique(group_labels)';
%     if length(lgd_entries) > 1
%         lgd_entries = compose("%s: %s (n=%d)",lgd_entries',string(group_names),groupSize')';
%     else
%         lgd_entries = compose("%s: %s (n=%d)",lgd_entries,string(group_names)',groupSize);
%     end
%     lgd_entries = lgd_entries(:);
% end
% legend(t.Children(end),plotHandles{1,1},lgd_entries);

FigureHandles = struct;
FigureHandles.figure = fig;
FigureHandles.tiledlayout = t;
FigureHandles.group_plot_handles = plotHandles;
FigureHandles.individuals_plot_handles = individualPlotHandles;
end