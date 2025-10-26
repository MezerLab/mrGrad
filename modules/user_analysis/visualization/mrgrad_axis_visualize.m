function [fig, ax] = mrgrad_axis_visualize(RG, axis_name, subject_idx, varargin)
%--------------------------------------------------------------------------
% INPUTS:
%--------------------------------------------------------------------------
%
%   RG          mrGrad output struct
%
%   subject_idx Logical index for subjects to visualize (can show 1 or more
%               subjects)
%
%   axis_name   axis for visualization (e.g., 'axis1')
%
%--------------------------------------------------------------------------
% Optional Name-Value arguments:
%--------------------------------------------------------------------------
%   'ax_view'   followed by an array of [azimuth, elevation] angles of the
%               camera's line of sight for the current axes. (see MATLAB's
%               'view' function);
%
%   'pc_lines'  followed by 0 or 1 (default), to display the axis line of
%               segmentation
%
%   'cmap'      m x 3 array of RGB triplets, with m at least as large as
%               the number of segments
%
%   'planes'    followed by 0 (default) or 1, to display the intersecting
%               planes of segmentation (work best for a single ROI)
%
%   'centroids' followed by 0 (default) or 1, to display the ROIs centroids
%
%   'segment'   followed by a segment number, to highlight the specific
%               single segment
%
%   'alpha'     specify the datapoints transparency
%
%   'markerSize' specify the datapoints size (default: 400)
%
%   'newFig'    followed by 1 (default) or 0. specify if to generate a new
%               figure (not relevant for multiple subjects)
%
%   'subplotInfo'   followed by field names for subplot titles (e.g.
%                   {'subject_names'})
%
% (C) Elior Drori, Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------
if ~isa(RG,'struct')
    error('input invalid');
end
rg_names = fieldnames(RG)';

if isscalar(axis_name)
    if ~ismember(axis_name,(1:3))
        error("Please provide a valid axis name: 'axis1' | 'axis2' | 'axis3' (or 1 | 2 | 3)");
    else
        axis_name = char("axis"+axis_name);
    end
end

% -------------------------------------------------------------------------
% Validate all ROIs have the chosen axis
axis_exist = cellfun(@(r) ismember(axis_name,RG.(r).axis_names), rg_names);
if any(~axis_exist)
    error('All ROIs should have the chosen axis');
end

% -------------------------------------------------------------------------
% Validate all ROIs share the same subjects
reqFields = {'subject_ids'};
for f = reqFields
    vals = cellfun(@(r) RG.(r).(f{1}), rg_names, 'un', 0);
    if ~all(cellfun(@(x) isequal(x, vals{1}), vals))
        error('All ROI structs must share identical %s values.', f{1});
    end
end

subjectIDs = RG.(rg_names{1}).subject_ids;
nSubjects_all = numel(subjectIDs);

if ~exist('subject_idx','var')
    subject_idx = [];
end
subject_idx = logical(subject_idx);

% analysis subject (up to 30)
if isempty(subject_idx)
    subject_idx = true(nSubjects_all, 1);
elseif ~isequal(numel(subjectIDs), numel(subject_idx))
    error('length of subject index should match the number of subjects.');
end

% Allow maximum 30 subjects
maxPlots = 30;
if nnz(subject_idx) > maxPlots
    trueIdx = find(subject_idx);
    subject_idx(trueIdx(31:end)) = false;
    warning('number of subjects is large. Plotting only first %d subject ids.',maxPlots);
end

subjectIDs(subject_idx);
nSubjects  = nnz(subject_idx);
subject_idx = find(subject_idx);
% -------------------------------------------------------------------------


[~, show_pc_lines, varargin] = argParse(varargin, 'pc_lines');
if isempty(show_pc_lines); show_pc_lines = true; end

[~, show_planes, varargin] = argParse(varargin, 'planes');
if isempty(show_planes); show_planes = false; end

[~, show_centroids, varargin] = argParse(varargin, 'centroids');
if isempty(show_centroids); show_centroids = false; end

[~, ax_view, varargin] = argParse(varargin, 'ax_view');
if isempty(ax_view); ax_view = [-45, 45]; end

[~, cmap, varargin] = argParse(varargin, 'cmap');
if isempty(cmap); cmap = 'autumn'; end

[found, Alpha, varargin] = argParse(varargin, 'alpha');
if ~found; Alpha = 0.2; end

[found, markerSize, varargin] = argParse(varargin, 'markerSize');
if ~found; markerSize = 400; end


fig = figure(WindowState="maximized");
nrows = floor(sqrt(nSubjects)); % subplot #rows
ncols = ceil(nSubjects/nrows); % subplot #cols

tiledlayout(nrows, ncols, 'TileSpacing', 'none');

ax = cell(nSubjects,1);
for ii = 1:nSubjects
    subIdx = subject_idx(ii);
    
    ax{ii} = nexttile;
    hold on
    grid on
    axis equal
    view(ax_view);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(subjectIDs{subIdx});


    for rr = 1:length(rg_names)
        rg = RG.(rg_names{rr});
        sub_rg = rg.individual_data{subIdx}(ismember(rg.axis_names,axis_name));
        nSegments = length(sub_rg.segment_inds);

       %--------------------
       % SCATTER DATA POINTS
       %--------------------
        h = cellfun(@(x) scatter3(x(:,1), x(:,2), x(:,3), 's'),...
            sub_rg.segment_inds);

        C = getColormap(cmap, nSegments);
        arrayfun(@(i) set(h(i), 'MarkerEdgeColor', 'none'), 1:nSegments);
        arrayfun(@(i) set(h(i), 'MarkerFaceColor', C(i,:)), 1:nSegments);
        arrayfun(@(i) set(h(i), 'MarkerFaceAlpha', Alpha), 1:nSegments);
        arrayfun(@(i) set(h(i), 'SizeData', markerSize), 1:nSegments);

        %------------------
        % EXTENDED PC LINES
        %------------------
        if show_pc_lines
            f = 0.2; % factor of extension
            
            l = sub_rg.analysisPC_line(:,[1,size(sub_rg.analysisPC_line,2)]);
            x1=l(1,1);
            y1=l(2,1);
            z1=l(3,1);
            x2=l(1,2);
            y2=l(2,2);
            z2=l(3,2);

            x3 = x2 + f * (x2 - x1);
            y3 = y2 + f * (y2 - y1);
            z3 = z2 + f * (z2 - z1);

            x0 = x1 + f * (x1 - x2);
            y0 = y1 + f * (y1 - y2);
            z0 = z1 + f * (z1 - z2);

            lnew = [x0 x3; y0 y3; z0 z3];
            h2 = line(lnew(1,:),lnew(2,:),lnew(3,:));
            h2.Color = [0 0 0];
            h2.LineWidth = 1.5;
        end

        %------------------
        % CENTROID
        %------------------
        if show_centroids
            C = sub_rg.data_centroid;
            pp = plot3(C(1),C(2),C(3),'.k');
            pp.MarkerSize = 16;
        end

        %------------------
        % PLANES
        %------------------
        if show_planes
            % save axes lims
            x_lim = xlim;
            y_lim = ylim;
            z_lim = zlim;
        
            h4 = cellfun(@(gz) surf(sub_rg.planes.gx, sub_rg.planes.gy, gz,...
                'FaceAlpha',0.15, 'EdgeColor', 'none'), sub_rg.planes.gz);

            % restore lims
            xlim(x_lim);
            ylim(y_lim);
            zlim(z_lim);
        end

    end
end
end

function cmap = getColormap(userCmap, n)
    % userCmap can be a string like 'parula' or a function handle like @summer
    if isa(userCmap, 'function_handle')
        cmap = userCmap(n);
    elseif ischar(userCmap) || isstring(userCmap)
        cmap = feval(userCmap, n);
    else
        error('userCmap must be a string (e.g., ''parula'') or a function handle (e.g., @summer).');
    end
end
