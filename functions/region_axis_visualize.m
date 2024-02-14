function ff=region_axis_visualize(RG,groupRow,ROIsCols,sub_ind,ex_pc,varargin)
%--------------------------------------------------------------------------
% INPUTS:
%--------------------------------------------------------------------------
%
%   RG          cell of regiongram structs
%   
%   groupRow         Group index (RG row index)
%               e.g. 2
%
%   ROIsCols    ROIs index (RG column indices)
%               e.g. 1:4
%
%   sub_ind     subject indices (1 or more)
%               e.g. 1:length(RG{Row,1}.individual_data);
%
%   ex_pc       axis of segmentation
%               e.g. 1
%
%--------------------------------------------------------------------------
% Optional Name-Value arguments:
%--------------------------------------------------------------------------
%   'ax_view'   followed by an array of [azimuth, elevation] angles of the
%               camera's line of sight for the current axes. (see MATLAB's
%               'view' function);
%   'fig_ax'    0 or 1 to specify figure axis on / axis off
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
% (C) Elior Drori, Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

RG = RG(groupRow,:);


[found, show_pc_lines, varargin] = argParse(varargin, 'pc_lines');
if ~found; show_pc_lines = 1; end

[found, show_planes, varargin] = argParse(varargin, 'planes');
if ~found; show_planes = 0; end

[found, show_centroids, varargin] = argParse(varargin, 'centroids');
if ~found; show_centroids = 0; end

[found, fig_ax, varargin] = argParse(varargin, 'fig_ax');
if ~found; fig_ax = 1; end

[~, ax_view, varargin] = argParse(varargin, 'ax_view');

[~, cmap, varargin] = argParse(varargin, 'cmap');

[~, segment, varargin] = argParse(varargin, 'segment');

[found, Alpha, varargin] = argParse(varargin, 'alpha');
if ~found; Alpha = 0.2; end

[found, markerSize, varargin] = argParse(varargin, 'markerSize');
if ~found; markerSize = 400; end

[found, newFig, varargin] = argParse(varargin, 'newFig');
if ~found; newFig = 1; end


if fig_ax
    fig_color = 'w';
else
    fig_color = 'k';
end    

n = length(sub_ind);
r = floor(sqrt(n)); % subplot #rows
c = ceil(n/r); % subplot #cols

if newFig
    ff=figure('units','normalized','outerposition',[0 0 1 1]);
    ff.Color = fig_color;
end
a=1;
% a=15; warning('downsampling for faster visualization');
for ii=1:n % choose number of subjects
    
    if n > 1
        subplot(r,c,ii)
    end
    for jj=ROIsCols
    rg = RG{jj};
    ex = rg.individual_data{sub_ind(ii)}(ex_pc);

    bcmap = [178,24,43; 214,96,77; 244,165,130; 253,219,199; 209,229,240;...
         146,197,222; 67,147,195; 33,102,172]/255;
    bcmap = repmat(bcmap,15,1);

    if ~isempty(cmap) % if user provide colormap, use it unless it is too short
        maxN_segments = max(cellfun(@(x) x.N_segments(ex_pc),RG(groupRow,ROIsCols)));
        if size(cmap,1) < maxN_segments
            warning('Ignoring user provided colormap. Not enough RGB triplets provided.')
        else
            bcmap = cmap;
        end
    end
    
    lcmap = [0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5];
    lcmap(ex.analysisPC,:) = [1 0 0];
    
    % To notate specific segments
    scmap = [214,96,77;...
           50,50,50]/255;

    h = cell(ex.N_segments,1);
    h2= cell(3,1);
%------------------
% DATA POINTS
%------------------
for j=1:ex.N_segments % DATA POINTS
    inpoints = ex.segment_inds{j};
    h{j}=scatter3(inpoints(1:a:end,1),inpoints(1:a:end,2),inpoints(1:a:end,3),'s');
    if isempty(segment)
        h{j}.MarkerFaceColor = bcmap(j,:);
    else
        if ismember(j,segment)
            h{j}.MarkerFaceColor = scmap(1,:);
        else
            h{j}.MarkerFaceColor = scmap(2,:);
        end
    end
    h{j}.MarkerEdgeColor = 'none';
    h{j}.MarkerFaceAlpha = Alpha;
    h{j}.MarkerEdgeAlpha = Alpha;
    h{j}.SizeData = markerSize;
    hold on;
end

%------------------
% PC LINES
%------------------
%     for j=1:3%ex_pc %
%         l=rg.individual_data{ii}(j).analysisPC_line;
%         h2{j}=line(l(1,:),l(2,:),l(3,:));
%         h2{j}.Color = lcmap(j,:);
%         h2{j}.LineWidth = 4;
%         hold on;
%     end

%------------------
% EXTENDED PC LINES
%------------------
if show_pc_lines
    f = 0.2; % factor of extension
    for j=1
        l = ex(j).analysisPC_line(:,[1,size(ex(j).analysisPC_line,2)]);
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
        h2{j}=line(lnew(1,:),lnew(2,:),lnew(3,:));
        h2{j}.Color = [0 0 0];%[254,153,41]/255;
        h2{j}.LineWidth = 1.5;
        hold on;
    end
end

%------------------
% CENTROID
%------------------
if show_centroids
    C = rg.individual_data{ii}(j).data_centroid;
    pp=plot3(C(1),C(2),C(3),'.k');
    pp.MarkerSize = 16;
end
%------------------
% PC NORMAL VECTORS
%------------------
%     for j=1:3
%         c=RGG.individual_data{ii}(j).PC_coeff(:,j)' + RGG.individual_data{ii}(j).data_centroid;
%         p2=plot3(c(1),c(2),c(3),'*r');
%         p2.MarkerSize = 12;
%         hold on;
%     end
%------------------

axis equal
if ~fig_ax
    axis off
end

if ~isempty(ax_view)
    view(ax_view);
end


    end
    %------------------
    % PLANES
    %------------------
    if show_planes
        % save axes lims
        x_lim = xlim;
        y_lim = ylim;
        z_lim = zlim;
    
        for jj=ROIsCols
            rg = RG{jj};
            ex = rg.individual_data{sub_ind(ii)}(ex_pc);

            % PLANES

                for j=1:rg.N_segments+1
                    hs=surf(ex.planes.gx, ex.planes.gy, ex.planes.gz{j},'FaceAlpha',0.15, 'EdgeColor', 'none');
                    hold on;
                end
        end
        % restore lims
        xlim(x_lim);
        ylim(y_lim);
        zlim(z_lim);
    end
    
    if ~newFig
        ff = gca;
    end
    
end
