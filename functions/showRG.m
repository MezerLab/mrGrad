function h = showRG(RG,varargin)
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
%   RG          cell of mrGrad structs
%--------------------------------------------------------------------------
% Optional Name-Value arguments:
%--------------------------------------------------------------------------
%   'PC':       Followed by an array of axis numbers to show (default 1:3)
%
%   'ind':      A binary array in length (NROIs x Ngroups), specifying
%               whether to display individual subjects trajectories for each
%               mean function. e.g. [1 0 0 0]
%
%   'errorType': 'STD' (1) / 'SEM' (2) / 'none' (0)
%
%   'MeanFunction': @mean (default) or @median
%
%   'noBL'      0 or 1 to ommit group-average function (display only
%               individual trajectories) 
%
%   'markershapes' 0 or 1 to add marker shapes to mean plots
%                  (alternitavely - specify specific markershapes, e.g.
%                  {'-o','-s'}).
%
%   'legend'    0 or 1 for showing default legend, or a cell array
%               containing legend entries
%
%   'labels'    0 or 1 for showing axis labels
%
%   'gradnumber' 0 or 1 for indicating individual gradients numbering
%
%   'ylim'      specify ylim (1x2 array for all subplots, or (nPCs x 2) for
%               specific ylim for each subplot.
%
%   'yticks'    specify yticks (1 x nPCs cell array)
%
%   'yticklabels'  specify yticklabels
%
%   'ylbl'      specify ylabel
%
%   'xlbl'      specify xlabel
%
%   'cmap'      specify colors (one RGB row per function)
%
%   'square' 0 or 1(default) for square figure axis
%
%   'axis_spec' default is 'tight'
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

[found, ind, varargin] = argParse(varargin, 'ind');
if ~found; ind = zeros(1,numel(RG)); end
if numel(ind)==1
    ind = ind*ones(size(RG));
end

[found, PC, varargin] = argParse(varargin, 'PC');
if ~found; PC = 1:numel(RG{1}.Y); end

[found, use_err, varargin] = argParse(varargin, 'errorType');
if ~found; use_err = 1; end

[found, MeanFunction, varargin] = argParse(varargin, 'MeanFunction');
if ~found; MeanFunction = @mean; end

[found, markershapes, varargin] = argParse(varargin, 'markershapes');
if ~found; markershapes = false; end
if ischar(markershapes)
    markershapes = {markershapes};
end

[found, leg, varargin] = argParse(varargin, 'legend');
if ~found; leg = false; end
if iscell(leg)
    leg_ent = leg;
    leg = true;
end

[found, noBL, varargin] = argParse(varargin, 'noBL');
if ~found; noBL = false; end

[found, lbl_flag, varargin] = argParse(varargin, 'labels');
if ~found; lbl_flag = true; end

[found, gradnumber, varargin] = argParse(varargin, 'gradnumber');
if ~found; gradnumber = false; end

[found, y_lbl, varargin] = argParse(varargin, 'ylbl');
if ~found
    clearvars y_lbl;
end

[found, x_lbl, varargin] = argParse(varargin, 'xlbl');
if ~found
    clearvars x_lbl;
end


[found, y_lim, varargin] = argParse(varargin, 'ylim');
if ~found || isempty(y_lim)
    clearvars y_lim;
end

[found, y_ticks, varargin] = argParse(varargin, 'yticks');
if ~found || isempty(y_ticks)
    clearvars y_ticks;
end

[found, y_ticklabels, varargin] = argParse(varargin, 'yticklabels');
if ~found %|| isempty(y_ticklabels)
    clearvars y_ticklabels;
end

[found, square_flag, varargin] = argParse(varargin, 'square');
if ~found; square_flag = true; end

[found, axis_spec, varargin] = argParse(varargin, 'axis_spec');
if ~found; axis_spec = 'tight'; end

[found, cmap, varargin] = argParse(varargin, 'cmap');
if ~found
    cmap = lines(numel(RG));
end
%--------------------------------------------------------------------------
Nrgs = numel(RG);
axlbls = {'a','p';'v','d';'m','l'};

if isfield(RG{1},'qMRI')
    qMRI = RG{1}.qMRI;
elseif isfield(RG{1},'parameter')
    qMRI = RG{1}.parameter;    
else
    qMRI = '';
end
if isfield(RG{1},'units')
    units = RG{1}.units;
else
    units = '';
end

%--------------------------------------------------------------------------
%% FIG
if size(PC)==1
    pc=PC;
    j=1;
    h=make_plot(pc);
else
    h=figure;
    j=0;
    for pc = PC
        j=j+1;
        subplot(1,length(PC),find(PC==pc))
        make_plot(pc);
    end
end

if ~exist('boundedline.m','file')
    warning('For better visuzalization download the function ''boundedline.m'' and add to MATLAB path (https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m)');
end

%% FUNCTIONS
%==========================================================================
% BOUNDED LINES
%==========================================================================
    function h = bl(RG,jj,pc,use_err,cmap,markershapes)
        Nsubs = length(RG{jj}.individual_data);
        if Nsubs==1
            warning('subject group has only N=1 subjects');
        end
%         x = 1:size(RG{jj}.Y{pc},1);
        nsegs = size(RG{jj}.Y{pc},1);
        x = linspace(0,1,nsegs);
        y = MeanFunction(RG{jj}.Y{pc},2,"omitnan");
        if any(isnan(RG{jj}.Y{pc}(:)))
            warning('Ignoring %d NaN values in data (%s axis %d)',nnz(isnan(RG{jj}.Y{pc}(:))),RG{jj}.ROI_label,pc);
        end
        
        around_zero = 0;
        if around_zero
            f = mean(y);
            y = y - f;
        end
                   
        err = {std(RG{jj}.Y{pc},0,2,"omitnan");...
               std(RG{jj}.Y{pc},0,2,"omitnan")/sqrt(size(RG{jj}.Y{pc},2));...
               mad(RG{jj}.Y{pc},0,2)};
           
        Transparency = .25;
        
        if any(use_err==0) || strcmpi(use_err,'none')
            err = zeros(size(err{1}));
            Transparency = 0;
        elseif any(use_err==1)||strcmpi(use_err,'STD')
            err = err{1};
        elseif any(use_err==2)||strcmpi(use_err,'SEM')
            err = err{2};
        elseif any(use_err==3)||strcmpi(use_err,'MAD')
            err = err{3};
        end
        
        bl_shape = '-';
        
%         if ~markershapes % if markershapes not specified choose -o
%             bl_shape = '-o';
%         else
%             bl_shape = '-';
%         end
        
        if exist('boundedline.m','file')==2 && Nsubs > 1
            h = boundedline(x,y,double(err),...
                bl_shape,'alpha','cmap',cmap(jj,:),'transparency',Transparency);
            h.LineWidth = 2;
            h.Color = cmap(jj,:);
            hold on;
        else
            h = errorbar(x,y,err,'-','color',cmap(jj,:));
            h.LineWidth = 2;
            hold on;   
        end
        
        if isa(markershapes,'cell')
            shapes = markershapes;
            markershapes = true;
        else
            shapes = {'-o','-^','-s','-*','-+','-x','-d','-p','-h'};
        end
        if markershapes
            
            
            if jj> length(shapes)
                h = plot(x,y,'-');
            else
                h = plot(x,y,shapes{jj});
            end
            
            h.MarkerSize = 6;
            h.MarkerEdgeColor = cmap(jj,:);
            h.MarkerFaceColor = cmap(jj,:);
            h.LineWidth = 1;
            h.Color = cmap(jj,:);
        end
end
%==========================================================================
% INDIVIDUAL TRACTS
%==========================================================================
function h = individual_tracts(RG,jj,pc,gradnumber)

        nsegs = size(RG{jj}.Y{pc},1);
        X = linspace(0,1,nsegs);
        Y = RG{jj}.Y{pc};
        Y1 = Y(:,:);

        h = plot(X,Y1,'-');
        
        for ss=1:size(Y1,2)
            h(ss).Color = cmap(jj,:);%[.1 .1 .1]*8;
            [~,i] = max(abs(Y1(:,ss)-median(Y1,2)));

            if gradnumber
                text(X(end),Y1(end,ss),num2str(ss)); % number gradients
                text(X(i),Y1(i,ss),num2str(ss),HorizontalAlignment='center',VerticalAlignment='middle'); % number gradients
            end
        end
end

%==========================================================================
% MAIN FUNCTION - MAKE PLOT
%==========================================================================
function ax = make_plot(pc)
p1 = cell(Nrgs,1);
p2 = cell(Nrgs,1);
    %----------------------------------------------------------------------
    % INDIVIDUAL TRACTS ---------------------------------------------------
    if length(ind)<Nrgs
        ind = [ind zeros(1,Nrgs-length(ind))];
    end
    for jj = 1:Nrgs
        if ind(jj)
                p1{jj} = individual_tracts(RG,jj,pc,gradnumber);
            hold on;
        end
    end
    %----------------------------------------------------------------------
    % BOUNDED LINES -------------------------------------------------------
    % means +/-1 SD
    if ~noBL
        for jj = 1:Nrgs
            p2{jj} = bl(RG,jj,pc,use_err,cmap,markershapes);
            hold on;
        end
    end
    %----------------------------------------------------------------------
    % INDICATE SPECIFIC INDIVIDUALS ---------------------------------------
%         for jj = 1:Nrgs
%             X = RG{jj}.X{pc};
%             Y = RG{jj}.Y{pc};
%             hold on;
%             if jj ==1
%                 t = 13;
%             else
%                 t=[];
%             end
% 
%             Y2 = Y(:,t);
%             if ~isempty(Y2)
%             p2{jj} = plot(X,Y2,'-');
%             for ss=1:size(Y2,2)
%                 p2{jj}(ss).LineStyle = '--';
%                 p2{jj}(ss).LineWidth = 2;
%                 p2{jj}(ss).Color = cmap(jj,:) * 0.7;
%             end
%             end
%         end
    %----------------------------------------------------------------------
    hold off;
    grid on;
    box on;
    %----------------------------------------------------------------------
    % LEGEND --------------------------------------------------------------
    if leg && pc==PC(end) % legend
        
        if ~exist('leg_ent','var')

            % Define legend entries based on descriptive fields that differ
            % between RG cells
            
            % find char fields
            field_names = fieldnames(RG{1});
            char_fields = cellfun(@(x) isa(RG{1}.(x),"char") || isa(RG{1}.(x),"string") && length(RG{1}.(x))==1, field_names);
            char_fields = field_names(char_fields);
            % find char fields that differ between cells
            legend_fields = cellfun(@(y) length(unique(cellfun(@(x) string(x.(y)), RG)))>1 ,char_fields);
            legend_fields = string(char_fields(legend_fields));

            if ~isempty(legend_fields)
                group_n = cellfun(@(x) length(x.individual_data),RG);
                grouping_strings = cellfun(@(rg) strjoin(arrayfun(@(f) string(rg.(f)),legend_fields),", "),RG);
                leg_ent = grouping_strings + arrayfun(@(n) sprintf(" (N=%d)",n),group_n);
            end
        end
                
        if ~exist('leg_ent','var')
            leg_ent = arrayfun(@(x) sprintf('func_%d',x),1:Nrgs,'un',0);
        end
        l=legend([p2{:}],leg_ent);
        
        l.TextColor = [89 89 89]/255;
        l.Interpreter = 'none';
        l.Location = 'northeast';
    end

    %----------------------------------------------------------------------
    % Y-LABEL -------------------------------------------------------------
    if lbl_flag && pc==PC(1)
        if ~exist('y_lbl','var')
            y_lbl = sprintf('%s [%s]', qMRI, units);
            y_lbl = strrep(y_lbl,'_',' ');
        end
        h1=ylabel(y_lbl);
        h1.FontSize = 20;
%         h1.Interpreter = 'None';
    end
    ax = gca;
    ax.FontSize = 26;

    %----------------------------------------------------------------------
    % X Tick Labels (instead of X label) ----------------------------------
    Option = 2;
    
    if Option == 1
        nsegs = size(RG{1}.Y{pc},1);
        xticks(linspace(0,1,nsegs))
%         xticks(1:nsegs);
        xticklabels([]);
        if lbl_flag
            c = arrayfun(@(x) '', 1:nsegs,'UniformOutput',false);
            c{1} = axlbls{pc,1};
            c{end} = axlbls{pc,2};
            c{floor(nsegs/2)+1} = '[segments]';
            xticklabels(c);
        end
    elseif Option == 2
        nsegs = size(RG{1}.Y{pc},1);
        xticks(linspace(0,1,nsegs))
%         xticks(1:nsegs);
        xticklabels([]);
        if lbl_flag
            if exist('x_lbl','var')
                x_lbl = string(x_lbl);
                if length(x_lbl)>1
                    xlabel(x_lbl{pc});
                else
                    xlabel(x_lbl);
                end                    
            else
            x_labels = {'axis 1'; 'axis 2'; 'axis 3'};
            xlabel(x_labels{pc});
            end
        end
    end
    % ---------------------------------------------------------------------

    axis tight;
    if square_flag
        axis square;
    end
    axis(axis_spec);

    %----------------------------------------------------------------------
    % Y-LIM , Y-TICKS, Y-TICK LABELS---------------------------------------
    if exist('y_lim','var')
        if size(y_lim,1) > 1
            current_ylim = y_lim(pc,:);
        else
            current_ylim = y_lim;
        end
        ylim(current_ylim);
    end
    
    if exist('y_ticks','var')
        if iscell(y_ticks)
            current_yticks = y_ticks{pc,:};
        else
            current_yticks = y_ticks;
        end
        yticks(current_yticks);
    else
        current_yticks = yticks;
    end
    
    if exist('y_ticklabels','var')
        if iscell(y_ticklabels)
            current_yticklabels = y_ticklabels{pc,:};
        else
            current_yticklabels = y_ticklabels;
        end
        tmp = round(current_yticks,5);
        tmp(~ismember(tmp,round(current_yticklabels,5)))=nan;
        tmp = arrayfun(@(x) num2str(x),tmp,'UniformOutput',false);
        tmp(strcmp(tmp,'NaN')) = {''};
        
        if j~=1
            tmp(:) = {''};
        end
        
        current_yticklabels = tmp;
        yticklabels(current_yticklabels);
    end
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
end
end