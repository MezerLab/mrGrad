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
%   'cmap'      specify colors (one RGB row per function)
%
%   'square' 0 or 1(default) for square figure axis
%
%   'axis_spec' default is 'tight'
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

[found, ind, varargin] = argParse(varargin, 'ind');
if ~found; ind = zeros(size(RG)); end
if numel(ind)==1
    ind = ind*ones(size(RG));
end

[found, PC, varargin] = argParse(varargin, 'PC');
if ~found; PC = 1:3; end

[found, use_err, varargin] = argParse(varargin, 'errorType');
if ~found; use_err = 1; end


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
if ~found || isempty(y_lbl)
    clearvars y_lbl;
end

[found, x_lbl, varargin] = argParse(varargin, 'xlbl');
if ~found || isempty(x_lbl)
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

% what groups are we dealing with? (for legend)
if isfield(RG{1},'ROI_label')
    leg_roi=1;
else
    leg_roi=0;
end
if isfield(RG{1},'group_name')
    leg_group=1;
%     group_diff = numel(unique(cellfun(@(x) x.group_name,RG,'un',0))) >1;
%     if group_diff
%         leg_group=1;
%     else
%         leg_group=0;
%     end
else
    leg_group=0;
end




% if isfield(RG{1},{'group_name', 'ROI_label'})
%     if numel(RG)==1
%         group_diff=1;
%         roi_diff=1;
%     else
%         group_diff = numel(unique(cellfun(@(x) x.group_name,RG,'un',0))) >1;
%         roi_diff = numel(unique(cellfun(@(x) x.ROI_label,RG,'un',0))) > 1;
%     end
% else
%     group_diff=0;
%     roi_diff=0;
% end

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
        x = RG{jj}.X{pc};
        y = RG{jj}.Y_mean{pc};
        
        around_zero = 0;
        if around_zero
            f = mean(y);
            y = y - f;
        end
        
%         err = {RG{jj}.Y_std{pc};...
%                RG{jj}.Y_SEM{pc}};
           
        err = {std(RG{jj}.Y{pc},0,2);...
               std(RG{jj}.Y{pc},0,2)/sqrt(size(RG{jj}.Y{pc},2))};
           
        Transparency = .25;
        
        if any(use_err==0) || strcmpi(use_err,'none')
            err = zeros(size(err{1}));
            Transparency = 0;
        elseif any(use_err==1)||strcmpi(use_err,'STD')
            err = err{1};
        elseif any(use_err==2)||strcmpi(use_err,'SEM')
            err = err{2};
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

            X = RG{jj}.X{pc};
            Y = RG{jj}.Y{pc};
            Y1 = Y(:,:);

            h = plot(X,Y1,'-');
            
            for ss=1:size(Y1,2)
                h(ss).Color = cmap(jj,:);%[.1 .1 .1]*8;
                if gradnumber
                    text(X(end),Y1(end,ss),num2str(ss)); % number gradients
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
        if leg_roi && leg_group
            gr = cellfun(@(x) x.group_name,RG,'un',0);
            group_n = cellfun(@(x) num2str(length(x.individual_data)),RG,'un',0);
            gr = cellfun(@(a,b) [a,' (N=',b,')'],gr,group_n,'un',0);
            ri = cellfun(@(x) x.ROI_label,RG,'un',0);
            leg_ent = cellfun(@(a,b) [a,': ',b],ri,gr,'un',0);
        elseif leg_roi
            ri = cellfun(@(x) x.ROI_label,RG,'un',0);
            leg_ent = ri;
        elseif leg_group
            gr = cellfun(@(x) x.group_name,RG,'un',0);
            gr = cellfun(@(a,b) [a,' (N=',b,')'],gr,group_n,'un',0);
            leg_ent = gr;
        end
        end
        
%         if any([group_diff,roi_diff])
%             if ~group_diff
%                 ri = cellfun(@(x) x.ROI_label,RG,'un',0);
%                 leg_ent = ri;
%             else
%             gr = cellfun(@(x) x.group_name,RG,'un',0);
%             ri = cellfun(@(x) x.ROI_label,RG,'un',0);
%             
%             leg_ent = cellfun(@(a,b) [a,': ',b],ri,gr,'un',0);
%             
%         end
        
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
        xticks(1:nsegs);
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
        xticks(1:nsegs);
        xticklabels([]);
        if lbl_flag
            if exist('x_lbl','var')
                if ischar(x_lbl)
                    xlabel(x_lbl);
                elseif iscell(x_lbl) && numel(x_lbl)>1
                    xlabel(x_lbl{pc});
                elseif iscell(x_lbl)
                    xlabel(x_lbl{1});
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