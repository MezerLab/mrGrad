function tbls = mrGrad_for_lme_modeling(RG,Axes, Grouping)
% A table for linear mixed-effects modelling (fitlme input)
% -------------------------------------------------------------------------
% INPUTS:
% 
%   RG        mrGrad output with one or more groups and one or more ROIs
%               either .mat or .csv mrgrad output
%   Axes      one or more ROI's axes [1-3]
% -------------------------------------------------------------------------
% OUTPUT:
% 
%   tbl     a table with variables {SubjectID, ClinicalGroup, Sex, Age,
%           RoiName, Axis, Position, Response}
% -------------------------------------------------------------------------
% E.D.
% -------------------------------------------------------------------------

if isa(RG,'string') || isa(RG,'char')
    if endsWith(RG,'.csv')
        RG = readtable(RG,ReadVariableNames=true,delimiter=",");
    elseif endsWith(RG,'.mat')
        RG = load(RG).RG;
    end
end

if isa(RG,'table')
    T = RG;
else
    T = mrgrad_rg2table(RG);
end

if exist('Grouping', 'var') && ~isempty(Grouping)
    if length(Grouping)==height(T)
        T.('GroupName') = categorical(Grouping);
    else
        warning('Grouping variable doesn''t match table height. Ignoring.');
    end
end

if ~ismember('subject_id', T.Properties.VariableNames)
    T.subject_id = T.Properties.RowNames;
end

for v = {'seg_list','map_list'}
    idx_rm = startsWith(T.Properties.VariableNames, v);
    T(:,idx_rm) = [];
end


T.Properties.RowNames = {};
T_all = T;

% separate tables for each image parameter (response variable)
expr = '^(.*)_axis\d+_seg\d+_(\w+)$';
tokens = regexp(T_all.Properties.VariableNames', expr, 'tokens', 'once');
params = cellfun(@(t) t{2}, tokens(~cellfun('isempty', tokens)), 'un', 0);
parameter_names = unique(params, 'stable');

formattedMask = ~cellfun('isempty', tokens);

% parameters = cellfun(@(x) regexp(x, '.*_axis\d+_seg\d+_(.*)$', 'tokens', 'once'), T.Properties.VariableNames','un',0);
% parameters = 
% parameter_names = unique(parameters, 'stable');
tbls = struct;
for jj = 1:length(parameter_names)
    param = string(parameter_names{jj});

    pexpr = '^(.*)_axis\d+_seg\d+_' + param + '$';
    idx = ~cellfun(@isempty, regexp(T_all.Properties.VariableNames', pexpr, 'once'));
    idx = idx | ~formattedMask;

    T = T_all(:,idx);

% convert mrgrad table to a long table for fitlme

% 1) Pick the measurement columns ("..._axis#_seg#")
vn      = T.Properties.VariableNames;
isMeas  = contains(vn,'_axis') & contains(vn,'_seg') & endsWith(vn,"_"+param);

% 2) Stack wide -> long
%    - 'value' will hold the numeric measurements
%    - 'var'   will hold the original variable name like "Right_Putamen_axis1_seg8"
longT = stack(T, vn(isMeas), ...
    'NewDataVariableName','Response', ...
    'IndexVariableName','var');

longT.var = regexprep(string(longT.var), '_'+param+'$', '');

% 3) Parse roi / axis / seg from the original variable name
%    Works for any ROI text before "_axis"
m = regexp(string(longT.var), '^(?<roi>.+?)_axis(?<axis>\d+)_seg(?<seg>\d+)$', 'names');
m = cat(1,m{:});
longT.RoiName  = [m.roi]';
longT.Axis = str2double([m.axis]');
longT.Position  = str2double([m.seg]');

% 4) Clean up and reorder columns
longT.var = [];



mrgrad_vars  = {'RoiName','Axis','Position','Response'};
other_vars = setdiff(longT.Properties.VariableNames, mrgrad_vars, 'stable');
if ~isempty(other_vars)
    longT = movevars(longT, mrgrad_vars, 'After',other_vars{end});
end

% 5) Convert position to 0-1 for each (roi x axis) combination
roi_names = unique(longT.RoiName);
for rr = 1:length(roi_names)
    idx_roi = longT.RoiName==roi_names(rr);
    axis_vals = unique(longT.Axis(idx_roi));
    for aa = 1:length(axis_vals)
        idx_axis = longT.Axis==axis_vals(aa);
        idx = idx_roi & idx_axis;
        longT.Position(idx) = longT.Position(idx)./max(longT.Position(idx));
    end
end

% % 6) Add Hemisphere variable (assuming roi names start with 'L' / 'R'
% longT.Hemisphere = cellfun(@(x) x(1), cellstr(longT.RoiName));

% 7) Convert GroupName and other variables if needed
% longT.GroupName = categorical(string(longT.GroupName));
varNames = {'GroupName', 'Sex', 'Hemisphere'};
for ii = 1:length(varNames)
    if ismember(varNames{ii},other_vars)
        longT.(varNames{ii}) = categorical(longT.(varNames{ii}));
    end
end

if exist('Axes','var')
    longT(~ismember(longT.Axis, Axes),:) = [];
end

% longT.SubjectName = grp2idx(longT.SubjectName);

tbl = longT;
tbl.Response = double(tbl.Response);
if numel(unique(tbl.RoiName)) > 1
    warning('more than one ROI in table:');
    disp(unique(tbl.RoiName));
end

tbls.(param) = tbl;
end

return

