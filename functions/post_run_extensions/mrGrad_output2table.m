function tbl = mrGrad_output2table(RG,Axes)
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
        T = RG;
    elseif endsWith(RG,'.mat')
        RG = load(RG).RG;
    end
end

if ~isa(RG,'table')
    T = mrgrad_rg2table(RG);
end

% convert mrgrad table to a long table for fitlme

% 1) Pick the measurement columns ("..._axis#_seg#")
vn      = T.Properties.VariableNames;
isMeas  = contains(vn,'_axis') & contains(vn,'_seg');

% 2) Stack wide -> long
%    - 'value' will hold the numeric measurements
%    - 'var'   will hold the original variable name like "Right_Putamen_axis1_seg8"
longT = stack(T, vn(isMeas), ...
    'NewDataVariableName','Response', ...
    'IndexVariableName','var');

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
longT = movevars(longT, mrgrad_vars, 'After',other_vars{end});

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

% 6) Add Hemisphere variable (assuming roi names start with 'L' / 'R'
longT.Hemisphere = cellfun(@(x) x(1), cellstr(longT.RoiName));

% 7) Convert GroupName and other variables if needed
longT.GroupName = categorical(string(longT.GroupName));
varNames = {'Axis', 'GroupName', 'Sex', 'RoiName', 'Hemisphere'};
for jj = 1:length(varNames)
    if ismember(varNames{jj},other_vars)
        longT.(varNames{jj}) = categorical(longT.(varNames{jj}));
    end
end

if exist('Axes','var')
    longT(~ismember(longT.Axis, Axes),:) = [];
end

% longT.SubjectName = grp2idx(longT.SubjectName);

tbl = longT;
if numel(unique(tbl.RoiName)) > 1
    warning('more than one ROI in table:');
    disp(unique(tbl.RoiName));
end

return



% for gg = 1:size(RG,1)
%     for rr = 1:size(RG,2)
%         rg = RG{gg,rr};
%         for ax = Axes
%             nsegs = size(rg.Y{ax},1);
%             nsubs = size(rg.Y{ax},2);
%             for ii = 1:nsubs
%                 t = table;
%                 t.SubjectID = repmat(rg.subject_names(ii),nsegs,1);
%                 t.ClinicalGroup = repmat(string(rg.group_name),nsegs,1);
%                 t.Sex = repmat(rg.sex(ii),nsegs,1);
%                 t.Age = repmat(rg.age(ii),nsegs,1);
%                 t.RoiName = repmat({rg.ROI_label},nsegs,1);
%                 t.Axis = repmat(ax,nsegs,1);
%     %             t.hemisphere_id = double(startsWith(t.RoiName,'Right'));
%                 t.Position = linspace(0,1,nsegs)';
% %                 t.Position = (1:nsegs)';
%                 t.Response = rg.Y{ax}(:,ii);
%                 tbl = [tbl;t];
%             end
%         end
%     end
% end
% tbl.SubjectID_original = string(tbl.SubjectID);
% tbl.SubjectID = grp2idx(tbl.SubjectID);
% tbl.Hemisphere = cellfun(@(x) x(1), tbl.RoiName,'un',0);
% % tbl.ClinicalGroup = grp2idx(tbl.ClinicalGroup)-1;
% % tbl.Hemisphere = double(startsWith(tbl.RoiName,'Right'));
% % tbl.Sex = grp2idx(tbl.Sex)-1;
% 
% tbl = tbl(:,{'SubjectID_original','SubjectID','ClinicalGroup','Sex','Age','Hemisphere','RoiName','Axis','Position','Response'});
% 
% tbl.SubjectID = categorical(tbl.SubjectID);
% 
% tbl.ClinicalGroup = categorical(tbl.ClinicalGroup);
% tbl.Hemisphere = categorical(tbl.Hemisphere);
% tbl.RoiName = categorical(tbl.RoiName);
% tbl.Sex = categorical(tbl.Sex);
% tbl.Axis = categorical(tbl.Axis);
% 
% if numel(unique(tbl.RoiName)) > 2
%     warning('more than one ROI in table:');
%     disp(unique(tbl.RoiName));
% end

