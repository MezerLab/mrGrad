function tbl = mrGrad_output2table(RG,Axes)
% A table for linear mixed-effects modelling (fitlme input)
% -------------------------------------------------------------------------
% INPUTS:
% 
%   RG        mrGrad output with one or more groups and one or more ROIs
%   Axes      one or more ROI's axes [1-3]
% -------------------------------------------------------------------------
% OUTPUT:
% 
%   tbl     a table with variables {SubjectID, ClinicalGroup, Sex, Age,
%           RoiName, Axis, Position, Response}
% -------------------------------------------------------------------------
% E.D.
% -------------------------------------------------------------------------

tbl = table;
for gg = 1:size(RG,1)
    for rr = 1:size(RG,2)
        rg = RG{gg,rr};
        for ax = Axes
            nsegs = size(rg.Y{ax},1);
            nsubs = size(rg.Y{ax},2);
            for ii = 1:nsubs
                t = table;
                t.SubjectID = repmat(rg.subject_names(ii),nsegs,1);
                t.ClinicalGroup = repmat({rg.group_name},nsegs,1);
                t.Sex = repmat(rg.sex(ii),nsegs,1);
                t.Age = repmat(rg.age(ii),nsegs,1);
                t.RoiName = repmat({rg.ROI_label},nsegs,1);
                t.Axis = repmat(ax,nsegs,1);
    %             t.hemisphere_id = double(startsWith(t.RoiName,'Right'));
                t.Position = linspace(0,1,nsegs)';
%                 t.Position = (1:nsegs)';
                t.Response = rg.Y{ax}(:,ii);
                tbl = [tbl;t];
            end
        end
    end
end
tbl.SubjectID_original = string(tbl.SubjectID);
tbl.SubjectID = grp2idx(tbl.SubjectID);
tbl.Hemisphere = cellfun(@(x) x(1), tbl.RoiName,'un',0);
% tbl.ClinicalGroup = grp2idx(tbl.ClinicalGroup)-1;
% tbl.Hemisphere = double(startsWith(tbl.RoiName,'Right'));
% tbl.Sex = grp2idx(tbl.Sex)-1;

tbl = tbl(:,{'SubjectID_original','SubjectID','ClinicalGroup','Sex','Age','Hemisphere','RoiName','Axis','Position','Response'});

tbl.SubjectID = categorical(tbl.SubjectID);
tbl.ClinicalGroup = categorical(tbl.ClinicalGroup);
tbl.Hemisphere = categorical(tbl.Hemisphere);
tbl.RoiName = categorical(tbl.RoiName);
tbl.Sex = categorical(tbl.Sex);
tbl.Axis = categorical(tbl.Axis);

if numel(unique(tbl.RoiName)) > 2
    warning('more than one ROI in table:');
    disp(unique(tbl.RoiName));
end

