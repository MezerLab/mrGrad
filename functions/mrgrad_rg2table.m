function T = mrgrad_rg2table(RG,axis)

if ~isa(RG,"cell")
    RG = {RG};
end

T = cell(size(RG));

if ~exist('axis','var')
    axis = [];
end

for jj = 1:numel(RG)
    rg = RG{jj};

    if isempty(axis)
        rg_ax = 1:length(rg.Y);
    else
        rg_ax = axis;
    end

    t = array2table(cat(1, rg.Y{rg_ax})');
    
    axis_str = cellstr(string(rg.ROI_label) +"_" + rg.y_lbls(rg_ax));
    seg_str = cellfun(@(y) "_seg"+ string(1:size(y,1))',rg.Y(rg_ax),'un',0);
    t_colnames = cellfun(@(a,b) a+b, axis_str,seg_str,'un',0);
    t_colnames = cat(1,t_colnames{:});
    t.Properties.VariableNames = t_colnames;
    t.GroupName = repmat(string(rg.group_name),height(t),1);
    t = movevars(t,"GroupName","Before",t.Properties.VariableNames{1});
    t.SubjectName = string(rg.subject_names);
    t = movevars(t,"SubjectName","Before","GroupName");

    T{jj} = t;
end

% concatenate all tables
if numel(T)==1
    T = T{1};
else
    % join all rows (subjects of all groups)
    h = height(T)+1;
    for j = 1:width(T)
        T{h,j} = vertcat(T{:, j});
    end
    T = T(end,:);
    
    % make sure all ID columns are identical before joining all columns (all ROIs)
    
    allEqual = all(cellfun(@(x) isequal(x(:,{'SubjectName','GroupName'}), T{1}(:,{'SubjectName','GroupName'})), T));
    if allEqual
        T(2:end) = cellfun(@(t) removevars(t,{'SubjectName','GroupName'}),T(2:end),'un',0);
    end
    T = horzcat(T{:});
end