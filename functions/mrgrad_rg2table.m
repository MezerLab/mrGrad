function T = mrgrad_rg2table(RG,axis)

if ~isa(RG,"cell")
    RG = {RG};
end

T = cell(size(RG));

if ~exist('axis','var')
    axis = [];
end

% -------------------------------------------------------------------------
% get descriptive fields (same for all rgs)
n_observations = length(RG{1}.subject_names);
s = RG{1};
field_names = fieldnames(s);
if isfield(field_names,'user_input_fields')
    s = s.user_input_fields;
    field_names = fieldnames(s);
end

mrgrad_reserved_fields = ["Y", "Y_mean", "Y_std", "Y_SEM", "X", ...
    "N_segments", "parameter", "units", "sampling_method", ...
    "method", "y_lbls", "ROI_label", "individual_data", "group_name", "subject_names"];

field_names = setdiff(field_names,mrgrad_reserved_fields,'stable');

is_descrip = arrayfun(@(v) mrgrad_valid_desciptive_field(s.(v)), field_names) ...
    & cellfun(@(x) length(s.(x))==n_observations  & ~ischar(s.(x)), ...
        field_names);

sub_descrips = field_names(is_descrip);
% -------------------------------------------------------------------------

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

    % add descriptive fields
%     n_observations = length(rg.subject_names);
%     s = rg;
%     field_names = fieldnames(s);
%     if isfield(field_names,'user_input_fields')
%         s = s.user_input_fields;
%         field_names = fieldnames(s);
%     end
% 
%     mrgrad_reserved_fields = ["Y", "Y_mean", "Y_std", "Y_SEM", "X", ...
%         "N_segments", "parameter", "units", "sampling_method", ...
%         "method", "y_lbls", "ROI_label", "individual_data", "group_name", "subject_names"];
% 
%     field_names = setdiff(field_names,mrgrad_reserved_fields,'stable');
% 
%     is_descrip = arrayfun(@(v) mrgrad_valid_desciptive_field(s.(v)), field_names) ...
%         & cellfun(@(x) length(s.(x))==n_observations  & ~ischar(s.(x)), ...
%             field_names);
% 
%     sub_descrips = field_names(is_descrip);

    % add descriptive fields
    for ii = 1:length(sub_descrips)
        t.(sub_descrips{ii}) = rg.(sub_descrips{ii});
    end
    t = movevars(t,sub_descrips,"After","GroupName");
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
    

    % Make sure all ID columns and descriptive fields are identical, remove
    % them, and join all columns (all ROIs)
    vars = ["SubjectName", "GroupName", sub_descrips'];
    allEqual = all(cellfun(@(x) isequal(x(:,vars), T{1}(:,vars)), T));
    if allEqual
        T(2:end) = cellfun(@(t) removevars(t,vars),T(2:end),'un',0);
    end

    T = horzcat(T{:});
end