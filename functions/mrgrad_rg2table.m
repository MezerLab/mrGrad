function T = mrgrad_rg2table(RG)

rg_fields = fieldnames(RG);
T = cell(1,length(rg_fields));


% -------------------------------------------------------------------------
% get subject descriptive fields (same for all rgs)
n_observations = length(RG.(rg_fields{1}).subject_ids);
s = RG.(rg_fields{1}).user_input_fields;
field_names = fieldnames(s);

sub_descrips = field_names(...
    structfun(@(f) ~ischar(f) & length(f)==n_observations, s) ...
    );
% -------------------------------------------------------------------------

for jj = 1:length(rg_fields)
    rg = RG.(rg_fields{jj});

    tbl_names = fieldnames(rg.Results);
    for tt = 1:length(tbl_names)

        base = tbl_names{tt};
        segnames = rg.Results.(tbl_names{tt}).y.Properties.VariableNames;

        % add seg names after axis names
        [~,~, endIdx] = regexp(base, rg.axis_names, 'match', 'start', 'end');
        axis_ind = ~cellfun(@isempty, endIdx);
        endIdx = endIdx{axis_ind};

        roi_name = regexprep(lower(rg.roi_name),'[\\/:*?"<>|-]|[\x00-\x1F]','_');
        varnames = roi_name+"_"+base(1:endIdx)+"_"+segnames + base(endIdx+1:end);    
    
        rg.Results.(tbl_names{tt}).y.Properties.VariableNames = varnames;
    end
    t = cellfun(@(t) rg.Results.(t).y, tbl_names, 'un', 0);
    t = [t{:}];

    % add descriptive fields
    for ii = 1:length(sub_descrips)
        t.(sub_descrips{ii}) = s.(sub_descrips{ii});
    end
    t = movevars(t,sub_descrips,"Before", t.Properties.VariableNames{1});
    T{jj} = t;
end

% concatenate all tables
if numel(T)==1
    T = T{1};
else
%     % join all rows (subjects of all groups) -- mrGrad v1.2
%     h = height(T)+1;
%     for j = 1:width(T)
%         T{h,j} = vertcat(T{:, j});
%     end
%     T = T(end,:);
    
    % Make sure all ID columns and descriptive fields are identical, remove
    % them, and join all columns (all ROIs)
    vars = sub_descrips';
    allEqual = all(cellfun(@(x) isequal(x(:,vars), T{1}(:,vars)), T));
    if allEqual
        T(2:end) = cellfun(@(t) removevars(t,vars),T(2:end),'un',0);
    end

    T = horzcat(T{:});

    % add subject id column
    T.subject_id = T.Properties.RowNames;
    T = movevars(T,"subject_id","Before", T.Properties.VariableNames{1});

    % split map_list if more than one parameter (as of mrGrad v-2)
    n_parameters = size(T.map_list,2);
    if n_parameters > 1
        for pp = 1:n_parameters
            T.("map_list"+pp) = T.map_list(:,pp);
        end
        T = movevars(T,"map_list"+(1:n_parameters),"After","map_list");
        T.map_list = [];
    end

end