function [groups_data, group_names] = generate_rg_inputs_HCP(param,varargin)

if any(strcmpi(param,{'T2woverT1w','T2overT1'}))
    param = 'T2wDividedByT1w';
elseif any(strcmpi(param,{'T1woverT2w','T1overT2'}))
    param = 'T1wDividedByT2w';
end

group_names = {'young'};
groups_data = cell(size(group_names));
gg=1;
group = group_names{gg};
    [subjects,age,sex] = HCP_subjects;
    Nsubs = length(subjects);
    
    n = Nsubs; % chosen subjects
    subjects = subjects(1:n);
    age = age(1:n);
    sex = sex(1:n);
    
    [map_list, seg_list, subjects, idx] = HCPY_MapAndSeg(subjects, param,'noEdge',varargin{:});
    age(idx)=[];
    sex(idx)=[];

    groups_data{gg}.group_name = group;
    groups_data{gg}.subjects = subjects;
    groups_data{gg}.age = age;
    groups_data{gg}.sex = sex;
    groups_data{gg}.map_list = map_list;
    groups_data{gg}.seg_list = seg_list;
    groups_data{gg}.tv_list = [];    
