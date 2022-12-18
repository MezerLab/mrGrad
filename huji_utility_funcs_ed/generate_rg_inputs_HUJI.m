function [groups_data, group_names] = generate_rg_inputs_HUJI(param,varargin)

% [found,~, varargin] = argParse(varargin, 'bsatlas');
% if found
%     bsatlas_val = 'bsatlas';
% else
%     bsatlas_val = ' ';
% end

group_names = {'young','old'};
groups_data = cell(size(group_names));
for gg = 1:2
    group = group_names{gg};
    [subjects,age,sex] = HUJI_subjects(group,param,varargin{:});
    Nsubs = length(subjects);

    n = Nsubs; % chosen subjects
    subjects = subjects(1:n);
    age = age(1:n);
    sex = sex(1:n);

    [map_list, seg_list] = huji_MapAndSeg(subjects, param,varargin{:});
%         tv_list = huji_MapAndSeg(subjects, 'MTV');

    groups_data{gg}.group_name = group;
    groups_data{gg}.subject_names = subjects;
    groups_data{gg}.age = age;
    groups_data{gg}.sex = sex;
    groups_data{gg}.map_list = map_list;
    groups_data{gg}.seg_list = seg_list;
%         groups_data{gg}.tv_list = tv_list;
end
end
    