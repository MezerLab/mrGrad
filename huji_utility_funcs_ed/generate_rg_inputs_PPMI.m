function [groups_data, group_names] = generate_rg_inputs_PPMI(param,varargin)   
group_names = {'Control','PD'};
groups_data = cell(size(group_names));

analysisDir = '/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/PPMI_trio';

for gg = 1:2
    group = group_names{gg};
%     [subjects,age,sex] = PPMI_subjects(group);
    [subjects,age,sex] = PPMI_subjects(group,'over',55,'under',76);

    map_list = cellfun(@(x) fullfile(analysisDir,x,'1mm',[param,'.nii.gz']),subjects,'un',0);
    seg_list = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST/first_all_none_firstseg_noEdge.nii.gz'),subjects,'un',0);

    groups_data{gg}.group_name = group;
    groups_data{gg}.subjects = subjects;
    groups_data{gg}.age = age;
    groups_data{gg}.sex = sex;
    groups_data{gg}.map_list = map_list;
    groups_data{gg}.seg_list = seg_list;
    groups_data{gg}.tv_list = [];
end
