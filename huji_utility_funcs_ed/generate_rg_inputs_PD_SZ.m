function [groups_data, group_names] = generate_rg_inputs_PD_SZ(param,varargin)
    
    group_names = {'PD','CTL','DYT','GBA'};
    groups_data = cell(size(group_names));
    for gg = 1:length(group_names)
        group = group_names{gg};
        [subjects,age,sex] = PD_SZ_subjects(group,param);
        Nsubs = length(subjects);
        
        n = Nsubs; % chosen subjects
        subjects = subjects(1:n);
        age = age(1:n);
        sex = sex(1:n);
        
        [map_list, seg_list] = PD_SZ_MapAndSeg(subjects,param,varargin{:});
%         tv_list = PD_SZ_MapAndSeg(subjects, 'MTV');
        
        groups_data{gg}.group_name = group;
        groups_data{gg}.subjects = subjects;
        groups_data{gg}.age = age;
        groups_data{gg}.sex = sex;
        groups_data{gg}.map_list = map_list;
        groups_data{gg}.seg_list = seg_list;
%         groups_data{gg}.tv_list = tv_list;
    end

%     %% add HUJI Old as Control (mrQ_fixbias)
%     gg = gg+1;
%     [subjects,age,sex] = HUJI_subjects('old',param);    
%     [map_list, seg_list] = huji_MapAndSeg(subjects,param,'noedge',noEdge,'mrQver',2);
%     groups_data{gg}.group_name = 'Old HC';
%     groups_data{gg}.subjects = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
%     
%     %% add HUJI Young as Control (mrQ_fixbias)
%     gg = gg+1;
%     [subjects,age,sex] = HUJI_subjects('young',param);    
%     [map_list, seg_list] = huji_MapAndSeg(subjects,param,'noedge',noEdge,'mrQver',2);
%     groups_data{gg}.group_name = 'Young HC';
%     groups_data{gg}.subjects = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
% 
%     %% add HUJI Old as Control (mrQ_fixbias - New)
%     gg = gg+1;
%     [subjects,age,sex,analysisDir] = HUJI_subjects('old',param);
%     map_list = cellfun(@(x) fullfile(analysisDir,x,'mrQ_tests_ed/mrQ_fixbias/OutPutFiles_1/BrainMaps/R1_map.nii.gz'),subjects,'un',0);
%     seg_list = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST_v2/first_all_none_firstseg.nii.gz'), subjects,'un', 0);
%     
%     groups_data{gg}.group_name = 'Old HC new';
%     groups_data{gg}.subjects = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
%     
%     %% add HUJI Young as Control (mrQ_fixbias - New)
%     gg = gg+1;
%     [subjects,age,sex,analysisDir] = HUJI_subjects('young',param);
%     map_list = cellfun(@(x) fullfile(analysisDir,x,'mrQ_tests_ed/mrQ_fixbias/OutPutFiles_1/BrainMaps/R1_map.nii.gz'),subjects,'un',0);
%     if noEdge
%         seg_list = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST_v2/first_all_none_firstseg_noEdge.nii.gz'), subjects,'un', 0);
%     else
%         seg_list = cellfun(@(x) fullfile(analysisDir,x,'FSLFIRST_v2/first_all_none_firstseg.nii.gz'), subjects,'un', 0);
%     end
%     
%     groups_data{gg}.group_name = 'Young HC new';
%     groups_data{gg}.subjects = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
%     
%     
%     group_names = [group_names,'Old HC','Young HC','Old HC new','Young HC new'];
    
end
    