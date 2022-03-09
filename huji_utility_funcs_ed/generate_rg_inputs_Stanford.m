function [groups_data, group_names] = generate_rg_inputs_Stanford(param,varargin)
%--------------------------------------------------------------------------
% chosen subjects from WH + Gotlib to match huji
%--------------------------------------------------------------------------
[found, noedge, varargin] = parseparam(varargin, 'noEdge');
if ~found; noedge = true; end


group_names = {'young','old'};
groups_data = cell(size(group_names));



Only_WH = 1;


gg = 1; % YOUNG
group = group_names{gg};

    if ~Only_WH
        warning('Using Gotlib Subjects to match HUJI age and sex');
    %----------------------------------------------------------------------
    % This list is optimized to match huji age and sex
    % (see sample_match_huji_wh.m)
    subjects = {'002_KM';'009_AL';'012_JR';'041_LL';...
                '2614';'2722';'2770';'2796';'2825';'2830';'2841';'2848';...
                '2852';'2854';'2908';'2933';'2945';'2955';'2967';'3039'};
    age = [31;25;31;24;28;24;29;27;25;25;29;26;23;29;24;22;30;32;25;26];
    sex = {'M';'F';'M';'M';'F';'F';'F';'F';'M';'F';'F';'M';'M';'M';'F';'M';'M';'M';'F';'F'};

    [map_list, seg_list] = stanfordWH_MapAndSeg(subjects(1:4),param,'noedge',noedge);
    tv_list = stanfordWH_MapAndSeg(subjects(1:4), 'MTV','noedge',noedge);

    [map_list1, seg_list1] = Gotlib_MapAndSeg(subjects(5:end),param,'noedge',noedge);
    tv_list1 = Gotlib_MapAndSeg(subjects(5:end), 'MTV','noedge',noedge);

    map_list = [map_list; map_list1];
    seg_list = [seg_list; seg_list1];
    tv_list = [tv_list; tv_list1];
    
    elseif Only_WH
        [subjects,age,sex] = stanfordWH_subjects('over',18,'under',32);
        [map_list, seg_list] = stanfordWH_MapAndSeg(subjects,param,'noedge',noedge);
        tv_list = stanfordWH_MapAndSeg(subjects, 'MTV','noedge',noedge);
    end

    groups_data{gg}.group_name = group;
    groups_data{gg}.subjects = subjects;
    groups_data{gg}.age = age;
    groups_data{gg}.sex = sex;
    groups_data{gg}.map_list = map_list;
    groups_data{gg}.seg_list = seg_list;
    groups_data{gg}.tv_list = tv_list;    

gg = 2; % OLD
    group = group_names{gg};
    [subjects,age,sex] = stanfordWH_subjects(group);
    
    [map_list, seg_list] = stanfordWH_MapAndSeg(subjects,param,'noedge',noedge);
    tv_list = stanfordWH_MapAndSeg(subjects,'MTV','noedge',noedge);

    if ~Only_WH
        [subjects1,age1,sex1] = stanfGotlib_subjects(group,'CTL');
        [map_list1, seg_list1] = Gotlib_MapAndSeg(subjects1,param,'noedge',noedge);
        tv_list1 = Gotlib_MapAndSeg(subjects1,'MTV','noedge',noedge);

        subjects = [subjects;subjects1];
        age = [age; age1];
        sex = [sex; sex1];
        map_list = [map_list; map_list1];
        seg_list = [seg_list; seg_list1];
        tv_list = [tv_list; tv_list1];
    end
    groups_data{gg}.group_name = group;
    groups_data{gg}.subjects = subjects;
    groups_data{gg}.age = age;
    groups_data{gg}.sex = sex;
    groups_data{gg}.map_list = map_list;
    groups_data{gg}.seg_list = seg_list;
    groups_data{gg}.tv_list = tv_list;    
    
%--------------------------------------------------------------------------
%% use all WH + gotlib subjects
%--------------------------------------------------------------------------

% add_gotlib_CTL = true;
% 
% group_names = {'young','old'};
% groups_data = cell(size(group_names));
% for gg = 1:2
%     group = group_names{gg};
%     [subjects,age,sex] = stanfordWH_subjects(group);
%     Nsubs = length(subjects);
%     
%     [map_list, seg_list] = stanfordWH_MapAndSeg(subjects,param);
%     tv_list = stanfordWH_MapAndSeg(subjects, 'MTV');
%     
%     if add_gotlib_CTL
%         warning('adding Gotlib''s control subjects');
%         [subjects1,age1,sex1] = stanfGotlib_subjects(group,'CTL');
%         [map_list1, seg_list1] = Gotlib_MapAndSeg(subjects1,param);
%         tv_list1 = Gotlib_MapAndSeg(subjects1, 'MTV');
%         
%         subjects = [subjects;subjects1];
%         age = [age; age1];
%         sex = [sex; sex1];
%         map_list = [map_list; map_list1];
%         seg_list = [seg_list; seg_list1];
%         tv_list = [tv_list; tv_list1];
%     end
% 
%     groups_data{gg}.subjects = subjects;
%     groups_data{gg}.age = age;
%     groups_data{gg}.sex = sex;
%     groups_data{gg}.map_list = map_list;
%     groups_data{gg}.seg_list = seg_list;
%     groups_data{gg}.tv_list = tv_list;    
% 
% end

