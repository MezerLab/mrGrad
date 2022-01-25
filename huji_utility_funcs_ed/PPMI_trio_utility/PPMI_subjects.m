function [subjects,age,sex,group,analysisDir] = PPMI_subjects(varargin)
% this function returns PPMI subjects from a subgroup of individuals (PD
% and Control) that went through the same T1w, T2w, DTI protocols in
% SIEMENS TRIO scanner. for each subject one session was chosen s.t. the
% mean age the two clinical groups is optimized to match.

analysisDir = '/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/PPMI_trio';
meta = readtable(fullfile(analysisDir,'subjects_chosen.csv'),'Delimiter',',');

[found, overval, varargin] = argParse(varargin, 'OVER');
if ~found; overval = 0; end
[found, underval, varargin] = argParse(varargin, 'UNDER');
if ~found; underval = 120; end


varargin = upper(varargin);
varargin = cellfun(@(x) strrep(x,'CTL','CONTROL'),varargin,'un',0);
varargin = cellfun(@(x) strrep(x,'HC','CONTROL'),varargin,'un',0);

% remove corrupt subjects
if ~ismember('ALL',varargin)
    meta(meta{:,'PROBLEM'}==1,:)=[];
end

% group
idx_group = ismember(upper(meta.ResearchGroup),varargin);
if nnz(idx_group)==0 % if no group specified, return all
    idx_group = true(size(idx_group));
end

% sex
idx_sex = ismember(upper(meta.Sex),varargin);
if nnz(idx_sex)==0 % if no group specified, return all
    idx_sex = true(size(idx_sex));
elseif nnz(idx_sex) < length(idx_sex)
    warning('APPLYING SEX LIMITATION');
end

% age
idx_age = meta.Age >= overval & meta.Age <= underval;
if nnz(idx_age)
    warning('APPLYING AGE RANGE LIMITATION: [%d-%d]',overval,underval);
end



idx = idx_group & idx_sex & idx_age;
t = meta(idx,:);
subjects = t.SubjectPath;
age = t.Age;
sex = t.Sex;
group = t.ResearchGroup;


