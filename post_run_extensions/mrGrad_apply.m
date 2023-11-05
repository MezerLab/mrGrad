function RG = mrGrad_apply(RG_path,Data,varargin)
% This function applies the mrGrad axis computations on a new set of maps.
% i.e., if you already run mrGrad on a set of maps (e.g. R1), you can use
% this function to compute spatial profiles of the same (or a subset) of
% subjects in another map (e.g., R2*) without the need to compute the axes
% again using mrGrad.
% NOTE:
% (1) The two maps must be in the same space
% (2) Both the original mrGrad struct and the DATA input must have the
% field subject_names, and the latter must be a subset of the former.
%--------------------------------------------------------------------------
% INPUTS:
%--------------------------------------------------------------------------
%
%   RG_path: the path to the original mrGrad cell array
%
%   Data: 1xN-groups cell containing structs with fields:
%           * 'map_list': (cell) paths to subjects qMRI images (nii).
%           * 'subject_names': a list of subject names
%--------------------------------------------------------------------------
%   OPTIONAL Name-Value Arguments
%
%   'stat':   followed by the wanted statistic name for the qMRI function:
%             'median' (default) / 'mean'
%
%   'param': MRI parameter name (e.g. 'R1')
%
%   'units': MRI parameter units (e.g. 'sec^{-1}')
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

[found, stat, varargin] = argParse(varargin, 'stat');
if ~found; stat = 'median'; end
[found, param, varargin] = argParse(varargin, 'param');
if ~found; param = 'unknown_parameter'; end
[found, units, varargin] = argParse(varargin, 'units');
if ~found; units = 'unknown_units'; end

RG = load(RG_path);
RG = RG.RG;

Ngroups = size(RG,1);
Nrois = size(RG,2);

if ~isequal(Ngroups, length(Data))
    error('subject groups do not match');
end

%--------------------------------------------------------------------------
for gg = 1:Ngroups
    for rr = 1:Nrois
        maps = Data{gg}.map_list;
        subject_names = Data{gg}.subject_names;
        Nsubjects = length(subject_names);
        rg = RG{gg,rr};
        Naxes = length(rg.Y);
        
        rg.Y = {};
        rg.Y_mean = {};
        rg.Y_std  = {};
        rg.Y_SEM  = {};
        
        if ~isfield(rg,'subject_names')
            error('original mrGrad structs have no subject_names fields. comparison is not possible');
        end
        subjects_orig = rg.subject_names;

        % make sure all new subjects appear in old RG
        idx = ~ismember(subject_names,subjects_orig);
        if any(idx)
            error('axes were not previously calculated for one or more new subjects');
        end
        
        %------------------------------------------------------------------
        % remove original subjects that don't have new maps
        %------------------------------------------------------------------
        idx_rm = ~ismember(subjects_orig,subject_names);
        rg.individual_data(idx_rm) = [];
        rg.subject_names(idx_rm) = [];
        
        if isfield(rg,'age')
            rg.age(idx_rm) = [];
        end
        if isfield(rg,'sex')
            rg.sex(idx_rm) = [];
        end
        
        %------------------------------------------------------------------
        % loop over subjects and compute new image statistics (profiles)
        %------------------------------------------------------------------
        
        for ii = 1:Nsubjects
            
            sub = subject_names{ii};
            sub_index = find(ismember(rg.subject_names,{sub}));

            qmap = readFileNifti(maps{ii});
            strides = keep_strides(qmap);
            qmap = qmap.data;
            sz = size(qmap);
            
            roi_inds = rg.individual_data{sub_index}.all_Inds;
            roi_inds = sub2ind(sz,roi_inds(:,1),roi_inds(:,2),roi_inds(:,3));
            tmp = zeros(size(qmap));
            tmp(roi_inds) = qmap(roi_inds);
            qmap = tmp;
            %----------------------------------------------------------------------
            % make sure data is in positive strides (L>R P>A I>S)
            %----------------------------------------------------------------------
            dims=1:3;
            dimsflip = dims(strides<0);
            for d=dimsflip
                qmap = flip(qmap,d);
            end
            %----------------------------------------------------------------------
            % compute new map's stats
            %----------------------------------------------------------------------
            for ax = 1:Naxes
                Inds = rg.individual_data{sub_index}(ax).segment_inds_linear;
                if strcmpi(stat,'mean')
                    func = cellfun(@(x) mean(qmap(x)), Inds);
                    sub_err = cellfun(@(x) std(qmap(x)), Inds);
                    error_type = 'std';
                elseif strcmpi(stat,'median')    
                    func = cellfun(@(x) median(qmap(x)), Inds);
                    sub_err = cellfun(@(x) mad(qmap(x)), Inds);
                    error_type = 'mad';
                end
                
                % Change the individual's data that is dependent on image
                % values (fields: function, error,err_type). The rest of
                % the fields are independent of image values.
                rg.individual_data{sub_index}(ax).function = func;
                rg.individual_data{sub_index}(ax).error = sub_err;
                rg.individual_data{sub_index}(ax).err_type = error_type;
            end
            fprintf('%d/%d\n',ii,Nsubjects);
        end
        %----------------------------------------------------------------------
        % save spatial data for all subjects in one table y
        %----------------------------------------------------------------------
        y = cell(1,Naxes);
        for ii = 1:Nsubjects
            y = cellfun(@(a,b) double([a b]), y, {rg.individual_data{ii}.function}, 'un',0);
        end
        %----------------------------------------------------------------------
        % update the group-level stats
        %----------------------------------------------------------------------
        rg.Y = y;
        rg.Y_mean = cellfun(@(x) nanmean(x,2), rg.Y, 'un', 0);
        rg.Y_std  = cellfun(@(x) nanstd(x,0,2), rg.Y, 'un', 0);
        rg.Y_SEM  = cellfun(@(x) nanstd(x,0,2)/sqrt(size(x,2)), rg.Y, 'un', 0);
        rg.parameter = param;
        rg.units = units;
        rg.method = stat;
        RG{gg,rr} = rg;
    end
end

