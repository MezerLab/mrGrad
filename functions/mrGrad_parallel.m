%% mrGrad: MRI Region Gardients
function RG = mrGrad(Data,varargin)
%--------------------------------------------------------------------------
% INPUTS:
%--------------------------------------------------------------------------
%
%   DATA    1xN-groups cell containing structs with fields:
%           * 'map_list': (cell) paths to subjects qMRI images (nii).
%           * 'seg_list': (cell) paths to subjects segmentation files or
%               binary mask files *in the same resolution* (nii).
%           Optional:
%           * 'group_name': specify subject group name (e.g. 'Older adults')
%           * 'subject_names': a list of subject names
%           * 'age': a list of subjects age
%           * 'sex': a list of subjects sex
%           The optional fields will be copied to the output structs for
%           documentaion purposes.
%           
%
%   'ROI':  followed by FreeSurfer/FSL labels (e.g. [11 50 12 51] for
%           l-caudate, r-caudate, l-putamen r-putamen);
%           in case 'seg_list' contains paths to binary masks, ROI shoud be 1
%           (or as labeld in the file). In case the provided labels do not refer 
%           to freesurfer's look-up table, please provide the 'roi_names'
%           argument.
%
%   OPTIONAL Name-Value Arguments
%
%   'Nsegs':    followed by a scalar or a 1x3 vector specifying number of
%               segments along axes. default: 7 
%
%   'segmentingMethod': followed by 'spacing' (default) or 'VoxN', to
%                    specify whether to use equal spaced segments or equal
%                    voxel count in segments.
%
%   'stat':   followed by the wanted statistic name for the qMRI function:
%             'median' (default) / 'mean'
%
%   'PC'        followed by PC number(s): 1,2,3 (default: all)
%
%   'max_change': a nROIs x nPCs array, contianing image axis numbers (1=X;
%                 2=Y; 3=Z).
%                 for multiple subjects we want the directionality (sign)
%                 of the individual data-driven ROI axes to be consistent
%                 (e.g. A >> P and not P >> A). 'max_change' is specifying
%                 the image axis for which the ROI axis has consistent
%                 change across subjects. for example: For putamen and
%                 caudate max_change is set by default to [2 3 1] (denoting
%                 [Y Z X]), in order to achieve consistent directionality
%                 of axis 1 (AP) with image Y axis (A>>P), axis 2 (VD) with
%                 image Z axis (V>>D), and axis 3 (ML) with image X axis (M>>L).
%                 
%   'erode':    followed by 0 (default) or 1 - remove outer surface of ROI
%               to avoid partial voluming
%
%   'invert':   followed by 0 (default) or 1 to invert data (e.g. T1 >> 1/T1)
%
%   'normalize': followed by 1 or 0 (default). Substract individual
%                baselines from gradients. basline = median value in the
%                ROI of the individual subject. This intends to eliminate
%                absolute differences between subjects, if wanted.
%
%   'roi_names': followed by a cell array with ROI names strings. This is
%                not needed if ROI labels refer to freesurfer's look-up table
%
%   'param': MRI parameter name (e.g. 'R1')
%
%   'units': MRI parameter units (e.g. 'sec^{-1}')
%
%   SOFTWARE REQUIREMENTS:
%
%        * MATLAB          - http://www.mathworks.com/products/matlab/
%        * Vistasoft       - https://github.com/vistalab/vistasoft   
%        * boundedline-pkg - https://github.com/kakearney/boundedline-pkg (recommended)
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------
fprintf('mrGrad\n(C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021\n')

% obligatory input

if isa(Data,'struct')
    Data = {Data};
end
Ngroups = numel(Data);

% make sure obligatory input and files exist
for gg = 1:Ngroups
    obfields = {'map_list','seg_list'};
    if any(~isfield(Data{gg},obfields))
        error('DATA fields ''map_list'',''seg_list'' are obligatory');
    end
    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.map_list);
    if any(idx)
        error('one or all input image files not exist.')
    end
    idx = cellfun(@(x) ~exist(x,'file'),Data{gg}.seg_list);
    if any(idx)
        error('one or all input segmentation files not exist.')
    end
end
[found, ROI, varargin] = argParse(varargin, 'ROI');
if ~found
    error('please specify ROI label');
end
NROIs = numel(ROI);

% optional arguments

[found, stat, varargin] = argParse(varargin, 'stat');
if ~found; stat = 'median'; end

[found, erode, varargin] = argParse(varargin, 'erode');
if ~found; erode = 0; end

[found, PC, varargin] = argParse(varargin, 'PC');
if ~found; PC = 1:3; end

[found, Nsegs, varargin] = argParse(varargin, 'Nsegs');
if ~found; Nsegs = 7; end
if numel(Nsegs) < numel(PC)
    Nsegs = repmat(Nsegs,1,length(PC));
end
if numel(Nsegs) < numel(PC)
    error('mismatch between Nsegs and NPCs');
end

[found, invert, varargin] = argParse(varargin, 'invert');
if ~found; invert = false; end

[found, param, varargin] = argParse(varargin, 'param');
if ~found; param = 'unknown_parameter'; end

[found, units, varargin] = argParse(varargin, 'units');
if ~found; units = 'unknown_units'; end

[found, roi_names, varargin] = argParse(varargin, 'roi_names');
if ~found
    roi_names = ROI_name(ROI);
end

[found_m, max_change, varargin] = argParse(varargin, 'max_change');
if found_m
    if size(max_change,1)==1
        max_change = repmat(max_change,NROIs,1);
    end
    if size(max_change,1)~=NROIs
        error('MAX_CHANGE should have a row for each ROI');
    else
        msg = 'using user''s input priors for PCs directionality sign';
        disp(msg)
    end
end

% make sure functions of SPM not run over matlab's nanstd
tmp = which('nanstd');
if contains(tmp,'spm')
    error('SPM functions cause interference. Please remove SPM package from matlab''s path');
end

%--------------------------------------------------------------------------
% Keep unused arguments for upcoming functions
varforward = varargin;
%--------------------------------------------------------------------------
% LOOP OVER SUBJECT GROUPS AND ROIS
%--------------------------------------------------------------------------

RG = cell(Ngroups,NROIs);
j=0;

gcp();

for gg = 1:Ngroups
    % SUBJECT GROUP NAME
    if isfield(Data{gg},'group_name')
        group = Data{gg}.group_name;
    else
        group = sprintf('SubjectGroup_%d',gg);
    end
    
    maps = Data{gg}.map_list;
    segmentations = Data{gg}.seg_list;
    
    
    for rr = 1:NROIs
        roi = ROI(rr);
        j=j+1;
        fprintf('\n(%d/%d) %s %s %s\n',j,numel(RG),roi_names{rr},param,group);
        
        % Unify the sign of axes directions across subjects, according to prior
        % knowledge if exists, about the image axis that consistently changes with
        % the PC axis.
        if found_m
            maxchange = max_change(rr,:);
        else
            msg = ['using priors for PCs directionality sign in ROI ',num2str(roi)];
            if     ismember(roi,[10,49]) % thalamus
                maxchange = [2 3 3];
            elseif ismember(roi,[11,50]) % caudate
                maxchange = [2 3 1]; % [y z x]
            elseif ismember(roi,[12,51]) % putamen
                maxchange = [2 3 1]; % [y z x]
            elseif ismember(roi,[13,52]) % pallidum
                maxchange = [2 3 1];
            else
                maxchange = [2 3 2];
                msg=['no default directionslity specs. for ROI ',num2str(roi),'. Agreement between subjects might be compromised'];
            end
            disp(msg)
        end

        %--------------------------------------------------------------------------
        % RUN OVER MULTIPLE SUBJECTS' DATA
        %--------------------------------------------------------------------------
        Nsubs = length(maps);
        Allsubs_rg_data = cell(Nsubs,1);
        fprintf('Computing ROI axes and gradients for %d subject...',Nsubs)
        clearvars stridesWarnFlag
        
        parfor ii = 1:Nsubs
%             fprintf('%d\n',ii); % uncomment for debugging
            %----------------------------------------------------------------------
            % load subject's qMRI data
            %----------------------------------------------------------------------
            mask = ROImask(segmentations{ii},roi,erode);
            qmap = readFileNifti(maps{ii});
            strides = keep_strides(qmap);
            qmap = qmap.data;

            %----------------------------------------------------------------------
            % Outliers removal
            %----------------------------------------------------------------------

            % warn about outliers
            remove_Outliers = false;
            if remove_Outliers
                Outlier = isoutlier(qmap(mask>0));
                mask(mask>0) = ~Outlier;
                warning('%d outliers removed.',nnz(Outlier));
            end
            qmap = double(mask.* qmap);

            %----------------------------------------------------------------------
            % To obtain 1/param (e.g. T2w/T1w from T1w/T2w or R1 from T1)
            %----------------------------------------------------------------------
            if invert
                warning('inverting map to obtain 1/parameter');
                qmap(qmap~=0) = 1./qmap(qmap>0);
            end
            %----------------------------------------------------------------------
            % make sure data is in positive strides (L>R P>A I>S)
            %----------------------------------------------------------------------
            dims=1:3;
            dimsflip = dims(strides<0);
            for d=dimsflip
                qmap = flip(qmap,d);
                mask = flip(mask,d);
                if ii==1
                    warning('images of some/all subjects are flipped to match positive strides.')
                end
            end

            %----------------------------------------------------------------------
            % MAIN FUNCTION EXECUTION mrgrad_per_sub.m
            %----------------------------------------------------------------------
            % single subject mrgrads in (up to) 3 PCs
            singlsb_rgs = arrayfun(@(x,y)  mrgrad_per_sub(qmap,mask,'PC',x,'Nsegs', y,...
                'stat',stat,'maxchange',maxchange,'subID',ii,varforward{:}),PC,Nsegs);
            Allsubs_rg_data{ii} = singlsb_rgs;
            %-------------------------------------
        end

        y = cell(1,length(PC));
        for jj = 1:length(y)
            y{jj} = [];
        end
        for ii=1:Nsubs
            y = cellfun(@(a,b) double([a b]), y, {Allsubs_rg_data{ii}.function}, 'UniformOutput',false);
        end
        %% Packing output
        rg.Y = y;
        rg.Y_mean = cellfun(@(x) nanmean(x,2), rg.Y, 'UniformOutput', false);
        rg.Y_std  = cellfun(@(x) nanstd(x,0,2), rg.Y, 'UniformOutput', false);
        rg.Y_SEM  = cellfun(@(x) nanstd(x,0,2)/sqrt(size(x,2)), rg.Y, 'UniformOutput', false);
        rg.X      = arrayfun(@(x) (1:x)', Nsegs, 'UniformOutput', false);
        rg.N_segments = Nsegs;
        rg.parameter = param;
        rg.units = units;
        rg.method = stat;
        rg.y_lbls = {'axis1','axis2','axis3'};
        rg.ROI_label = roi_names{rr};
        rg.individual_data = Allsubs_rg_data;
        
        for v = {'group_name','subject_names','age','sex'}
            if isfield(Data{gg},v)
                rg.(v{:}) = Data{gg}.(v{:});
            end
        end

        %------------------------
        % Flip PA to AP, LM to ML
        %------------------------
        % flip PA to AP in all striata (11,12,50,51)
        ax = 1;
        if maxchange(ax)==2 % y-coordinate
            rg = RG_flip(rg, ax);
        end

        % flip LM to ML in left-hemisphere striata (11,12)
        ax = 3;
        if ismember(roi,[11,12]) && maxchange(ax)==1 % x-coordinate
            rg = RG_flip(rg, ax);
        end
        RG{gg,rr} = rg;
        fprintf(2,' done!\n');
    end
end
% delete(gcp);
fprintf('\nAll done!\n');
end
% function strides = keep_strides(nifti_struct)
%     strides = diag(nifti_struct.qto_xyz(1:3,1:3));
%     strides = strides./abs(strides);
% 
% end