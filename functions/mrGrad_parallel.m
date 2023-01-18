%% mrGrad: MRI Region Gardients
function RG = mrGrad_parallel(Data,varargin)
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
%   'segmentingMethod': followed by 'equidistance' (default) or
%               'equivolume', to specify whether to use equally-spaced
%               segments or segemnts of equal voxel count.
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
clear -global mrgrad_defs
global mrgrad_defs
setGlobalmrgrad(varargin{:});
mrgrad_defs.fname = mfilename;

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


NROIs = numel(mrgrad_defs.ROI);


% make sure functions of SPM not run over matlab's nanstd
nanstd_path = which('nanstd');
if contains(nanstd_path,'spm')
    error('SPM functions cause interference. Please remove SPM package from matlab''s path');
end

%--------------------------------------------------------------------------
% LOOP OVER SUBJECT GROUPS AND ROIS
%--------------------------------------------------------------------------

RG = cell(Ngroups,NROIs);
j=0;

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
        roi = mrgrad_defs.ROI(rr);
        j=j+1;
        fprintf('\n(%d/%d) %s %s %s\n',j,numel(RG),mrgrad_defs.roi_names{rr},mrgrad_defs.param,group);
        
        % Unify the sign of axes directions across subjects, according to prior
        % knowledge if exists, about the image axis that consistently changes with
        % the PC axis (e.g., all subjects spatial functions will be A>>P and not P>>A).
        if ~isempty(mrgrad_defs.max_change)
            maxchange_roi = mrgrad_defs.max_change(rr,:);
        else
            [maxchange_roi,msg] = get_roi_priors(roi);
            disp(msg);
        end
        
        %--------------------------------------------------------------------------
        % RUN OVER MULTIPLE SUBJECTS' DATA
        %--------------------------------------------------------------------------
        Nsubs = length(maps);
        Allsubs_rg_data = cell(Nsubs,1);
        fprintf('Computing ROI axes and gradients for %d subject...',Nsubs)
        
        
        erode_flag = mrgrad_defs.erode_flag;
        invert_flag = mrgrad_defs.invert_flag;
        fname = mrgrad_defs.fname;
        stat = mrgrad_defs.stat;
        PC = mrgrad_defs.PC;
        Nsegs = mrgrad_defs.Nsegs;
        sampling_method = mrgrad_defs.segmentingMethod;
        BL_normalize = mrgrad_defs.BL_normalize;
        isfigs = mrgrad_defs.isfigs;

        parfor ii = 1:Nsubs
%             fprintf('%d\n',ii); % uncomment for debugging
            %----------------------------------------------------------------------
            % load subject's qMRI data
            %----------------------------------------------------------------------
            mask = ROImask(segmentations{ii},roi,erode_flag);
            im = readFileNifti(maps{ii});
            strides = keep_strides(im);
            im = im.data;

            % Make sure mask and image have the same dimensions
            if ~isequal(size(mask),size(im))
                error('Input image and mask/segmentation'' dimensions must agree.');
            end
            %----------------------------------------------------------------------
            % Outliers removal
            %----------------------------------------------------------------------
            % warn about outliers
            remove_Outliers = false;
            if remove_Outliers
                Outlier = isoutlier(im(mask>0));
                mask(mask>0) = ~Outlier;
                warning('%d outliers removed.',nnz(Outlier));
            end
            %----------------------------------------------------------------------
            % mask the image
            %----------------------------------------------------------------------
            im = double(mask.* im);

            %----------------------------------------------------------------------
            % To obtain 1/param (e.g. T2w/T1w from T1w/T2w or R1 from T1)
            %----------------------------------------------------------------------
            if invert_flag
                warning('inverting map to obtain 1/parameter');
                im(im~=0) = 1./im(im>0);
            end
            %----------------------------------------------------------------------
            % make sure data is in positive strides (L>R P>A I>S)
            %----------------------------------------------------------------------
            dims = 1:3;
            dimsflip = dims(strides<0);
            for d = dimsflip
                im = flip(im,d);
                mask = flip(mask,d);
                
                % Warn once about strides' change
                if isequal(fname,'mrGrad') || ii==1
                    fprintf('\n');
                    warning('mrGrad:Strides','images of some/all subjects'' images are flipped to match positive strides.')
                    warning('off','mrGrad:Strides');
                end
            end
            %----------------------------------------------------------------------
            % MAIN FUNCTION EXECUTION mrgrad_per_sub.m
            %----------------------------------------------------------------------
            % single subject mrgrads in (up to) 3 PCs
            singlsb_rgs = arrayfun(@(x,y)...
                mrgrad_per_sub(im,mask,'PC',x,'Nsegs',y,'sampling_method',sampling_method,...
                'stat',stat,'maxchange',maxchange_roi,'BL_normalize',BL_normalize,...
                'subID',ii,'isfigs',isfigs),...
                PC,Nsegs);
            Allsubs_rg_data{ii} = singlsb_rgs;
            %-------------------------------------
        end

        y = cell(1,length(mrgrad_defs.PC));
        for jj = 1:length(y)
            y{jj} = [];
        end
        for ii=1:Nsubs
            y = cellfun(@(a,b) double([a b]), y, {Allsubs_rg_data{ii}.function},'un',0);
        end
        %% Packing output
        
        rg.Y = y;
        rg.Y_mean = cellfun(@(x) nanmean(x,2),rg.Y,'un',0);
        rg.Y_std  = cellfun(@(x) nanstd(x,0,2),rg.Y,'un',0);
        rg.Y_SEM  = cellfun(@(x) nanstd(x,0,2)/sqrt(size(x,2)),rg.Y,'un',0);
        rg.X      = arrayfun(@(x) (1:x)',mrgrad_defs.Nsegs,'un',0);
        rg.N_segments = mrgrad_defs.Nsegs;
        rg.parameter = mrgrad_defs.param;
        rg.units = mrgrad_defs.units;
        rg.sampling_method = mrgrad_defs.segmentingMethod;
        rg.method = mrgrad_defs.stat;
        rg.y_lbls = {'axis1','axis2','axis3'};
        rg.ROI_label = mrgrad_defs.roi_names{rr};
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
        ax = 1; % AP
        if maxchange_roi(ax)==2 % y-coordinate
            rg = RG_flip(rg, ax);
        end

        % flip LM to ML in left-hemisphere striata (11,12)
        ax = 3;
        if ismember(roi,[11,12]) && maxchange_roi(ax)==1 % x-coordinate
            rg = RG_flip(rg, ax);
        end
        RG{gg,rr} = rg;
        fprintf(2,' done!\n');
    end
end
clear -global mrgrad_defs
fprintf('\nAll done!\n');
end
function strides = keep_strides(nifti_struct)
    strides = diag(nifti_struct.qto_xyz(1:3,1:3));
    strides = strides./abs(strides);

end