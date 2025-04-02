%% mrGrad: MRI Region Gardients
function [RG, T] = mrGrad(Data,varargin)
%--------------------------------------------------------------------------
% INPUTS:
%--------------------------------------------------------------------------
%
%   DATA    1xN-groups cell containing structs with fields:
%           * 'map_list': (cell) paths to subjects qMRI images (nii).
%           * 'seg_list': (cell) paths to subjects segmentation files or
%               binary mask files *in the same resolution* (nii).
%           
%           User can provide the data separated into research groups or
%           provide the whole cohort as one group.
%
%           Optional:
%           * 'group_name': specify subject group name (e.g. 'Older adults')
%           * 'subject_names': a list of subject names (ideally unique names)
%           * 'age': a list of subjects age
%           * 'sex': a list of subjects sex
%           * and any other descriptive fields
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
%   'apply_alternative_axes': If users choose, axes can be computed on an
%               alternative ROI A, then be applied for sampling the ROI of
%               analysis B. In this case, users should provide a struct
%               argument with field 'seg_list' containing the alternative
%               segmentation paths (or a 1xN-groups cell) (can be the same as the main seg_list),
%               and a field 'ROI' containing the alternative ROI labels.
%               NOTE that this option should be used carefully as the used
%               axes would not carry anatomical meaning in the sampled ROI.
%               EXAMPLE: mrGrad(...'apply_alternative_axes',AlternativeAxes)
%                           where:
%                           AlternativeAxes.seg_list{1} = {paths/subject_group_1}
%                           AlternativeAxes.seg_list{2} = {paths/subject_group_2}
%                           AlternativeAxes.ROI = [5];
%
%   'ignore_missing': boolean. If true, ignore subjects with missing data
%                     and run anyway. (Default: false)
%
%   'outfile': full path to save output
%
%   'output_mode':  'minimal' / 'default' / 'extended'
%                   'minimal' mode saves only the summarized parameter
%                   value results; 'default' mode saves also the individual
%                   subjects' axis segmentation information; 'extended'
%                   mode also saves the new axis segmentation masks per
%                   subject.
%
%   SOFTWARE REQUIREMENTS:
%
%        * MATLAB          - http://www.mathworks.com/products/matlab/
%        * boundedline-pkg - https://github.com/kakearney/boundedline-pkg (recommended)
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------
fprintf('mrGrad\n(C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021\n')
mrgrad_defs = setGlobalmrgrad(varargin{:});
mrgrad_defs.fname = mfilename;

Parallel = isequal(mrgrad_defs.fname,'mrGrad_parallel');
if ~Parallel
    fprintf(2,' mrGrad() may run slowly — consider using mrGrad_parallel() for faster performance.');
end


% check obligatory input
Data = mrgrad_check_input(Data,mrgrad_defs);


Ngroups = numel(Data);
NROIs = numel(mrgrad_defs.ROI);

% % make sure functions of SPM not run over matlab's nanstd
% nanstd_path = which('nanstd');
% if contains(nanstd_path,'spm')
%     error('SPM functions cause interference. Please remove SPM package from matlab''s path');
% end

%--------------------------------------------------------------------------
% LOOP OVER SUBJECT GROUPS AND ROIS
%--------------------------------------------------------------------------
warning('on','mrGrad:Strides');
RG = cell(Ngroups,NROIs);
j=0;

for gg = 1:Ngroups
    maps = Data{gg}.map_list;
    segmentations = Data{gg}.seg_list;
    group = Data{gg}.group_name;
    Nsubs = length(Data{gg}.map_list);

    for rr = 1:NROIs
        roi = mrgrad_defs.ROI(rr);
        roi_name = mrgrad_defs.roi_names{rr};
        j=j+1;
        fprintf('\n(%d/%d) %s %s %s\n',j,numel(RG),mrgrad_defs.roi_names{rr},mrgrad_defs.param,group);
        
        % Unify the sign of axes directions across subjects, according to prior
        % knowledge if exists, about the image axis that consistently changes with
        % the PC axis (e.g., all subjects spatial functions will be A>>P and not P>>A).
        if ~isempty(mrgrad_defs.max_change)
            maxchange_roi = mrgrad_defs.max_change(rr,:);
        elseif isfield(mrgrad_defs,'Alternative_ROI')
            maxchange_roi = [2 3 2];
            msg = ['no default directionslity specs. for Alternative ROI ',num2str(mrgrad_defs.Alternative_ROI(rr)),'. Agreement between subjects might be compromised'];
            disp(msg);
        elseif ismember(roi_name,ROI_name(roi))
            [maxchange_roi,msg] = get_roi_priors(roi);
            disp(msg);
        else
            maxchange_roi = [2 3 2];
            msg = ['no default directionslity specs. for ROI ',num2str(roi),'. Agreement between subjects might be compromised'];
            disp(msg);            
        end
        
        %--------------------------------------------------------------------------
        % RUN OVER MULTIPLE SUBJECTS' DATA
        %--------------------------------------------------------------------------
        Allsubs_rg_data = cell(Nsubs,1);
        fprintf('Computing ROI axes and gradients for %d subject...',Nsubs)
        
        
        stat = mrgrad_defs.stat;
        PC = mrgrad_defs.PC;
        Nsegs = mrgrad_defs.Nsegs;
        sampling_method = mrgrad_defs.segmentingMethod;
        BL_normalize = mrgrad_defs.BL_normalize;
        isfigs = mrgrad_defs.isfigs;

        
        for ii = 1:Nsubs
%             fprintf('%d\n',ii); % uncomment for debugging
            %----------------------------------------------------------------------
            % load subject's qMRI data
            %----------------------------------------------------------------------
            
            if any(cellfun(@(x) ~exist(x,"file"),[maps(ii);segmentations(ii)]))
                continue
            end

            mask = ROImask(segmentations{ii},roi,mrgrad_defs.erode_flag);

            image_info = niftiinfo(maps{ii});
            [strides,im_dims] = keep_strides(image_info);

            im = niftiread(maps{ii});
            im = single(im);
            imsize = size(im);

            % Make sure mask and image have the same dimensions
            if ~isequal(size(mask),imsize)
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
            if mrgrad_defs.invert_flag
                warning('Inverting map to obtain 1/parameter');
                im(im~=0) = 1./im(im>0);
            end
            %----------------------------------------------------------------------
            % if an alternative ROI was given for Axes calculation, make mask
            %----------------------------------------------------------------------
            alternative_mask = [];
            if isfield(mrgrad_defs,'Alternative_ROI')
                alternative_roi = mrgrad_defs.Alternative_ROI(rr);
                alternative_seg = mrgrad_defs.Alternative_seg_list{gg}{ii};
                alternative_mask = ROImask(alternative_seg,alternative_roi,0);
            end
            %----------------------------------------------------------------------
            % make sure data is in positive strides (L>R P>A I>S)
            %----------------------------------------------------------------------

            % change image strides order to [1,2,3]
            [~, im_perm] = sort(im_dims);
            im = permute(im,im_perm);
            mask = permute(mask,im_perm);
            alternative_mask = permute(alternative_mask,im_perm);

            % flip negative strides to achieve [+1,+2,+3]
            dimsflip = im_dims(strides < 0);
            for d = dimsflip
                im = flip(im,d);
                mask = flip(mask,d);
                alternative_mask = flip(alternative_mask,d);
            end
            % Warn once about strides' change
            if ~isequal(strides,[1,2,3]) && ...
                    (ii==1 || ~Parallel)
                warning('mrGrad:Strides','\nImages of some/all subjects are flipped to match positive strides.')
                warning('off','mrGrad:Strides');
            end

            %----------------------------------------------------------------------
            % MAIN FUNCTION CALL mrgrad_per_sub.m
            %----------------------------------------------------------------------
            % single subject mrgrads in (up to) 3 PCs
            singlsb_rgs = arrayfun(@(x,y)...
                mrgrad_per_sub(im,mask,'PC',x,'Nsegs',y,'sampling_method',sampling_method,...
                'stat',stat,'maxchange',maxchange_roi,'BL_normalize',BL_normalize,...
                'subID',ii,'isfigs',isfigs,'apply_alternative_axes',alternative_mask),...
                PC,Nsegs);

            % Keep original stride info for generating segmentation files
            [singlsb_rgs.analysis_image_size] = deal(size(im));
            [singlsb_rgs.original_image_size] = deal(imsize);
            singlsb_rgs(1).original_nifti_info = image_info;
            
            [singlsb_rgs.analysis_strides] = deal([1,2,3]);
            [singlsb_rgs.original_strides] = deal(strides);

            Allsubs_rg_data{ii} = singlsb_rgs;
            %-------------------------------------
        end

        % Combine gradient data from all subjects (allowing missing data)
        y = arrayfun(@(n_segs) nan(n_segs,Nsubs),Nsegs,'un',0);
        for ii = 1:Nsubs
            if isfield(Allsubs_rg_data{ii},"function")
                for jj = 1:length(Allsubs_rg_data{ii})
                    y{jj}(:,ii) = Allsubs_rg_data{ii}(jj).function;
                end
            end
        end

        %% Packing output  
        rg.Y = y;
        rg.Y_mean = cellfun(@(x) mean(x,2,"omitnan"),rg.Y,'un',0);
        rg.Y_std  = cellfun(@(x) std(x,0,2,"omitnan"),rg.Y,'un',0);
        rg.Y_SEM  = cellfun(@(x) std(x,0,2,"omitnan")/sqrt(size(x,2)),rg.Y,'un',0);
        rg.X      = arrayfun(@(x) (1:x)',mrgrad_defs.Nsegs,'un',0);
        rg.N_segments = mrgrad_defs.Nsegs;
        rg.parameter = mrgrad_defs.param;
        rg.units = mrgrad_defs.units;
        rg.sampling_method = mrgrad_defs.segmentingMethod;
        rg.method = mrgrad_defs.stat;
        rg.y_lbls = cellstr("axis" + PC);
        rg.ROI_label = mrgrad_defs.roi_names{rr};
        if ~strcmpi(mrgrad_defs.output_mode,"minimal")
            rg.individual_data = Allsubs_rg_data;
        end
        description_fields = fieldnames(Data{gg})';
        for v = description_fields
            rg.(v{:}) = Data{gg}.(v{:});
        end

        %------------------------
        % Flip PA to AP, LM to ML
        %------------------------
        % flip PA to AP in all striata (11,12,50,51)
        ax = 1; % AP
        if maxchange_roi(ax)==2 ... % y-coordinate
            && ismember(ax,PC) % only if pc1 is included in analysis
            rg = RG_flip(rg, ax);
        end

        % flip LM to ML in left-hemisphere striata (11,12)
        ax = 3;
        if ismember(roi,[11,12]) && maxchange_roi(ax)==1 ... % x-coordinate
            && ismember(ax,PC) % only if pc3 is included in analysis
            rg = RG_flip(rg, find(ismember(PC,ax)));
        end
        RG{gg,rr} = rg;
        fprintf(2,' done!\n');
    end
end

% Generate Summary Results
T = mrgrad_rg2table(RG);

if ~isempty(mrgrad_defs.outfile)

    fprintf('\nSaving Summary Outputs to disc... ');
    % save .mat file
    out_mat = string(mrgrad_defs.outfile);
    save(out_mat,'RG',"-mat");

    % save .csv file
    out_csv = regexprep(out_mat,".mat$",".csv");
    writetable(T,out_csv);

    % make sure files saved
    Saved = all(arrayfun(@(x) exist(x,"file"), [out_mat,out_csv]));
    if Saved
        fprintf(' done!\n');
    else
        fprintf(2,'\nAn error occurred while saving the output files. The results were not saved!\n');
    end

    % extended mode: Generate segmentation masks
    if mrgrad_defs.output_mode == "extended"
        fprintf('\nExtended output mode: saving result segmentations to disc... ');
        seg_output_dir = fileparts(mrgrad_defs.outfile);
        if ~Parallel
            fprintf('\n This may take a while. Consider using mrGrad_parallel() for faster performance.\n' )
        end
        for jj = 1:numel(RG)
            output_dir = mrGrad_seg(RG{jj},seg_output_dir,true,Parallel);
        end
        if ~isempty(dir(fullfile(output_dir,"**/*.nii.gz")))
            fprintf(' done!\n Segmentation files were saved in: %s\n',output_dir);
        else
            fprintf(2,'\nAn error occurred while saving the output files. The results were not saved!\n');
        end
    end

end

clear mrgrad_defs
fprintf('\nAll done!\n');
end