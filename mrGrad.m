function [RG, T] = mrGrad(Data,varargin)
%% mrGrad v2.0 — MRI Regional Gradients
%
%--------------------------------------------------------------------------
% DESCRIPTION:
%--------------------------------------------------------------------------
%   Computes spatial gradients of MRI parameters within specified brain
%   subcortical regions (ROIs). The function supports multiple MRI
%   modalities per subject and provides summary statistics across participants.
%
%--------------------------------------------------------------------------
% INPUTS:
%--------------------------------------------------------------------------
%   Data : struct containing subject-level data and metadata.
%       Required fields:
%
%       - seg_list : {N x 1} cell array
%             Paths to segmentation or binary mask files for each subject.
%
%       - map_list : {N x M} cell array
%             Paths to MRI image files (e.g., T1, T2, R1) corresponding to
%             each subject. Multiple modalities (M) are supported, provided
%             all are coregistered with seg_list.
%
%       Optional fields:
%       - subject_ids : {N x 1} cell array or numeric vector
%             Unique identifiers per subject (e.g., {'sub-001', 'sub-002'}).
%
%       - Any other fields with the same length as seg_list/map_list
%             (e.g., 'Age', 'Sex', 'Group'). These fields are retained
%             in the output for reference.
%
%--------------------------------------------------------------------------
% OPTIONAL NAME–VALUE PAIRS:
%--------------------------------------------------------------------------
%
%   'ROI' : vector
%       FreeSurfer/FSL label indices (e.g., [11 50 12 51]).
%       If nonstandard labels are used, specify 'roi_names'.
%
%   'roi_names' : cellstr
%       ROI names corresponding to each label.
%
%   'Axes' or 'PC' : vector
%       ROI principal component indices to compute (e.g., [1 2 3]).
%       Default: all.
%
%   'n_segments' : scalar or vector
%       Number of bins (segments) per axis. Default: 7.
%
%   'segmenting_method' : char
%       'equidistance' (default) or 'equivolume'.
%
%   'stat' : char
%       Statistic used per segment ('median' [default] or 'mean').
%
%   'max_change' : numeric array [nROIs x nPCs]
%       Specifies the dominant anatomical axis (1 = X, 2 = Y, 3 = Z)
%       aligned with each principal component (PC) for every ROI.  
%       Ensures consistent **directionality** of PC axes across subjects,
%       preventing sign flips (e.g., anterior -> posterior vs. posterior -> anterior).
%
%       Example:
%           max_change = [2 3 1];
%           -> PC1 aligns with Y-axis (A -> P)
%           -> PC2 aligns with Z-axis (V -> D)
%           -> PC3 aligns with X-axis (M -> L)
%
%       If multiple ROIs are analyzed, use an [nROIs x nPCs] matrix.
%       If omitted, mrGrad attempts to infer defaults or reverts to [2 3 2],
%       printing a warning if subject-level alignment may differ.
%
%       Note:
%           Only affects axis sign (orientation), not gradient magnitude.
%           Recommended when combining hemispheres or enforcing anatomical
%           consistency across datasets.
%
%   'erode' : logical
%       Remove outer voxel layer to minimize partial-volume effects.
%       Default: false.
%
%   'parameter_names' : string/cell
%       Names of MRI parameters (e.g., {'R1', 'MTsat'}).
%
%   'units' : string/cell
%       Units of corresponding parameters (e.g., {'1/s', 'p.u.'}).
%
%   'apply_alternative_axes' : struct
%       Defines alternative axes based on a separate ROI segmentation.
%       Fields:
%           seg_list : cell array of alternative segmentations.
%           ROI      : label(s) for axis derivation.
%       Note: Use cautiously - axes may not retain anatomical meaning.
%
%   'allow_missing' : logical
%       Skip missing subjects instead of throwing an error. Default: false.
%
%   'output_name' : char
%       Base name for saved summary files (.mat and .csv).
%
%   'output_dir' : char
%       Directory for saving results.
%
%   'output_mode' : char
%       'minimal'  - summary tables only (.mat, .csv)
%       'default'  - includes per-subject data
%       'extended' - also saves per-subject segmentation masks (.nii.gz)
%
%   'Parallel' : logical
%       Use MATLAB Parallel Toolbox (parfor). Recommended for large cohorts.
%
%--------------------------------------------------------------------------
% OUTPUTS:
%--------------------------------------------------------------------------
%   RG : struct [1 x nROIs]
%       ROI-specific results containing:
%           * Gradient values per segment and axis
%           * Mean/SD/SEM statistics
%           * Subject-level data
%           * Metadata (parameter names, ROI labels)
%
%   T : table
%       Summary table combining all ROI and subject data.
%
%--------------------------------------------------------------------------
% DEPENDENCIES:
%--------------------------------------------------------------------------
%   MATLAB (https://www.mathworks.com)
%   boundedline-pkg (recommended)
%       https://github.com/kakearney/boundedline-pkg
%
%--------------------------------------------------------------------------
% AUTHORS:
%   Elior Drori, Aviv A. Mezer
%   Mezer Lab, The Hebrew University of Jerusalem
%   Copyright (C) 2021
%
%--------------------------------------------------------------------------

fprintf('\nmrGrad Toolbox (v2.0)\n(C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021\n')
mrgrad_defs = mrgrad_set(varargin{:});
mrgrad_defs.fname = mfilename;

% check obligatory input
[Data, mrgrad_defs] = mrgrad_check_input(Data,mrgrad_defs);
num_regions = numel(mrgrad_defs.ROI);
num_subjects = length(Data.seg_list);
parameter_str = strjoin(strrep(unique(cellstr(mrgrad_defs.parameter_names)),'unknown_parameter','unknown parameter'),", ");

if mrgrad_defs.parallel
    gcp();
else
    fprintf(2,' May run slowly on large cohorts — for faster performance set Parallel=true.');
end
%--------------------------------------------------------------------------
% LOOP OVER SUBJECT GROUPS AND ROIS
%--------------------------------------------------------------------------
warning('on','mrGrad:Strides');
RG = struct;
invalidPattern = '[\\/:*?"<>|-]|[\x00-\x1F]';
rg_names = regexprep(mrgrad_defs.roi_names,invalidPattern,'_');
for rr = 1:num_regions
    roi = mrgrad_defs.ROI(rr);
    roi_name = mrgrad_defs.roi_names{rr};
    fprintf('\n(%d/%d) %s %s\n',rr,num_regions,mrgrad_defs.roi_names{rr},parameter_str);

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
    fprintf('Computing ROI axes and parameter gradients for %d subjects...', num_subjects)

    Allsubs_rg_data = cell(num_subjects,1);
    seg_list = Data.seg_list;
    map_list = Data.map_list;
    if mrgrad_defs.parallel && num_subjects > 1
        gcp();
        parfor ii = 1:num_subjects
            Allsubs_rg_data{ii} = mrgrad_inner_loop(ii, rr, seg_list, map_list, mrgrad_defs, maxchange_roi);
        end
    else
        for ii = 1:num_subjects
            Allsubs_rg_data{ii} = mrgrad_inner_loop(ii, rr, seg_list, map_list, mrgrad_defs, maxchange_roi);
        end
    end

    %--------------------------------------------------------------------------
    % SUMMARIZE ALL SUBJECTS' RESULTS
    %--------------------------------------------------------------------------

    PC = mrgrad_defs.PC;
    n_segments = mrgrad_defs.n_segments;

    % allow integration of subjects with no data
    idx_no_data = cellfun(@isempty, Allsubs_rg_data);
    if any(idx_no_data)
        pseudo_data = struct;
        for pc = 1:length(PC)
            pseudo_data(pc).function = single(nan(n_segments(pc),length(mrgrad_defs.parameter_names)));
        end
        Allsubs_rg_data(idx_no_data) = repmat({pseudo_data},nnz(idx_no_data),1);
    end

    rg = struct;
    rg.roi_name = mrgrad_defs.roi_names{rr}; % change to region_name
    rg.axis_names = cellstr("axis" + PC); % change to axis_names
    rg.n_segments = mrgrad_defs.n_segments;
    rg.sampling_method = mrgrad_defs.segmentingMethod;
    rg.parameter_names = cellstr(mrgrad_defs.parameter_names(:)');
    rg.parameter_units = cellstr(mrgrad_defs.units(:)');
    rg.summary_stat = mrgrad_defs.stat;
    rg.subject_ids = Data.subject_ids;

    % summarize data in tables
    rg.Results = struct;
    summary_rows = rg.subject_ids;
    for kk = 1:length(rg.axis_names)
        summary_vars = "seg"+(1:n_segments(kk));
        for pp = 1:length(mrgrad_defs.parameter_names)
            summary_name = strjoin({rg.axis_names{kk},mrgrad_defs.parameter_names{pp}},"_");
            summary_mat = cellfun(@(x) x(kk).function(:,pp)',Allsubs_rg_data,'un',0);
            summary_mat = cat(1,summary_mat{:});

            rg.Results.(summary_name).x = (1:n_segments(kk))/n_segments(kk);
            rg.Results.(summary_name).y = ...
                array2table(summary_mat,VariableNames=summary_vars,...
                RowNames=summary_rows);


            rg.Results.(summary_name).y_mean = mean(rg.Results.(summary_name).y{:,:},1,"omitnan");
            rg.Results.(summary_name).y_std = std(rg.Results.(summary_name).y{:,:},[],1,"omitnan");
            n_actual = sum(~isnan(rg.Results.(summary_name).y{:,:}),1);
            rg.Results.(summary_name).y_sem = rg.Results.(summary_name).y_std ./ sqrt(n_actual);
        end
    end

    % keep individual subjects' roi-axis segmentation results
    Allsubs_rg_data(idx_no_data) = {[]};
    if ~strcmpi(mrgrad_defs.output_mode,"minimal")
        rg.individual_data = Allsubs_rg_data;
    end

    % keep user's descriptive fields (if valid)
    mrgrad_reserved_fields = fieldnames(rg)';
    description_fields = setdiff(fieldnames(Data)',mrgrad_reserved_fields,'stable');
    is_valid = cellfun(@(v) mrgrad_valid_descriptive_field(Data.(v)), description_fields);
    description_fields(~is_valid) = [];
    for v = description_fields
        user_data = Data.(v{:});
        rg.user_input_fields.(v{:}) = user_data;
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

    RG.(rg_names{rr}) = rg;
    fprintf(2,' done!\n');
end

T = mrgrad_rg2table(RG);

% Save Results
fprintf('\nSaving summary outputs to disc... ');
if ~exist(mrgrad_defs.output_dir,"dir")
    mkdir(mrgrad_defs.output_dir);
end

% save .mat file
out_mat = fullfile(mrgrad_defs.output_dir,mrgrad_defs.output_name);
save(out_mat,'RG',"-mat");

% save .csv file
out_csv = regexprep(out_mat,".mat$",".csv");
writetable(T,out_csv);

% make sure files were saved
Saved = all(arrayfun(@(x) exist(x,"file"), [out_mat,out_csv]));
if Saved
    fprintf(' done!\n Summary files were saved in: %s\n',mrgrad_defs.output_dir);
else
    fprintf(2,'\nAn error occurred while saving the output files. The results were not saved!\n');
end

% extended mode: Generate segmentation masks
if mrgrad_defs.output_mode == "extended"
    fprintf('\nExtended output mode: saving result segmentations to disc... ');
    if ~mrgrad_defs.parallel
        fprintf('\n This may take a while. Consider choosing Parallel=true for faster performance.\n' )
    end
    for jj = 1:length(mrgrad_defs.roi_names)
        seg_output_dir = mrGrad_seg(RG.(rg_names{jj}),mrgrad_defs.output_dir,true,mrgrad_defs.parallel);
    end
    if ~isempty(dir(fullfile(seg_output_dir,"**/*.nii.gz")))
        fprintf(' done!\n Segmentation files were saved in: %s\n',seg_output_dir);
    else
        fprintf(2,'\nAn error occurred while saving the output files. The results were not saved!\n');
    end
else
    fprintf(' NOTE: Individual subjects'' NIfTI outputs were not generated. Choose output_mode=''extended'' to generate files.\n');
end

clear mrgrad_defs
fprintf('\nAll done!\n');
end