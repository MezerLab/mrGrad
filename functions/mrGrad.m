function [RG, T] = mrGrad(Data,varargin)
%% mrGrad v2.0: MRI Region Gardients
%
%--------------------------------------------------------------------------
% INPUTS:
%--------------------------------------------------------------------------
%
%   Data :  struct with the following fields:
%
%           - 'seg_list' : (cell/strings of size N x 1)
%                          Paths to segmentation or binary mask files, in
%                          length N (number of subjects)
%
%           - 'map_list' : (cell/strings of size N x M)
%                          Paths to subjects MRI images (nii).
%                          Can have multiple (M) image types (e.g., T1, T2)
%                          of the same size and subject space as seg_list.
%
%           - 'subject_ids': (cell/strings/numeric of size N)
%                            Optional. Unique identifiers for each subject
%                            record (e.g., {'sub-001', 'sub-002'}).
%
%           - Any other fields with the same length as map_list/seg_list
%             (e.g., 'Age', 'Sex', 'ResearchGroup'). These fields will be
%             copied to the output struct for reference.
%
%--------------------------------------------------------------------------
% OPTIONAL Name-Value Pair Arguments:
%--------------------------------------------------------------------------
%
%   'ROI'  : (vector) list of FreeSurfer/FSL label indices (e.g. [11 50 12 51] for
%           l-caudate, r-caudate, l-putamen r-putamen);
%           In case the provided labels do not refer to FreeSurfer's
%           look-up table, please provide the 'roi_names' argument.
%
%   'roi_names'          : (cell array of strings/char) names corresponding to ROI labels.
%                          Required if labels do not match FreeSurfer LUT.
%
%   'Axes' / 'PC'        : (vector) list of principal component (PC) indices
%                          (axes) to compute [1 2 3]. Default: all.
%
%   'n_segments'         : (scalar or vector) number of segments per axis
%                          Default: 7.
%
%   'segmenting_method'  : (string/char) 'equidistance' (default) or 'equivolume'.
%
%   'stat'               : (string/char) statistic used for the axis function:
%                          'median' (default) or 'mean'.
%
%   'max_change' : (numeric array) [nROIs x nPCs], where each value specifies
%                  the image axis (1 = X, 2 = Y, 3 = Z) most aligned with a
%                  principal component (PC) axis for each ROI.
%
%                  This ensures consistent **directionality** of ROI axes
%                  across subjects. Without it, PC directions may vary in sign
%                  (e.g., anterior→posterior in one subject, posterior→anterior
%                  in another). With 'max_change', direction is enforced from
%                  lower to higher coordinate values along the given image axis.
%
%                  Example: max_change = [2 3 1] for a striatal ROI
%                  - PC1 (AP) aligns with Y-axis (A >> P)
%                  - PC2 (VD) aligns with Z-axis (V >> D)
%                  - PC3 (ML) aligns with X-axis (M >> L)
%
%   'erode'              : (logical) remove outer voxel layer to reduce partial volume effects.
%                          Default: false (0).
%
%   'parameter_names'    : (string/char/cell) names of the MRI parameters being
%                          analyzed (e.g., 'R1').
%
%   'units'              : (string/char/cell) unit of measurements (e.g., 'sec^{-1}').
%
%   'apply_alternative_axes' : (struct) define axes using a different ROI than the one analyzed.
%                               Fields:
%                               * seg_list : cell array of alternative segmentations.
%                               * ROI      : ROI label(s) for axis computation.
%                               NOTE that this option should be used carefully as the used
%                               axes would not carry anatomical meaning in the sampled ROI.
%
%   'allow_missing'     : (logical) skip missing subjects data instead of error. Default: false.
%
%   'output_name'        : (string/char) name for the saved summary file. If a full path is
%                          provided and 'output_dir' is also provided, mrGrad will
%                          use the 'output_dir' as the output path.
%
%   'output_dir'         : (string/char) directory for output files.
%
%   'output_mode'        : (string/char) one of:
%                          * 'minimal'  - only summarized results (.csv, .mat)
%                          * 'default'  - includes per-subject axis data (larger .mat file)
%                          * 'extended' - also saves per-subject axis segmentation masks (.nii.gz)
%
%   'Parallel'           : (logical) use MATLAB parallel pool (PARFOR loop)
%
%--------------------------------------------------------------------------
% OUTPUTS:
%--------------------------------------------------------------------------
%   RG : A struct with 1 x nROIs fields
%       Each field contains a struct with ROI-specific gradient results
%       * Gradient values per segment and axis
%       * Mean/STD/SEM per axis
%       * subject-level results
%       * Metadata (e.g., parameter name, ROI label)
%
%   T  : Combined summary table (compiled from all ROIs and groups)
%
%--------------------------------------------------------------------------
% DEPENDENCIES:
%--------------------------------------------------------------------------
%   - MATLAB (https://www.mathworks.com)
%   - boundedline-pkg (recommended, https://github.com/kakearney/boundedline-pkg)
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

if ~mrgrad_defs.parallel
    fprintf(2,' May run slowly on large cohorts — for faster performance set Parallel=true.');
end

% check obligatory input
[Data, mrgrad_defs] = mrgrad_check_input(Data,mrgrad_defs);
num_regions = numel(mrgrad_defs.ROI);
num_subjects = length(Data.seg_list);
parameter_str = strjoin(strrep(unique(cellstr(mrgrad_defs.parameter_names)),'unknown_parameter','unknown parameter'),", ");
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
    if mrgrad_defs.parallel
        gcp();
        parfor ii = 1:num_subjects
            Allsubs_rg_data{ii} = mrgrad_inner_loop(ii, rr, seg_list, map_list, mrgrad_defs);
        end
    else
        for ii = 1:num_subjects
            Allsubs_rg_data{ii} = mrgrad_inner_loop(ii, rr, seg_list, map_list, mrgrad_defs);
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
        summary_roi_name = regexprep(lower(rg.roi_name),{'-',' '},'_');
        summary_vars = "seg"+(1:n_segments(kk));

        for pp = 1:length(mrgrad_defs.parameter_names)
            summary_name = strjoin({rg.axis_names{kk},mrgrad_defs.parameter_names{pp}},"_");
            summary_mat = cellfun(@(x) x(kk).function(:,pp)',Allsubs_rg_data,'un',0);
            summary_mat = cat(1,summary_mat{:});

            rg.Results.(summary_name).x = (1:n_segments(kk))/n_segments(kk);
            rg.Results.(summary_name).y = ...
                array2table(summary_mat,VariableNames=summary_vars,...
                RowNames=summary_rows);
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