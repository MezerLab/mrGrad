function output_dir = mrGrad_seg(rg,output_dir,force_flag,Parallel)

force_flag = exist("force_flag","var") && force_flag;
Parallel = exist("Parallel","var") && Parallel;

seg_list = rg.user_input_fields.seg_list;
individual_data = rg.individual_data;
roi_name = lower(string(rg.roi_name));
axis_names = lower(string(rg.axis_names));
n_segments = rg.n_segments;
% parameter_name = lower(string(rg.parameter));
subject_names = string(rg.subject_ids);

invalidPattern = '[\\/:*?"<>|]|[\x00-\x1F]';

% make sure subject names are unique and ordered as input, for using as
% output directories
subject_names_unique = unique(subject_names,"stable");
if ~isequal(subject_names,subject_names_unique)
    maxDigits = max([floor(log10(length(subject_names)))+1,3]);
    fmt = ['sub-%0' num2str(maxDigits) 'd_%s'];
    subject_names = compose(fmt, (1:length(subject_names))', subject_names);
end
output_dir = fullfile(output_dir,"mrGradSeg");

if Parallel
    parfor ii = 1:length(seg_list)
        sub_info = individual_data{ii};

        % Get original image info
        if isfield(sub_info,'original_nifti_info')
            image_info = sub_info(1).original_nifti_info;
        elseif exist(seg_list{ii},'file')
            image_info = niftiinfo(seg_list{ii});
        else
            % data not exist for subject ii
            continue
        end
        [strides_orig,im_dims_orig] = keep_strides(image_info);
        image_size_orig = image_info.ImageSize;
        im_empty_orig = zeros(image_size_orig);

        % Get standard image [+1,+2,+3] info
        [~, im_perm] = sort(im_dims_orig);
        im_empty_std = permute(im_empty_orig, im_perm);  % permute to standard order

        image_size_std = size(im_empty_std);
        dimsflip = im_dims_orig(strides_orig < 0);

        sub_outdir = fullfile(output_dir,regexprep(subject_names{ii},invalidPattern,'_'));
        for Axis = 1:length(sub_info)
            filename = sprintf('mrGradSeg_%s_%s_%dsegments',roi_name,axis_names{Axis},n_segments(Axis));
            filename = regexprep(filename,invalidPattern,'');

            filepath = fullfile(sub_outdir,filename);
            if exist(filepath+".nii.gz",'file') && ~force_flag
                continue
            end

            coord_list_std = sub_info(Axis).segment_inds_linear;
            seg_std = mrgrad_coords2seg(coord_list_std,image_size_std);

            seg_restored = seg_std;
            % Flip back originally-negative strides
            for d = dimsflip
                seg_restored = flip(seg_restored,d);
            end

            % Invert the original permutation
            [~, im_perm_inv] = sort(im_perm);
            seg_restored = permute(seg_restored, im_perm_inv);

            if ~exist(sub_outdir,'dir')
                mkdir(sub_outdir);
            end

            image_info.Datatype = 'uint8';
            niftiwrite(seg_restored,filepath,image_info,Compressed=true);
        end
    end

else
    for ii = 1:length(seg_list)
        sub_info = individual_data{ii};

        % Get original image info
        if isfield(sub_info,'original_nifti_info')
            image_info = sub_info(1).original_nifti_info;
        elseif exist(seg_list{ii},'file')
            image_info = niftiinfo(seg_list{ii});
        else
            % data not exist for subject ii
            continue
        end
        [strides_orig,im_dims_orig] = keep_strides(image_info);
        image_size_orig = image_info.ImageSize;
        im_empty_orig = zeros(image_size_orig);

        % Get standard image [+1,+2,+3] info
        [~, im_perm] = sort(im_dims_orig);
        im_empty_std = permute(im_empty_orig, im_perm);  % permute to standard order

        image_size_std = size(im_empty_std);
        dimsflip = im_dims_orig(strides_orig < 0);

        sub_outdir = fullfile(output_dir,regexprep(subject_names{ii},invalidPattern,'_'));
        for Axis = 1:length(sub_info)
            filename = sprintf('mrGradSeg_%s_%s_%dsegments',roi_name,axis_names{Axis},n_segments(Axis));
            filename = regexprep(filename,invalidPattern,'');

            filepath = fullfile(sub_outdir,filename);
            if exist(filepath+".nii.gz",'file') && ~force_flag
                continue
            end

            coord_list_std = sub_info(Axis).segment_inds_linear;
            seg_std = mrgrad_coords2seg(coord_list_std,image_size_std);

            seg_restored = seg_std;
            % Flip back originally-negative strides
            for d = dimsflip
                seg_restored = flip(seg_restored,d);
            end

            % Invert the original permutation
            [~, im_perm_inv] = sort(im_perm);
            seg_restored = permute(seg_restored, im_perm_inv);

            if ~exist(sub_outdir,'dir')
                mkdir(sub_outdir);
            end

            image_info.Datatype = 'uint8';
            niftiwrite(seg_restored,filepath,image_info,Compressed=true);
        end
    end
end
