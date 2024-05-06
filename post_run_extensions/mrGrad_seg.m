
function mrGrad_seg(rg,output_dir,force_flag)

% REQUIREMENT
%        * Vistasoft       - https://github.com/vistalab/vistasoft   

force_flag = exist("force_flag","var") && force_flag;

seg_list = rg.seg_list;

parfor ii = 1:length(seg_list)
    seg = niftiread(seg_list{ii});


%     seg = readFileNifti(seg_list{ii});
    strides = keep_strides(seg_list{ii});
    dims = 1:3;
    dimsflip = dims(strides<0);
    for d = dimsflip
        seg = flip(seg,d);
    end

    dat = rg.individual_data{ii};
    
    for pc = 1:length(dat)
        filename = sprintf('Seg_mrGrad_%s_%s_%dsegments',rg.ROI_label,rg.y_lbls{pc},rg.N_segments(pc));
        outdir = fullfile(output_dir,rg.subject_names{ii},'analysis/segmentation/mrGrad');
        filepath = fullfile(outdir,filename);
        if exist(filepath,'file') && ~force_flag
            continue
        end

        coords = dat(pc).segment_inds_linear;
        newseg = single(zeros(size(seg)));
        for jj = 1:length(coords)
            newseg(coords{jj}) = jj;
        end
        for d = dimsflip
            newseg = flip(newseg,d);
        end
        
        if ~exist(outdir,'dir')
            mkdir(outdir);
        end
        info = niftiinfo(seg_list{ii});
        info.Datatype = 'single';
        niftiwrite(newseg,filepath,info,Compressed=true);
%         dtiWriteNiftiWrapper(newseg,seg.qto_xyz,filepath);
    end
end

