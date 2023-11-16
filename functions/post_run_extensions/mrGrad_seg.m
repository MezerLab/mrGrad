
function mrGrad_seg(rg,output_dir)

force_flag = 0;
seg_list = rg.seg_list;

for ii = 1:length(seg_list)
    seg = readFileNifti(seg_list{ii});
    strides = keep_strides(seg);
    dims = 1:3;
    dimsflip = dims(strides<0);
    for d = dimsflip
        seg.data = flip(seg.data,d);
    end

    dat = rg.individual_data{ii};
    
    for pc = 1:length(dat)
        filename = sprintf('Seg_mrGrad_%s_%s_%dsegments',rg.ROI_label,rg.y_lbls{pc},rg.N_segments(pc));
        outdir = fullfile(output_dir,rg.subject_names{ii},'analysis/Seg_mrGrad');
        filepath = fullfile(outdir,filename);
        if exist(filepath,'file') && ~force_flag
            continue
        end

        coords = dat(pc).segment_inds_linear;
        newseg = zeros(size(seg.data));
        for jj = 1:length(coords)
            newseg(coords{jj}) = jj;
        end
        for d = dimsflip
            newseg = flip(newseg,d);
        end
        
        if ~exist(outdir,'dir')
            mkdir(outdir);
        end
        dtiWriteNiftiWrapper(newseg,seg.qto_xyz,filepath);
    end
end

