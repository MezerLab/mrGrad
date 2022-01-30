function RG = RG_flip(RG,pcs)
% This function receives a RG struct of a dataset and flips the segment
% order along the requested pc(s)
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021

for pc = pcs
    if isfield(RG,'Y')
        RG.Y{pc} = flip(RG.Y{pc});
    end
    if isfield(RG,'Y_mean')
        RG.Y_mean{pc} = flip(RG.Y_mean{pc});
    end
    if isfield(RG,'Y_std')
        RG.Y_std{pc} = flip(RG.Y_std{pc});
    end
    if isfield(RG,'Y_SEM')
        RG.Y_SEM{pc} = flip(RG.Y_SEM{pc});
    end
    Nsubs = length(RG.individual_data);
    for ii = 1:Nsubs
        RG.individual_data{ii}(pc).segment_inds = flip(RG.individual_data{ii}(pc).segment_inds);
        RG.individual_data{ii}(pc).segment_inds_linear = flip(RG.individual_data{ii}(pc).segment_inds_linear);
        RG.individual_data{ii}(pc).planes.gz = flip(RG.individual_data{ii}(pc).planes.gz);
        RG.individual_data{ii}(pc).analysisPC_line = flip(RG.individual_data{ii}(pc).analysisPC_line, 2);
        RG.individual_data{ii}(pc).function = flip(RG.individual_data{ii}(pc).function);
    end
end