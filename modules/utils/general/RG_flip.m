function RG = RG_flip(RG,pcs)
% This function receives a RG struct and flips the segment order along the requested pc(s)
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021

for pc = pcs
    if isfield(RG,'Y') % mrgrad v-1
        RG.Y{pc} = flip(RG.Y{pc});
    end
    if isfield(RG,'Y_mean') % mrgrad v-1
        RG.Y_mean{pc} = flip(RG.Y_mean{pc});
    end
    if isfield(RG,'Y_std') % mrgrad v-1
        RG.Y_std{pc} = flip(RG.Y_std{pc});
    end
    if isfield(RG,'Y_SEM') % mrgrad v-1
        RG.Y_SEM{pc} = flip(RG.Y_SEM{pc});
    end
    if isfield(RG,'Results') % mrgrad v-2
        results_field_names = fieldnames(RG.Results);
        results_field_names = results_field_names(...
            contains(results_field_names,"axis"+pc));
        for jj = 1:length(results_field_names)

            RG.Results.(results_field_names{jj}).y{:,:} = flip(...
                RG.Results.(results_field_names{jj}).y{:,:}, 2);

            for f = ["y_mean", "y_std", "y_sem"]
                if isfield(RG.Results.(results_field_names{jj}),f)
                     RG.Results.(results_field_names{jj}).(f) = flip(...
                         RG.Results.(results_field_names{jj}).(f), 2);
                end
            end
        end
    end
    if isfield(RG,"individual_data")
        for ii = 1:length(RG.individual_data)
            if ~isempty(RG.individual_data{ii})
                RG.individual_data{ii}(pc).segment_inds = flip(RG.individual_data{ii}(pc).segment_inds);
                RG.individual_data{ii}(pc).segment_inds_linear = flip(RG.individual_data{ii}(pc).segment_inds_linear);
                RG.individual_data{ii}(pc).planes.gz = flip(RG.individual_data{ii}(pc).planes.gz);
                RG.individual_data{ii}(pc).analysisPC_line = flip(RG.individual_data{ii}(pc).analysisPC_line, 2);
                RG.individual_data{ii}(pc).function = flip(RG.individual_data{ii}(pc).function,1);
            end
        end
    end
end