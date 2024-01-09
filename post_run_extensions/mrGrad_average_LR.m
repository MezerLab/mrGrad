function RG_avg = mrGrad_average_LR(RG,hemisphere_notations)

% input: mrgrad cell contianing left, right ROIs
% output: average mrgrad struct


if ~exist('hemisphere_notations','var')
    hemisphere_notations = {'Left-','Right-'};
end
roi_labels = cellfun(@(x) x.ROI_label,RG(:),'un',0);
region_names = unique(cellfun(@(x) erase(x,hemisphere_notations),roi_labels,'un',0),'stable')';

ngroups = size(RG,1);
nregions = length(region_names);

RG_avg = cell(ngroups,nregions);
for gg = 1:ngroups
    for rr = 1:nregions
        rg = RG(gg,:);
        rg(cellfun(@(x) ~contains(x.ROI_label,region_names{rr}),rg)) = [];
        if numel(rg)~=2
            error('number of hemispheres per regions are different than 2')
        end
        RG_avg{gg,rr}.Y = cellfun(@(a,b) (a+b)/2, rg{1}.Y, rg{2}.Y,'un',0);
        RG_avg{gg,rr}.Y_mean = cellfun(@(x) mean(x,2,"omitnan"),RG_avg{gg,rr}.Y,'un',0);

        RG_avg{gg,rr}.Y_std = cellfun(@(x) std(x,0,2,"omitnan"),RG_avg{gg,rr}.Y,'un',0);
        RG_avg{gg,rr}.Y_SEM = cellfun(@(x) std(x,0,2,"omitnan")/sqrt(size(x,2)),RG_avg{gg,rr}.Y,'un',0);

        additional_fields = fieldnames(rg{1});
        additional_fields(ismember(additional_fields,fieldnames(RG_avg{gg,rr}))) = [];
        for ff = 1:length(additional_fields)
            RG_avg{gg,rr}.(additional_fields{ff}) = rg{1}.(additional_fields{ff});
        end
        RG_avg{gg,rr}.ROI_label = sprintf('LR-average-%s',region_names{rr});
        
        if isfield(RG_avg{gg,rr},'individual_data')
            RG_avg{gg,rr}.WARNING_NOTE = sprintf('"individual_data" field is of %s side only',rg{1}.ROI_label);
        end
    end
end


