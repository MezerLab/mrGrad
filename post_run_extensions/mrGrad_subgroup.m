function rg = mrGrad_subgroup(rg,idx,group_name)
% given an mrGrad rg struct output (an RG cell content), and a logical
% index, the function returns an rg struct of the indexed subgroup.


rg.Y = cellfun(@(x) x(:,idx),rg.Y,'un',0);
rg.Y_mean = cellfun(@(x) mean(x,2,"omitnan"),rg.Y,'un',0);
rg.Y_std = cellfun(@(x) std(x,0,2,"omitnan"),rg.Y,'un',0);
rg.Y_SEM = cellfun(@(x) std(x,0,2,"omitnan")/sqrt(size(x,2)),rg.Y,'un',0);
rg.individual_data = rg.individual_data(idx);
if exist("group_name","var")
    rg.group_name = group_name;
end

field_names = fieldnames(rg);
sub_descrips = field_names(cellfun(@(x) length(rg.(x))==length(idx),field_names));

% subjects' descriptive fields
sub_descrips(ismember(sub_descrips,{'Y','Y_mean','Y_std','Y_SEM','X','N_segments','y_lbls','group_name'})) = [];
for jj = 1:length(sub_descrips)
    rg.(sub_descrips{jj}) = rg.(sub_descrips{jj})(idx);
end
