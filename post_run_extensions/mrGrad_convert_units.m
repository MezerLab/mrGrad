function RG_conv = mrGrad_convert_units(RG,ConvertFunc,new_units)
    
RG_conv = RG;
for kk = 1:numel(RG)
    y_conv = cellfun(@(y) ConvertFunc(y) ,RG{kk}.Y,'un',0);
    RG_conv{kk}.Y = y_conv;

    RG_conv{kk}.Y_mean = cellfun(@(y1) mean(y1,2,"omitnan"),y_conv,'un',0);
    RG_conv{kk}.Y_std = cellfun(@(y1) std(y1,0,2,"omitnan"),y_conv,'un',0);
    RG_conv{kk}.Y_SEM = cellfun(@(y1) std(y1,0,2,"omitnan")./sqrt(size(y1,2)),y_conv,'un',0);
    RG_conv{kk}.units = new_units;
end