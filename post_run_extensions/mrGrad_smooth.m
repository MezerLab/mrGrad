function RG_smooth = mrGrad_smooth(RG,smoothing_window)

RG_smooth = RG;
for kk = 1:numel(RG)
    y = RG{kk}.Y;
    y1 = cellfun(@(y) smoothdata(y,1,"gaussian",smoothing_window),y,'un',0);
    RG_smooth{kk}.Y = y1;

    RG_smooth{kk}.Y_mean = cellfun(@(y1) mean(y1,2,"omitnan"),y1,'un',0);
    RG_smooth{kk}.Y_std = cellfun(@(y1) std(y1,0,2,"omitnan"),y1,'un',0);
    RG_smooth{kk}.Y_SEM = cellfun(@(y1) std(y1,0,2,"omitnan")./sqrt(size(y1,2)),y1,'un',0);
end
