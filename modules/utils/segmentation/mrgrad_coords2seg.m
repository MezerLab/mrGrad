function seg = mrgrad_coords2seg(coord_list,image_size)

seg = single(zeros(image_size));
for jj = 1:length(coord_list)
    seg(coord_list{jj}) = jj;
end
