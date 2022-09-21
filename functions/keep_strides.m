function strides = keep_strides(nifti_struct)
    strides = diag(nifti_struct.qto_xyz(1:3,1:3));
    strides = strides./abs(strides);

end