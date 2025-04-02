function [strides,dim_order,strides_str] = keep_strides(nifti_path)
    % input should be either a nifti path or a niftiinfo() struct
    
    if isa(nifti_path,"struct")
        mat = nifti_path.Transform.T(1:3,1:3)';
    else
        mat = niftiinfo(nifti_path).Transform.T(1:3,1:3)';
    end
    mat_new = zeros(size(mat));
    for jj = 1:3
        [~,i] = max(abs(mat(jj,:)));
        mat_new(jj,i) = mat(jj,i);
    end
    mat = mat_new;
    
    dim_order = zeros(1,3);
    SignNeg = false(1,3);
    for jj = 1:3
        dim_order(jj) = find(mat(jj,:));
        SignNeg(jj) = mat(jj,dim_order(jj)) < 0;
    end

    strides = dim_order;
    strides(SignNeg) = -1*dim_order(SignNeg);
    strides_str = strjoin(string(strides),",");
end