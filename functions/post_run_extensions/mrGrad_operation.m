function RG_op = mrGrad_operation(RG,FunctionHandle)
% apply an element-wise function on mrGrad output values 
% FunctionHandle    A function handle to apply on matrices RG.Y
%                   Example: FunctionHandle = @(A) A-mean(A,1,"omitnan");

RG_op = cell(size(RG));
for jj = 1:numel(RG)
        rg = RG{jj};

        RG_op{jj}.Y = cellfun(FunctionHandle, rg.Y,'un',0);
        RG_op{jj}.Y_mean = cellfun(@(x) mean(x,2,"omitnan"),RG_op{jj}.Y,'un',0);
        RG_op{jj}.Y_std = cellfun(@(x) std(x,0,2,"omitnan"),RG_op{jj}.Y,'un',0);
        RG_op{jj}.Y_SEM = cellfun(@(x) std(x,0,2,"omitnan")/sqrt(size(x,2)),RG_op{jj}.Y,'un',0);

        additional_fields = fieldnames(rg);
        additional_fields(ismember(additional_fields,fieldnames(RG_op{jj}))) = [];
        for ff = 1:length(additional_fields)
            RG_op{jj}.(additional_fields{ff}) = rg.(additional_fields{ff});
        end
        
        if isfield(RG_op{jj},'individual_data')
            RG_op{jj}.WARNING_NOTE = sprintf('"individual_data" field is not manipulated by [%s]',func2str(FunctionHandle));
        end
end


