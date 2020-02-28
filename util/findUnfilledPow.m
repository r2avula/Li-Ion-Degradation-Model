function [unfilled_pow_idx] = findUnfilledPow(energy_loss_samples,cur_soc_bin_idx)
pow_num = size(energy_loss_samples,2);
valid = zeros(1,pow_num);
for j = 1:pow_num
    if(~isempty(find(isinf(energy_loss_samples(cur_soc_bin_idx,j,:)) ==1, 1)))
        valid(j) = 1;
    end
end
found_indices = find(valid == 1);
new_num = length(found_indices);
if(new_num>0)
    unfilled_pow_idx = found_indices(unidrnd(new_num));
else
    unfilled_pow_idx = [];
end
end

