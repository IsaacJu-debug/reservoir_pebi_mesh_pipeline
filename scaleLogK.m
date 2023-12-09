function unstruct_k_map = scaleLogK(k_map, logk_mean, logk_std)
%scaleLogk Summary of this function goes here
    unstruct_k_map = exp(k_map * logk_std + logk_mean);
end

