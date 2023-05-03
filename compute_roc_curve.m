function [fpr, tpr] = compute_roc_curve(x0_axis, x0_freq, x1_axis, x1_freq, bin_size)
    %Compute the ROC curve from the frequency densities of the two classes
    
    min_tot = min([min(x0_axis), min(x1_axis)]) - bin_size;
    max_tot = max([max(x0_axis), max(x1_axis)]);
    
    th = [min_tot : bin_size : max_tot];
    fpr = zeros(size(th));
    tpr = zeros(size(th));
    
    i_th = 1;
    for th_i = th
        fpr(i_th) = sum(x0_freq(x0_axis > th_i));
        tpr(i_th) = sum(x1_freq(x1_axis > th_i));
        i_th = i_th + 1;
    end
    fpr = fpr(end:-1:1);
    tpr = tpr(end:-1:1);
    
end


