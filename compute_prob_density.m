function [x_hist, y_bar] = compute_prob_density(x, bin_size)
    %compute frequency density
    
    x_hist = [min(x) : bin_size : max(x)];
    y_bar = zeros(size(x_hist));
    i_bin = 1;
    for bin_i = x_hist
        y_bar(i_bin) = sum(x >= bin_i & x < bin_i+bin_size);
        i_bin = i_bin + 1;
    end
    y_bar = y_bar / length(x);
end
