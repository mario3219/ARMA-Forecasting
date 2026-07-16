function y = removeOutliers(x,n_samples,showRes,str)
    
    % Store original values for comparison
    x_b = x;

    % Append an index column to x1
    x = [(1:length(x))' x];
    x = sortrows(x,-2);           % Sort x1 according to data
                                  % from highest to lowest

    % Set highest and lowest values to NaN
    x(end-n_samples:end,2) = NaN;
    x(1:1+n_samples,2) = NaN;
    
    % Sort according to index, reforming same
    % ordering as original
    x = sortrows(x,1);
    x = x(:,2);

    % Interpolate Nan values
    x_new = fillmissing(x, 'linear');

    if showRes == 1
        low = min(x_b); high = max(x_b);
        figure; sgtitle(str);
        subplot(311); plot(x_b);    ylim([low high]);    title("Before")
        subplot(312); plot(x);      ylim([low high]);    title("Outliers removed (" +n_samples+" of highest and lowest removed)")
        subplot(313); plot(x_new);  ylim([low high]);    title("Interpolated")
    end
    y = x_new;
end