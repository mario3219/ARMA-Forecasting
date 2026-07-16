function y_filt = myMedianFilter(y, windowSize, showRes, str)

    halfWin = floor(windowSize/2);
    y_filt  = y;

    for i = 1:length(y)
        low  = max(1, i-halfWin);
        high = min(length(y), i+halfWin);

        y_filt(i) = median(y(low:high));
    end

    if showRes == 1
        figure; sgtitle(str)
        subplot(211); plot(y);      title("Before");
        subplot(212); plot(y_filt); title("After");
    end

end
