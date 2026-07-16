% Returns an array with lags that break confidence interval in ACF or PACF
% Only the 10 most significant

% First column:  significant lag
% Second column: by how much

function result = getSignLags(corr, N)

len = length(corr);
result = [];
conf = 2/sqrt(N);
for idx = 1:len
    if abs(corr(idx)) > conf
        result = [result; [idx-1 abs(corr(idx))] ];
    end
end
result = sortrows(result,-2);
if length(result)>5
    result = result(1:5,:);
end
