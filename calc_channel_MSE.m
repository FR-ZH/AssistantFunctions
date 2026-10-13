function MSE = calc_channel_MSE(EEG, scales, m, r_ratio)
%   Output:
%   MSE: channels × scales

if nargin < 2 || isempty(scales)
    scales = 10:10:50;
end
if nargin < 3 || isempty(m)
    m = 2;
end
if nargin < 4 || isempty(r_ratio)
    r_ratio = 0.15;
end

data = double(EEG.data);

if ndims(data) == 3
    data = reshape(data, size(data,1), []);
end

nChan = size(data,1);
nScale = numel(scales);
MSE = nan(nChan, nScale);

for ch = 1:nChan
    x = data(ch,:);
    x = x - mean(x, 'omitnan');

    if std(x, 'omitnan') == 0
        continue
    end

    for s = 1:nScale
        y = coarse_grain_signal(x, scales(s));
        r = r_ratio * std(y, 'omitnan');
        MSE(ch,s) = sample_entropy(y, m, r);
    end
end

end