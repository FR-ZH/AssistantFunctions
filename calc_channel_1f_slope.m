function slope = calc_channel_1f_slope(EEG, fitRange, excludeAlpha)
% Output:
%   slope: channels × 1

if nargin < 2 || isempty(fitRange)
    fitRange = [2 30];
end
if nargin < 3 || isempty(excludeAlpha)
    excludeAlpha = true;
end

data = double(EEG.data);
fs = EEG.srate;

if ndims(data) == 3
    data = reshape(data, size(data,1), []);
end

nChan = size(data,1);
slope = nan(nChan,1);

winSec = 4;
win = round(winSec * fs);
noverlap = round(0.5 * win);
nfft = max(2^nextpow2(win), win);

for ch = 1:nChan
    x = data(ch,:);
    x = x - mean(x, 'omitnan');
    
    [pxx, f] = pwelch(x, hamming(win), noverlap, nfft, fs);
    
    idx = f >= fitRange(1) & f <= fitRange(2);
    
    if excludeAlpha
        idx = idx & ~(f >= 8 & f <= 13);
    end
    
    logF = log10(f(idx));
    logP = log10(pxx(idx));
    
    good = isfinite(logF) & isfinite(logP);
    
    if sum(good) > 5
        beta = polyfit(logF(good), logP(good), 1);
        slope(ch) = beta(1);
    end
end

end