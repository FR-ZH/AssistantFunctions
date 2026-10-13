function alphaFreq = calc_channel_alpha_frequency(EEG, alphaBand)
% Output:
%   alphaFreq: channels × 1

if nargin < 2 || isempty(alphaBand)
    alphaBand = [8 13];
end

data = double(EEG.data);
fs = EEG.srate;

if ndims(data) == 3
    data = reshape(data, size(data,1), []);
end

nChan = size(data,1);
alphaFreq = nan(nChan,1);

winSec = 4;
win = round(winSec * fs);
noverlap = round(0.5 * win);
nfft = max(2^nextpow2(win), win);

for ch = 1:nChan
    x = data(ch,:);
    x = x - mean(x, 'omitnan');
    
    [pxx, f] = pwelch(x, hamming(win), noverlap, nfft, fs);
    
    idx = f >= alphaBand(1) & f <= alphaBand(2);
    
    if any(idx)
        alphaPxx = pxx(idx);
        alphaF = f(idx);
        [~, imax] = max(alphaPxx);
        alphaFreq(ch) = alphaF(imax);
    end
end

end