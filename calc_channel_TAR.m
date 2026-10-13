function TAR = calc_channel_TAR(EEG, thetaBand, alphaBand)
% Output:
%   TAR: channels × 1

if nargin < 2 || isempty(thetaBand)
    thetaBand = [4 8];
end
if nargin < 3 || isempty(alphaBand)
    alphaBand = [8 13];
end

data = double(EEG.data);
fs = EEG.srate;

if ndims(data) == 3
    data = reshape(data, size(data,1), []);
end

nChan = size(data,1);
TAR = nan(nChan,1);

winSec = 4;
win = round(winSec * fs);
noverlap = round(0.5 * win);
nfft = max(2^nextpow2(win), win);

for ch = 1:nChan
    x = data(ch,:);
    x = x - mean(x, 'omitnan');
    
    [pxx, f] = pwelch(x, hamming(win), noverlap, nfft, fs);
    
    thetaPow = bandpower_psd(pxx, f, thetaBand);
    alphaPow = bandpower_psd(pxx, f, alphaBand);
    
    TAR(ch) = thetaPow / alphaPow;
end

end