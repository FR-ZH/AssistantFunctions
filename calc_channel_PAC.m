function PAC = calc_channel_PAC(EEG, phaseBand, ampBand, nBins)
% Output:
%   PAC: channels × 1

if nargin < 2 || isempty(phaseBand)
    phaseBand = [4 8];
end
if nargin < 3 || isempty(ampBand)
    ampBand = [30 45];
end
if nargin < 4 || isempty(nBins)
    nBins = 18;
end

data = double(EEG.data);
fs = EEG.srate;

if ndims(data) == 3
    data = reshape(data, size(data,1), []);
end

nChan = size(data,1);
PAC = nan(nChan,1);

for ch = 1:nChan
    x = data(ch,:);
    x = x - mean(x, 'omitnan');
    
    if std(x, 'omitnan') == 0
        continue
    end
    
    xPhase = bandpass_eeg(x, fs, phaseBand);
    xAmp   = bandpass_eeg(x, fs, ampBand);
    
    phase = angle(hilbert(xPhase));
    ampEnv = abs(hilbert(xAmp));
    
    edges = linspace(-pi, pi, nBins + 1);
    meanAmp = nan(1, nBins);
    
    for b = 1:nBins
        idx = phase >= edges(b) & phase < edges(b+1);
        meanAmp(b) = mean(ampEnv(idx), 'omitnan');
    end
    
    P = meanAmp ./ nansum(meanAmp);
    
    if any(P <= 0) || any(~isfinite(P))
        PAC(ch) = NaN;
    else
        H = -sum(P .* log(P));
        PAC(ch) = (log(nBins) - H) / log(nBins);
    end
end

end