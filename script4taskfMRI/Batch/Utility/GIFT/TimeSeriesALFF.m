
function [ALFF,mALFF]=TimeSeriesALFF(TimeSeries,ASamplePeriod,LowCutoff,HighCutoff)

%-------------------------------------------------------------------------%

% Frequency Index
fprintf('\n\t Frequency Index ...');
sampleFreq 	 = 1/ASamplePeriod;
sampleLength = size(TimeSeries,1);
paddedLength = rest_nextpow2_one35(sampleLength); %2^nextpow2(sampleLength);
if (LowCutoff >= sampleFreq/2) % All high included
    idx_LowCutoff = paddedLength/2 + 1;
else % high cut off, such as freq > 0.01 Hz
    idx_LowCutoff = ceil(LowCutoff * paddedLength * ASamplePeriod + 1);
    % Change from round to ceil: idx_LowCutoff = round(LowCutoff *paddedLength *ASamplePeriod + 1);
end
if (HighCutoff>=sampleFreq/2)||(HighCutoff==0) % All low pass
    idx_HighCutoff = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
    idx_HighCutoff = fix(HighCutoff *paddedLength *ASamplePeriod + 1);
    % Change from round to fix: idx_HighCutoff	=round(HighCutoff *paddedLength *ASamplePeriod + 1);
end

% Detrend
TimeSeries=detrend(TimeSeries);

% Zero-Paddings
TimeSeries=[TimeSeries;zeros(paddedLength -sampleLength,size(TimeSeries,2))];

% FFT
fprintf('\n\t Performing FFT ...');
ALFF = 2*abs(fft(TimeSeries))/sampleLength;

% Mean ALFF
mALFF = mean(ALFF(idx_LowCutoff:idx_HighCutoff,:));
fprintf('\n\t Done\n');

return
end