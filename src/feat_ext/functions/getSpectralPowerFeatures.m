function feat_vector = getSpectralPowerFeatures(x,f,Fs)
% Function calculates the power spectral density in decibels
% of a given signal different frequencies
%   x  - Vector containing the signal of interest
%   f  - Vector containing the frequencies at which to calculate the PSD
%   Fs - The sampling frequency

[Pxx,~] = pwelch(double(x),[],[],f,Fs);
feat_vector = 10*log10(Pxx);
