function [p1_output, f] = getFft(rawS)
%%Discrete Fourier transform.
%
%   [p1_output, f] = getFft(rawS)
%
%   Description: Discrete Fourier transform of multiple photometry signals,
%   previously extracted from data files using extractRaw_fft.
%   The function will pad each signal such that the vector to be analyzed
%   is of length 2500 x the sampling frequency. The single-sided spectrum
%   P1 will be extracted and smoothed before being concatenated into an
%   output matrix, p1_output.
%
%   INPUT
%   'rawS' - structure generated using extractRaw_fft
%       must contain fields rawS(x).fp_sub and rawS(x).rawFs
%
%   OUTPUT
%   'p1_output' - matrix with concatenated smoothed fft outputs, where each
%       column contains the output for each recording in rawS
%   'f' - vector with frequency domain
%
%   Anya Krok, September 2022
%
p1_output = [];
h = waitbar(0, 'FFT photometry signals');
for x = 1:length(rawS)
    vec = [rawS(x).fp_sub]; 
    Fs = rawS(x).rawFs;
    needL = 2500*Fs;
    if isempty(vec) 
        p1_output(:,x) = nan(1+(needL/2),1);
        continue
    end
    vec = repmat(vec,[ceil(needL/length(vec)) 1]);
    vec = vec(1:needL);
    T = 1/Fs;               % Sampling period
    L = length(vec);        % Length of signal
    vec(isnan(vec)) = [];
    fftACh = fft(vec);      % Discrete Fourier Transform of photometry signal
    P2 = abs(fftACh/L);     % Two-sided spectrum P2
    P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;     % Frequency domain vector
    P1 = medfilt1(P1);      % Median filter initial FFT
    P1 = movmean(P1,500);   % Smooth FFT output
    p1_output(:,x) = P1; % Concatenate output
    waitbar(x/length(rawS),h);
end
close(h)