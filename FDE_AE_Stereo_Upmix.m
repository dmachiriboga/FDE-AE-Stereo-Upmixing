clc
clear

[stereo_signal, fs] = audioread('Stereo Test.mp3');

L = stereo_signal(:, 1); % Left channel
R = stereo_signal(:, 2); % Right channel



% Parameters for STFT
window_size = 1024;
overlap = window_size / 2;
nfft = window_size;

% Compute the STFT for both channels
X_L = stft(L, 'Window', hamming(window_size, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft);
X_R = stft(R, 'Window', hamming(window_size, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft);

% Presets
lambda = 0.85;  % Forgetting factor
sigma = 8; %slope
mu_0 = 0.1;
phi_0 = 0.15;

%Intialize vectors
phi_LR = zeros(size(X_L)); 
phi_LS = zeros(size(X_L)); 
phi_RS = zeros(size(X_L)); 
Phi = zeros(size(X_L)); 
u_c = zeros(size(X_L)); 
u_l = zeros(size(X_L)); 
u_r = zeros(size(X_L)); 
alpha = zeros(size(X_L)); 



% Compute the cross-correlation and coherence index

for m = 1:size(X_L, 2)
    for k = 1:size(X_L, 1)

        % Instantaneous cross-spectrum
        C_LR = X_L(k, m) * conj(X_R(k, m)); 
        C_LL = X_L(k, m) * conj(X_L(k, m));
        C_RR = X_R(k, m) * conj(X_R(k, m));
        
        % Cross-correlation function
        if m > 1
            phi_LR(k, m) = (1 - lambda) * phi_LR(k, m-1) + lambda * C_LR; %
            phi_LS(k, m) = (1 - lambda) * phi_LS(k, m-1) + lambda * C_LL;
            phi_RS(k, m) = (1 - lambda) * phi_RS(k, m-1) + lambda * C_RR;
        else
            phi_LR(k, m) = lambda * C_LR;
            phi_LS(k, m) = lambda * C_LL;
            phi_RS(k, m) = lambda * C_RR;
        end
        
        % Ambience Index
        if sqrt(phi_LS(k, m)*phi_RS(k, m)) == 0
            Phi(k, m) = 1 - abs(phi_LR(k, m))/(sqrt(phi_LS(k, m)*phi_RS(k, m)+0.001));
        else
            Phi(k, m) = 1 - abs(phi_LR(k, m))/(sqrt(phi_LS(k, m)*phi_RS(k, m)));
        end
        
        % Panning Index Value
        if C_LR == 0
             alpha(k, m) = 0; 
        else
            alpha(k, m) =  abs(asin( 2 * C_LR / (C_LL + C_RR)) / 2); 
        end

        % Panning Index
        u_c(k, m) = sin(2*alpha(k, m))/(sin(alpha(k, m))+cos(alpha(k, m)));  

        if C_LL-C_RR > 0
            u_l(k, m) =  cos(2*alpha(k, m)-pi)/sin(alpha(k, m)-pi/2);
        else 
            u_l(k, m) = 0;
        end

        if C_LL-C_RR < 0
            u_r(k, m) =  cos(2*alpha(k, m))/cos(alpha(k, m));
        else 
            u_r(k, m) = 0;
        end  
    end
end

% Apply the soft-decision function
Gamma = (1-mu_0)/2*tanh(sigma*pi*(Phi-phi_0))+(1+mu_0)/2;

% Reconstruct the signals
X_nL = X_L.* u_l;
X_C = X_L .* u_c + X_R .* u_c;
X_nR = X_R.* u_r;
X_LS = X_L .* Gamma;
X_RS = X_R .* Gamma;

% Inverse STFT
L_upmix = real(istft(X_nL, 'Window', hamming(window_size, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft));
C_upmix = real(istft(X_C, 'Window', hamming(window_size, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft));
R_upmix = real(istft(X_nR, 'Window', hamming(window_size, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft));
LS = real(istft(X_LS, 'Window', hamming(window_size, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft));
RS = real(istft(X_RS, 'Window', hamming(window_size, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft));

% All-pass delay and low pass
delay_time = 0.015; % Delay in seconds (5-20 ms)

% Compute the delay in terms of samples
delay_samples = round(delay_time * fs);

% Compute the all-pass filter coefficients
a = exp(-2 * pi * delay_samples / fs);
b = [1 -a];
a_filter = [1 -a]; 

% Delay ambience to decorrelate stereo image
LS_delay = filter(b, a_filter, LS);
RS_delay = filter(b, a_filter, RS);

% Low-pass filter
cutoff = 120; % Low-pass cutoff frequency in Hz
order = 4;   % Filter order

% Design a Butterworth low-pass filter at 120 Hz
[b_lp, a_lp] = butter(order, cutoff / (fs / 2), 'low');

% Apply low-pass filter to LS and RS signals
LS_upmix = filter(b_lp, a_lp, LS_delay);
RS_upmix = filter(b_lp, a_lp, RS_delay);


stereo(:,1) = L_upmix;
stereo(:,2) = R_upmix;
rear(:,1) = LS_delay;
rear(:,2) = RS_delay;

% Find the maximum absolute value in the audio signal
max_val = max(abs(C_upmix(:)));

% Normalize the audio signal if necessary
if max_val > 1
    C_upmix = C_upmix / max_val;
end

audiowrite("LR.mp3", stereo, fs)
audiowrite("LSRS.mp3", rear, fs)
audiowrite("center.mp3", C_upmix, fs)