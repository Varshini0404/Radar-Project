clc; clear; close all;

%% ============================================================
%  FMCW RADAR RANGE–DOPPLER MAP (FULL PROFESSIONAL VERSION)
%  1) Noiseless output
%  2) Noisy output for multiple SNRs
%  Includes: Windowing + Peak detection + Proper axis labels
%% ============================================================

%% =========================
%  STEP 1: Radar Parameters
% =========================
c  = 3e8;
fc = 77e9;
lambda = c/fc;

Rmax = 200;        % max range (m)
dR   = 1;          % range resolution (m)
Vmax = 70;         % max velocity (m/s)

B = c/(2*dR);
Tchirp = 5.5 * (2*Rmax/c);
S = B/Tchirp;

%% ===========================
%  STEP 2: Sampling parameters
% ===========================
Nr = 1024;         % samples per chirp
Nd = 128;          % chirps

t = linspace(0, Nd*Tchirp, Nr*Nd);

%% ======================
%  STEP 3: Target settings
% ======================
R0 = 110;          % initial range (m)
v  = -20;          % velocity (m/s) negative = approaching

R_t = R0 + v*t;
tau = 2*R_t/c;

%% ===============================
%  STEP 4: Tx and Rx (FMCW chirp)
% ===============================
Tx = cos(2*pi*(fc*t + (S*t.^2)/2));

t_rx = t - tau;
Rx = cos(2*pi*(fc*t_rx + (S*t_rx.^2)/2));

%% ======================
%  STEP 5: Beat signal
% ======================
Mix_clean = Tx .* Rx;

%% =========================
%  Common Axes for plots
% =========================
range_axis = (0:(Nr/2)-1) * (c/(2*B));     % meters
doppler_axis = linspace(-Vmax, Vmax, Nd); % m/s

%% ============================================================
%  FUNCTION: PROCESS + PLOT RANGE FFT AND RANGE DOPPLER MAP
%% ============================================================
function process_and_plot(Mix, Nr, Nd, range_axis, doppler_axis, title_tag)

  % Reshape
  Mix_reshaped = reshape(Mix, Nr, Nd);

  % 2D Windowing (Cleaner FFT)
  wr = hamming(Nr);
  wd = hamming(Nd);
  W = wr * wd';
  Mix_windowed = Mix_reshaped .* W;

  %% ----------------------
  % Range FFT (1D FFT)
  %% ----------------------
  sig_fft1 = fft(Mix_windowed, Nr);
  sig_fft1 = abs(sig_fft1/Nr);
  sig_fft1 = sig_fft1(1:Nr/2, :);

  figure;
  plot(range_axis, sig_fft1(:,1));
  grid on;
  title(["Range FFT (Single Chirp) - ", title_tag]);
  xlabel("Range (m)");
  ylabel("Amplitude");

  %% ----------------------
  % Range Doppler Map (2D FFT)
  %% ----------------------
  sig_fft2 = fft2(Mix_windowed, Nr, Nd);
  sig_fft2 = sig_fft2(1:Nr/2, :);
  sig_fft2 = fftshift(sig_fft2, 2);
  RDM = 20*log10(abs(sig_fft2) + 1e-6);

  % Peak detection
  [maxVal, idx] = max(RDM(:));
  [r_idx, d_idx] = ind2sub(size(RDM), idx);

  detected_range = range_axis(r_idx);
  detected_velocity = doppler_axis(d_idx);

  printf("\n==============================\n");
  printf("CASE: %s\n", title_tag);
  printf("Detected Range    = %.2f m\n", detected_range);
  printf("Detected Velocity = %.2f m/s\n", detected_velocity);
  printf("==============================\n");

  % Plot RDM
  figure;
  imagesc(range_axis, doppler_axis, RDM');  % transpose for correct axes
  colorbar;
  grid on;
  title(["Range–Doppler Map - ", title_tag]);
  xlabel("Range (m)");
  ylabel("Velocity (m/s)");

endfunction

%% ============================================================
%  PART A: NOISELESS OUTPUT
%% ============================================================
process_and_plot(Mix_clean, Nr, Nd, range_axis, doppler_axis, "No Noise");

%% ============================================================
%  PART B: NOISY OUTPUTS FOR MULTIPLE SNR VALUES
%% ============================================================

SNR_list = [30, 20, 10, 5];   % test values for report

for k = 1:length(SNR_list)

  SNR_dB = SNR_list(k);

  % Calculate noise power for AWGN
  signal_power = mean(Mix_clean.^2);
  noise_power = signal_power / (10^(SNR_dB/10));
  noise = sqrt(noise_power) * randn(size(Mix_clean));

  % Noisy signal
  Mix_noisy = Mix_clean + noise;

  % Process and plot
  process_and_plot(Mix_noisy, Nr, Nd, range_axis, doppler_axis, ...
                   ["AWGN SNR = ", num2str(SNR_dB), " dB"]);

endfor

