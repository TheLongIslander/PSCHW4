clear; clc; close all;

%% 1. Set Random Seed and Generate Bit Sequence
RUID = 208001821;
rng(RUID);
bb = randi([0, 1], 1, 1000);

%% 2. Generate Baseband Signals
T = 2;        
A = 1;       
Ts = 0.02;    
fs = 1/Ts; 

t = 0:Ts:T-Ts;

% Define square pulse p(t)
p_t = A * ones(size(t));

% Define raised cosine pulse p_s(t)
r = 5; 
p_s_t = sinc(t/T) .* cos(pi*r*t/T) ./ (1 - (2*r*t/T).^2);
p_s_t(abs(2*r*t/T) == 1) = 0; 
p_s_t = p_s_t / max(abs(p_s_t));


% Generate s(t) and s_s(t)
s = [];
s_s = [];
for i = 1:10 
    if bb(i) == 1
        s = [s p_t];
        s_s = [s_s p_s_t];
    else
        s = [s -p_t]; 
        s_s = [s_s -p_s_t];
    end
end

% Time vector for first 10 bits
t_full = 0:Ts:(10*T-Ts);

% Plot baseband signals
figure;
subplot(2,1,1);
plot(t_full, s, 'b');
title('Baseband Signal using Square Pulse');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t_full, s_s, 'r');
title('Baseband Signal using Raised Cosine Pulse');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% 3. Up-conversion
fc = 5;

% Perform up-conversion
u = s .* cos(2*pi*fc*t_full);
u_s = s_s .* cos(2*pi*fc*t_full);

% Normalize u_s to avoid very small values
if max(abs(u_s)) > 0
    u_s = u_s / max(abs(u_s));
else
    warning('u_s(t) is all zeros. Check baseband signal generation.');
end

% Compute FFT
U = fftshift(abs(fft(u)));
U_s = fftshift(abs(fft(u_s)));

f = linspace(-fs/2, fs/2, length(U));

% Plot Energy Spectral Density (ESD)
figure;
subplot(2,1,1);
plot(f, 10*log10(U.^2 + eps));
title('Energy Spectral Density of u(t)');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
grid on;
xlim([-fc-5 fc+5]);

subplot(2,1,2);
plot(f, 10*log10(U_s.^2 + eps));
title('Energy Spectral Density of u_s(t)');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
grid on;
xlim([-fc-5 fc+5]);


%% 4. Down-conversion and Filtering
% Multiply by cos(2*pi*fc*t) to down-convert
d = u .* cos(2*pi*fc*t_full);
d_s = u_s .* cos(2*pi*fc*t_full);

% Design LPF - Moving Average Filter
L1 = 2; 
L2 = 10;

h1 = ones(1, L1) / L1; 
h2 = ones(1, L2) / L2;

d_f1 = filter(h1, 1, d);
d_s_f1 = filter(h1, 1, d_s);

d_f2 = filter(h2, 1, d);
d_s_f2 = filter(h2, 1, d_s);

% Plot down-converted signals
figure;
subplot(2,1,1);
plot(t_full, d_f1);
title('Down-converted Signal with L=2');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t_full, d_s_f1);
title('Down-converted Raised Cosine Signal with L=2');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

figure;
subplot(2,1,1);
plot(t_full, d_f2);
title('Down-converted Signal with L=10');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t_full, d_s_f2);
title('Down-converted Raised Cosine Signal with L=10');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% 5. Detection and Decoding
% Correct sampling indices (mid-bit sampling)
num_bits = 10;
mid_sample_indices = round((T/2: T : num_bits * T) / Ts);

% Ensure indices are within range
mid_sample_indices = mid_sample_indices(mid_sample_indices <= length(d_f1));

% Simple threshold detection
decoded_bits_L2 = d_f1(mid_sample_indices) > 0; 
decoded_bits_L10 = d_f2(mid_sample_indices) > 0;

% Compute bit errors
bit_errors_L2 = sum(decoded_bits_L2 ~= bb(1:length(decoded_bits_L2)));
bit_errors_L10 = sum(decoded_bits_L10 ~= bb(1:length(decoded_bits_L10)));

% Display Results
disp('Decoded Bits with L=2: ');
disp(decoded_bits_L2);
disp(['Number of Bit Errors with L=2: ', num2str(bit_errors_L2)]);

disp('Decoded Bits with L=10: ');
disp(decoded_bits_L10);
disp(['Number of Bit Errors with L=10: ', num2str(bit_errors_L10)]);

%% 7. Save Figures
% Create a table with bit errors
bit_error_table = table(["L=2"; "L=10"], [bit_errors_L2; bit_errors_L10], ...
    'VariableNames', {'Filter_Length', 'Bit_Errors'});

disp('Bit Error Table:');
disp(bit_error_table);

%% Save Report
writetable(bit_error_table, 'bit_error_report.csv');
