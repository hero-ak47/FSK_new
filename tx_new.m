%% ================== THAM SO ==================
fs   = 48000;
f0   = 11000;
f1   = 11800;
Ns   = 660;

s         = 'Xin chao WicomLab';
data_bits = reshape(dec2bin(s, 8).' - '0', 1, []);
Ld = length(data_bits);

header_bits = [1 1 1 1 1 0 0 1 1 0 1 0 1];
footer_bits = [0 0 0 0 0 1 1 0 0 1 0 1 0];
Lh = length(header_bits);
Lf = length(footer_bits);
Lz = 30;

silence_time = 0.5;
silence_samp = zeros(1, round(silence_time * fs));

%% ================== GHEP FRAME ==================
zero_guard  = zeros(1, Lz);
frame_bits  = [header_bits, zero_guard, data_bits, zero_guard, footer_bits];
Nbits = length(frame_bits);

%% ================== PHASE ACCUMULATOR (CPFSK) ==================
% Mỗi bit chiếm đúng Ns mẫu, tần số chọn theo bit
freq_per_bit = (frame_bits == 1) * f1 + (frame_bits == 0) * f0;
freq_vec     = repelem(freq_per_bit, Ns);   % [1 x Nbits*Ns]
Ntotal       = length(freq_vec);

% Accumulate phase liên tục (CPFSK chuẩn)
delta_phase  = 2 * pi * freq_vec / fs;
phase        = cumsum(delta_phase);          % không mod ở đây để giữ continuity

%% ================== TẠO TÍN HIỆU ==================
sig_frame = sin(phase);

%% ================== CHÈN IM LẶNG CHO GUARD BAND ==================
% Chèn zero ở vùng guard TRONG MIỀN BIÊN DO (không phải miền phase)
% Điều này giữ nguyên phase accumulator nhưng tắt biên độ
guard_start_samp = Lh * Ns + 1;
guard_end_samp   = (Lh + Lz) * Ns;
data_start_samp  = guard_end_samp + 1;
data_end_samp    = (Lh + Lz + Ld) * Ns;
guard2_start     = data_end_samp + 1;
guard2_end       = (Lh + Lz + Ld + Lz) * Ns;

sig_frame(guard_start_samp : guard_end_samp) = 0;
sig_frame(guard2_start     : guard2_end)     = 0;

%% ================== GHEP VA XUAT ==================
tx_sig = [silence_samp, sig_frame, silence_samp];

fprintf('Frame: %d bits | %d mau\n', Nbits, length(sig_frame));
fprintf('Data : %d bits | van ban: "%s"\n', Ld, s);
fprintf('Total audio: %.3f giay\n', length(tx_sig)/fs);
disp('Data bits:'); disp(data_bits);

audiowrite('test01.wav', tx_sig / max(abs(tx_sig)), fs);  % normalize truoc khi ghi
save('test_data.mat', 'tx_sig');

sound(tx_sig, fs);

%% ================== KIEM TRA TRUC QUAN ==================
figure;
subplot(3,1,1);
plot((0:length(tx_sig)-1)/fs, tx_sig);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Tin hieu TX toan bo'); grid on;

subplot(3,1,2);
t_frame = (0:length(sig_frame)-1)/fs;
plot(t_frame, sig_frame);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Frame FSK (khong co silence)'); grid on;
xline(guard_start_samp/fs, 'r--', 'Guard1 start');
xline(data_start_samp/fs,  'g--', 'Data start');
xline(guard2_start/fs,     'r--', 'Guard2 start');

subplot(3,1,3);
% Spectrogram de kiem tra
spectrogram(sig_frame, 256, 200, 512, fs, 'yaxis');
title('Spectrogram Frame');
ylim([8 13]);  % kHz
