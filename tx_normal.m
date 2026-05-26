%% tx_truyen_thong.m
%% ================== THAM SO ==================
fs   = 48000;
f0   = 12000;
f1   = 12800;
Ns   = 660;

s         = 'Xin chao WicomLab ';
data_bits = reshape(dec2bin(s, 8).' - '0', 1, []);
Ld = length(data_bits);

header_bits = [1 1 1 1 1 0 0 1 1 0 1 0 1];
footer_bits = [0 0 0 0 0 1 1 0 0 1 0 1 0];
Lh = length(header_bits);
Lf = length(footer_bits);
Lz = 30;

%% ================== GHEP FRAME ==================
zero_guard = zeros(1, Lz);
frame_bits = [header_bits, zero_guard, data_bits, zero_guard, footer_bits];
Nbits      = length(frame_bits);

%% ================== TAO TIN HIEU FSK DON GIAN ==================
% Truyen thong: moi bit la 1 doan sin doc lap, KHONG lien tuc pha
% => khong CPFSK, chi la FSK thong thuong
t_bit = (0:Ns-1) / fs;

sig_frame = zeros(1, Nbits * Ns);
for k = 1:Nbits
    if frame_bits(k) == 1
        f = f1;
    else
        f = f0;
    end
    sig_frame((k-1)*Ns+1 : k*Ns) = sin(2*pi*f*t_bit);  % phase bat dau tu 0 moi bit
end

%% ================== KHONG CO SILENCE, KHONG CO GUARD ZERO ==================
% Phat thang frame, khong them silence truoc sau
% Guard bit van co nhung la bit 0 => sin(2*pi*f0*t) binh thuong
tx_sig = sig_frame;

%% ================== XUAT ==================
fprintf('=== TX TRUYEN THONG ===\n');
fprintf('Frame: %d bits | %d mau | %.3f giay\n', Nbits, length(tx_sig), length(tx_sig)/fs);

audiowrite('tx_truyen_thong.wav', tx_sig / max(abs(tx_sig)), fs);
save('tx_truyen_thong.mat', 'tx_sig');
sound(tx_sig, fs);

%% ================== KIEM TRA ==================
figure('Name','TX Truyen Thong');
subplot(2,1,1);
plot((0:length(tx_sig)-1)/fs, tx_sig);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Tin hieu TX truyen thong (khong silence, khong CPFSK)'); grid on;

subplot(2,1,2);
spectrogram(tx_sig, 256, 200, 512, fs, 'yaxis');
title('Spectrogram TX truyen thong'); ylim([8 13]);
