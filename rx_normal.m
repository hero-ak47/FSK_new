%% rx_truyen_thong.m
%% ================== THAM SO ==================
fs   = 48000;
f0   = 12000;
f1   = 12800;
Ns   = 660;
Lz   = 30;

header_bits = [1 1 1 1 1 0 0 1 1 0 1 0 1];
footer_bits = [0 0 0 0 0 1 1 0 0 1 0 1 0];
Lh = length(header_bits);
Lf = length(footer_bits);

header_wave = gen_fsk_simple(header_bits, f0, f1, fs, Ns);
footer_wave = gen_fsk_simple(footer_bits, f0, f1, fs, Ns);

%% ================== THU AM ==================
recObj = audiorecorder(fs, 16, 1);
disp('Dang thu am...');
recordblocking(recObj, 8);
disp('Ket thuc thu.');
rx = getaudiodata(recObj).';

% [rx, fs] = audioread('tx_truyen_thong.wav'); rx = rx.';

%% ================== BANDPASS FILTER ==================
[b, a] = butter(5, [f0-400, f1+400] / (fs/2), 'bandpass');
rx_filt = filtfilt(b, a, rx);

%% ================== PHAT HIEN FRAME ==================
[idx_h, idx_f, peak_h, peak_f] = detect_frame_edges(rx_filt, header_wave, footer_wave, 0.55);
fprintf('Header tai mau %d (peak=%.3f)\n', idx_h, peak_h);
fprintf('Footer tai mau %d (peak=%.3f)\n', idx_f, peak_f);

%% ================== CAT FRAME & DATA ==================
frame_full      = rx_filt(idx_h : idx_f + Lf*Ns - 1);
data_start_samp = (Lh + Lz) * Ns + 1;
data_end_samp   = length(frame_full) - (Lz + Lf) * Ns;
data_frame      = frame_full(data_start_samp : data_end_samp);
fprintf('So bit data uoc tinh: %d\n', floor(length(data_frame)/Ns));

%% ================== GIAI DIEU CHE ==================
bits_rx = fsk_demod_noncoherent(data_frame, fs, Ns, f0, f1);

%% ================== GIAI MA ASCII ==================
n_bytes  = floor(length(bits_rx) / 8);
bits_use = bits_rx(1 : n_bytes*8);
bytes    = reshape(bits_use, 8, []).';
text_out = char(bin2dec(char(bytes + '0')));

fprintf('\n=== KET QUA ===\n');
fprintf('%s\n', text_out);

%% ================================================================
%%  HAM HO TRO
%% ================================================================

function wave = gen_fsk_simple(bits, f0, f1, fs, Ns)
% FSK don gian: moi bit la sin doc lap, phase bat dau tu 0
    t    = (0:Ns-1) / fs;
    wave = zeros(1, length(bits)*Ns);
    for k = 1:length(bits)
        f = f0*(bits(k)==0) + f1*(bits(k)==1);
        wave((k-1)*Ns+1 : k*Ns) = sin(2*pi*f*t);
    end
end

function [idx_h, idx_f, peak_h, peak_f] = detect_frame_edges(rx, header_wave, footer_wave, thr)
    Lh_samp = length(header_wave);

    Rh      = conv(rx, fliplr(header_wave), 'valid');
    Rh_norm = abs(Rh) / max(abs(Rh));
    [peak_h, idx_h] = max(Rh_norm);

    if peak_h < thr
        error('Khong phat hien header (peak=%.3f)', peak_h);
    end

    search_start     = idx_h + Lh_samp;
    Rf2              = conv(rx(search_start:end), fliplr(footer_wave), 'valid');
    Rf2_norm         = abs(Rf2) / max(abs(Rf2));
    [peak_f, pos_f2] = max(Rf2_norm);

    if peak_f < thr
        error('Khong phat hien footer (peak=%.3f)', peak_f);
    end

    idx_f = search_start + pos_f2 - 1;
end

function bits = fsk_demod_noncoherent(data_seg, fs, Ns, f0, f1)
% Noncoherent thuan tuy: so nang luong 2 nhanh BPF, khong IQ, khong xoay pha
    bw = 300;
    [b0, a0] = butter(4, [(f0-bw) (f0+bw)]/(fs/2), 'bandpass');
    [b1, a1] = butter(4, [(f1-bw) (f1+bw)]/(fs/2), 'bandpass');
    seg0  = filtfilt(b0, a0, data_seg);
    seg1  = filtfilt(b1, a1, data_seg);
    Nbits = floor(length(data_seg) / Ns);
    bits  = zeros(1, Nbits);
    for k = 1:Nbits
        idx     = (k-1)*Ns + 1 : k*Ns;
        bits(k) = sum(seg1(idx).^2) > sum(seg0(idx).^2);
    end
end
