%% ================== THAM SO (khop voi TX) ==================
fs   = 48000;
f0   = 11000;
f1   = 11800;
Ns   = 660;
Lz   = 30;

header_bits = [1 1 1 1 1 0 0 1 1 0 1 0 1];
footer_bits = [0 0 0 0 0 1 1 0 0 1 0 1 0];
Lh = length(header_bits);
Lf = length(footer_bits);

% Tạo sóng tham chiếu header/footer (GIONG HET BEN PHAT - CPFSK)
header_wave = gen_cpfsk(header_bits, f0, f1, fs, Ns);
footer_wave = gen_cpfsk(footer_bits, f0, f1, fs, Ns);

%% ================== THU AM ==================
recObj = audiorecorder(fs, 16, 1);
disp('Dang thu am...');
recordblocking(recObj, 8);
disp('Ket thuc thu.');
rx = getaudiodata(recObj).';

% --- Hoac load file de test ---
% [rx, fs] = audioread('test01.wav');
% rx = rx.';

%% ================== BANDPASS FILTER ==================
% Loc bang thong bao quanh [f0-400, f1+400]
[b, a] = butter(5, [f0-400, f1+400] / (fs/2), 'bandpass');
rx_filt = filtfilt(b, a, rx);

%% ================== PHAT HIEN FRAME (MATCHED FILTER) ==================
[idx_h, idx_f, peak_h, peak_f] = detect_frame_edges(rx_filt, header_wave, footer_wave, 0.55);

fprintf('Header tai mau %d (peak=%.3f)\n', idx_h, peak_h);
fprintf('Footer tai mau %d (peak=%.3f)\n', idx_f, peak_f);

%% ================== CAT FRAME ==================
% idx_h: MAU DAU TIEN cua header trong rx
% idx_f: MAU DAU TIEN cua footer trong rx
frame_start = idx_h;
frame_end   = idx_f + Lf*Ns - 1;
frame_full  = rx_filt(frame_start : frame_end);

fprintf('Do dai frame: %d mau = %d bits\n', length(frame_full), floor(length(frame_full)/Ns));

%% ================== CAT PHAN DATA ==================
% Cau truc frame: header(Lh) | guard(Lz) | data(Ld) | guard(Lz) | footer(Lf)
% => data bat dau sau Lh+Lz bits, ket thuc truoc Lz+Lf bits
data_start_samp = (Lh + Lz) * Ns + 1;
data_end_samp   = length(frame_full) - (Lz + Lf) * Ns;

if data_end_samp <= data_start_samp
    error('Frame qua ngan, khong cat duoc phan data.');
end

data_frame = frame_full(data_start_samp : data_end_samp);
Ld_est = floor(length(data_frame) / Ns);
fprintf('So bit data uoc tinh: %d\n', Ld_est);

%% ================== GIAI DIEU CHE FSK (NONCOHERENT) ==================
bits_rx = fsk_demod(data_frame, fs, Ns, f0, f1);

%% ================== GIAI MA ASCII ==================
n_bytes  = floor(length(bits_rx) / 8);
bits_use = bits_rx(1 : n_bytes * 8);
bytes    = reshape(bits_use, 8, []).';
dec_val  = bin2dec(char(bytes + '0'));
text_out = char(dec_val);

fprintf('\n=== KET QUA ===\n');
fprintf('%s\n', text_out);

%% ================================================================
%%  HAM HO TRO
%% ================================================================

function wave = gen_cpfsk(bits, f0, f1, fs, Ns)
% Tao song CPFSK - GIONG HET BEN PHAT
    freq_per_bit = (bits == 1)*f1 + (bits == 0)*f0;
    freq_vec     = repelem(freq_per_bit, Ns);
    delta_phase  = 2 * pi * freq_vec / fs;
    phase        = cumsum(delta_phase);
    wave         = sin(phase);
end

function [idx_h, idx_f, peak_h, peak_f] = detect_frame_edges(rx, header_wave, footer_wave, thr)
% Tim vi tri header va footer bang cross-correlation
% Tra ve: chi so MAU DAU TIEN trong rx

    % Cross-correlation = matched filter
    % conv(rx, fliplr(ref)) hieu qua hon xcorr vi khong can tinh 2 phia
    Lh_samp = length(header_wave);
    Lf_samp = length(footer_wave);

    Rh = conv(rx, fliplr(header_wave), 'valid');
    Rf = conv(rx, fliplr(footer_wave), 'valid');

    Rh_norm = abs(Rh) / max(abs(Rh));
    Rf_norm = abs(Rf) / max(abs(Rf));

    [peak_h, pos_h] = max(Rh_norm);
    [peak_f, pos_f] = max(Rf_norm);

    % pos_h la vi tri KET THUC cua header (conv 'valid' offset)
    % => vi tri DAU TIEN = pos_h
    idx_h = pos_h;
    idx_f = pos_f + Lh_samp;  % footer chi co the xuat hien sau header

    % Kiem tra nguong tin cay
    if peak_h < thr
        error('Khong phat hien header (peak=%.3f < %.3f)', peak_h, thr);
    end
    if peak_f < thr
        error('Khong phat hien footer (peak=%.3f < %.3f)', peak_f, thr);
    end

    % Tim footer chi trong phan sau header de tranh false positive
    search_start = idx_h + Lh_samp;
    Rf2 = conv(rx(search_start:end), fliplr(footer_wave), 'valid');
    Rf2_norm = abs(Rf2) / max(abs(Rf2));
    [peak_f, pos_f2] = max(Rf2_norm);

    if peak_f < thr
        error('Khong phat hien footer sau header (peak=%.3f)', peak_f);
    end

    idx_f    = search_start + pos_f2 - 1;
    peak_f   = peak_f;

    if idx_f <= idx_h + Lh_samp
        error('Footer xuat hien qua som so voi header');
    end
end

function bits = fsk_demod(data_seg, fs, Ns, f0, f1)

    bw = 300;
    [b0, a0] = butter(4, [(f0-bw) (f0+bw)] / (fs/2), 'bandpass');
    [b1, a1] = butter(4, [(f1-bw) (f1+bw)] / (fs/2), 'bandpass');

    seg0 = filtfilt(b0, a0, data_seg);
    seg1 = filtfilt(b1, a1, data_seg);

    Nbits = floor(length(data_seg) / Ns);
    bits  = zeros(1, Nbits);

    %% ============================================================
    %  PHUONG PHAP 1: NONCOHERENT - SO SANH NANG LUONG (binh phuong)
    %  Uu: don gian, khong can biet phase, robust voi nhieu pha
    %  Nhuoc: khong dung toi thong tin pha
    % -------------------------------------------------------------
    % for k = 1:Nbits
    %     idx = (k-1)*Ns + 1 : k*Ns;
    %     E0 = sum(seg0(idx).^2);
    %     E1 = sum(seg1(idx).^2);
    %     bits(k) = (E1 > E0);
    % end

    %% ============================================================
    %  PHUONG PHAP 2: COHERENT - XOAY PHA (IQ demod)
    %  Uu: tan dung thong tin pha, ly thuyet tot hon trong AWGN
    %  Nhuoc: nhay cam voi lech pha kenh, can uoc luong phase offset
    % -------------------------------------------------------------
    n_all = (0 : length(data_seg)-1);

    % IQ demod nhanh f0
    I0 = data_seg .* cos(2*pi*f0*n_all/fs);
    Q0 = data_seg .* sin(2*pi*f0*n_all/fs);
    [lp_b, lp_a] = butter(4, (f1-f0)/(fs/2), 'low');   % LPF cat tai (f1-f0)/2
    I0 = filtfilt(lp_b, lp_a, I0);
    Q0 = filtfilt(lp_b, lp_a, Q0);
    z0 = I0 + 1j*Q0;   % baseband phuc nhanh f0

    % IQ demod nhanh f1
    I1 = data_seg .* cos(2*pi*f1*n_all/fs);
    Q1 = data_seg .* sin(2*pi*f1*n_all/fs);
    I1 = filtfilt(lp_b, lp_a, I1);
    Q1 = filtfilt(lp_b, lp_a, Q1);
    z1 = I1 + 1j*Q1;   % baseband phuc nhanh f1

    % Uoc luong phase offset toan cuc tu toan bo frame
    % (gia su phase offset gan nhu khong doi trong 1 frame)
    phi0_est = angle(sum(z0));   % pha trung binh nhanh f0
    phi1_est = angle(sum(z1));   % pha trung binh nhanh f1

% Xoay pha ve 0
z0_rot = z0 .* exp(-1j * phi0_est);
z1_rot = z1 .* exp(-1j * phi1_est);

% --- Lay mau dai dien tai giua moi symbol ---
mid     = round(Ns * 0.6);
idx_mid = (0:Nbits-1)*Ns + mid;

z_diff         = z1_rot - z0_rot;       % tinh hieu truoc
z_diff_samples = z_diff(idx_mid);       % roi moi lay mau

scatterplot(z_diff_samples);
title('Chom sao FSK: z1\_rot - z0\_rot (mau giua symbol)');

    % Quyet dinh theo phan thuc sau xoay pha
    for k = 1:Nbits
        idx = (k-1)*Ns + 1 : k*Ns;
        % Dung tong phan thuc lam metric (tuong duong MF)
        m0 = sum(real(z0_rot(idx)));
        m1 = sum(real(z1_rot(idx)));
        bits(k) = (m1 > m0);

    end

    
    

    

end
