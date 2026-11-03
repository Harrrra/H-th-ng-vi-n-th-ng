%% MÔ PHỎNG HỆ THỐNG OFDM (CP-OFDM kiểu 5G NR thu nhỏ)
%  Không cần toolbox nào (chỉ MATLAB cơ bản, R2016b trở lên).
%  Lưu file ở dạng UTF-8 để hiển thị đúng tiếng Việt.
%
%  Chuỗi xử lý:
%  bit -> QAM -> chèn pilot -> IFFT -> thêm CP -> kênh (đa đường + AWGN)
%      -> bỏ CP -> FFT -> ước lượng kênh (LS) -> cân bằng ZF -> bù pha -> giải điều chế -> BER
%
%  Kết quả: 6 hình
%   Hình 1: BER theo SNR (QPSK/16-QAM/64-QAM), kênh AWGN và kênh đa đường Rayleigh
%   Hình 2: Ảnh hưởng độ dài CP lên BER
%   Hình 3: Ước lượng kênh LS thường và LS có làm trơn
%   Hình 4: Chòm sao trước/sau cân bằng
%   Hình 5: Phổ tín hiệu OFDM
%   Hình 6: CCDF của PAPR, có và không có clipping

clear; close all; clc; rng(1);

%% ===== 1. THAM SỐ HỆ THỐNG =====
cfg.Nfft      = 64;                              % kích thước FFT
cfg.Ncp       = 16;                              % độ dài CP (1/4 Nfft)
cfg.scs       = 30e3;                            % khoảng cách sóng mang con (Hz), giống numerology mu=1 của 5G NR
cfg.used      = [-26:-1, 1:26];                  % 52 sóng mang dùng (bỏ DC và băng bảo vệ)
cfg.pilotIdx  = [-21 -7 7 21];                   % vị trí pilot (để bám pha)
cfg.pilotVal  = [1 1 1 -1].';                    % giá trị pilot (BPSK)
cfg.dataIdx   = setdiff(cfg.used, cfg.pilotIdx); % 48 sóng mang dữ liệu
cfg.pre       = 1 - 2*randi([0 1], numel(cfg.used), 1); % preamble BPSK đã biết, dùng ước lượng kênh
cfg.nSym      = 20;                              % số ký hiệu OFDM dữ liệu mỗi khung
cfg.M         = 16;                              % bậc điều chế
cfg.snr       = 20;                              % Es/N0 trên mỗi sóng mang (dB)
cfg.chan      = 'rayleigh';                      % 'awgn' hoặc 'rayleigh'
cfg.L         = 6;                               % số tap kênh đa đường
cfg.decay     = 2;                               % hằng số suy giảm theo hàm mũ của PDP
cfg.smooth    = 1;                               % cửa sổ làm trơn ước lượng kênh (1 = không làm trơn)
cfg.frames    = 300;                             % số khung tối đa mỗi điểm SNR
cfg.minErr    = 400;                             % dừng sớm khi đủ số bit lỗi

fprintf('Fs = %.2f MHz | băng thông dùng ~ %.2f MHz | CP = %.2f us\n', ...
    cfg.Nfft*cfg.scs/1e6, numel(cfg.used)*cfg.scs/1e6, cfg.Ncp/(cfg.Nfft*cfg.scs)*1e6);

Q = @(x) 0.5*erfc(x/sqrt(2));
mods = [4 16 64];
col  = {'b','r','k'};

%% ===== 2. HÌNH 1: BER THEO SNR =====
snrList = 0:2:30;
figure('Name','BER theo SNR','Position',[80 80 1100 450]);
for c = 1:2
    subplot(1,2,c); hold on; grid on;
    if c == 1, cfg.chan = 'awgn'; ttl = 'Kênh AWGN';
    else,      cfg.chan = 'rayleigh'; ttl = 'Kênh đa đường Rayleigh (LS + ZF)';
    end
    for m = 1:numel(mods)
        cfg.M = mods(m); k = log2(cfg.M);
        ber = sweep_snr(cfg, snrList);
        ber(ber == 0) = NaN;
        semilogy(snrList, ber, [col{m} 'o'], 'MarkerFaceColor', col{m});
        if c == 1   % đường lý thuyết cho AWGN
            snrLin = 10.^(snrList/10);
            th = (4/k)*(1-1/sqrt(cfg.M))*Q(sqrt(3*snrLin/(cfg.M-1)));
            semilogy(snrList, th, [col{m} '-']);
        end
        fprintf('%s | %2d-QAM xong\n', ttl, cfg.M);
    end
    set(gca,'YScale','log'); ylim([1e-5 1]); xlabel('Es/N0 (dB)'); ylabel('BER'); title(ttl);
    if c == 1
        legend('QPSK mô phỏng','QPSK lý thuyết','16-QAM mô phỏng','16-QAM lý thuyết', ...
               '64-QAM mô phỏng','64-QAM lý thuyết','Location','southwest');
    else
        legend('QPSK','16-QAM','64-QAM','Location','southwest');
    end
end

%% ===== 3. HÌNH 2: ẢNH HƯỞNG ĐỘ DÀI CP =====
cfg2 = cfg; cfg2.chan = 'rayleigh'; cfg2.M = 16; cfg2.snr = 30;
cfg2.L = 12; cfg2.decay = 5;                     % kênh có trải trễ dài (12 tap)
cpList = [0 2 4 8 12 16];
berCP = zeros(size(cpList));
for i = 1:numel(cpList)
    cfg2.Ncp = cpList(i);
    berCP(i) = run_link(cfg2);
end
figure('Name','Ảnh hưởng CP'); semilogy(cpList, berCP, 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
grid on; xlabel('Độ dài CP (mẫu)'); ylabel('BER');
title('16-QAM, SNR = 30 dB, kênh trải trễ 12 tap');
xline(cfg2.L-1, '--k', 'L-1 = 11');               % CP cần >= L-1 để hết ISI

%% ===== 4. HÌNH 3: LS THƯỜNG vs LS CÓ LÀM TRƠN =====
cfg3 = cfg; cfg3.chan = 'rayleigh'; cfg3.M = 16; cfg3.L = 6; cfg3.decay = 2;
snr3 = 0:3:30;
cfg3.smooth = 1; b1 = sweep_snr(cfg3, snr3);
cfg3.smooth = 3; b2 = sweep_snr(cfg3, snr3);
figure('Name','Ước lượng kênh'); semilogy(snr3, b1, 'bo-', snr3, b2, 'rs-', 'LineWidth', 1.5);
grid on; xlabel('Es/N0 (dB)'); ylabel('BER'); legend('LS','LS làm trơn (cửa sổ 3)','Location','southwest');
title('So sánh ước lượng kênh, 16-QAM');

%% ===== 5. HÌNH 4: CHÒM SAO =====
cfg4 = cfg; cfg4.chan = 'rayleigh'; cfg4.M = 16; cfg4.snr = 25; cfg4.nSym = 100; cfg4.frames = 1;
[~, info] = run_link(cfg4);
figure('Name','Chòm sao','Position',[100 100 950 420]);
subplot(1,2,1); plot(real(info.rawSym), imag(info.rawSym), '.'); axis equal; grid on;
title('Trước cân bằng (sau FFT)'); xlabel('I'); ylabel('Q');
subplot(1,2,2); plot(real(info.eqSym), imag(info.eqSym), '.'); axis equal; grid on;
title('Sau cân bằng kênh'); xlabel('I'); ylabel('Q');

%% ===== 6. HÌNH 5: PHỔ TÍN HIỆU =====
os = 8; Nos = cfg.Nfft*os;
x = gen_ofdm_os(cfg, 16, 400, os);
x = [x(end-cfg.Ncp*os+1:end,:); x];
s = x(:);
seg = 4096; nseg = floor(numel(s)/seg);
w = 0.5 - 0.5*cos(2*pi*(0:seg-1)'/seg);
P = zeros(seg,1);
for i = 1:nseg
    P = P + abs(fft(s((i-1)*seg+1:i*seg) .* w)).^2;
end
P = fftshift(P)/nseg;
f = (-seg/2:seg/2-1)'/seg * Nos * cfg.scs / 1e6;  % MHz
figure('Name','Phổ OFDM'); plot(f, 10*log10(P/max(P)), 'LineWidth', 1);
grid on; xlim([-4 4]); ylim([-80 5]); xlabel('Tần số (MHz)'); ylabel('Mật độ phổ chuẩn hoá (dB)');
title('Phổ tín hiệu OFDM (rò phổ ngoài băng do cửa sổ chữ nhật)');

%% ===== 7. HÌNH 6: PAPR =====
crList = [Inf 6 4 3];                              % tỉ số clipping (dB so với RMS)
figure('Name','PAPR'); hold on; grid on;
leg = cell(1, numel(crList));
for i = 1:numel(crList)
    papr = ofdm_papr(cfg, 16, 20000, 10^(crList(i)/20), 4);
    pdB = sort(10*log10(papr));
    ccdf = (numel(pdB):-1:1)/numel(pdB);
    semilogy(pdB, ccdf, 'LineWidth', 1.5);
    if isinf(crList(i)), leg{i} = 'Không clipping'; else, leg{i} = sprintf('Clipping %g dB', crList(i)); end
end
set(gca,'YScale','log'); ylim([1e-4 1]); xlabel('PAPR_0 (dB)'); ylabel('P(PAPR > PAPR_0)');
legend(leg, 'Location', 'southwest'); title('CCDF của PAPR, 16-QAM, oversampling x4');

disp('Hoàn tất.');

%% =====================================================================
%%                           HÀM CỤC BỘ
%% =====================================================================

function [ber, info] = run_link(cfg)
% Mô phỏng một liên kết OFDM, trả về BER và vài dữ liệu để vẽ chòm sao.
    N = cfg.Nfft; k = log2(cfg.M); nd = numel(cfg.dataIdx);
    sigma2 = 10^(-cfg.snr/10);                    % phương sai nhiễu trên mỗi sóng mang
    bu = mod(cfg.used, N) + 1;
    bd = mod(cfg.dataIdx, N) + 1;
    bp = mod(cfg.pilotIdx, N) + 1;
    errs = 0; tot = 0; info = struct();
    for f = 1:cfg.frames
        % --- Phát ---
        bits = randi([0 1], k*nd*cfg.nSym, 1);
        S = reshape(qam_mod(bits, cfg.M), nd, cfg.nSym);
        X = zeros(N, cfg.nSym + 2);
        X(bu, 1:2)       = repmat(cfg.pre, 1, 2);         % 2 ký hiệu preamble giống nhau
        X(bd, 3:end)     = S;
        X(bp, 3:end)     = repmat(cfg.pilotVal, 1, cfg.nSym);
        x = ifft(X) * sqrt(N);                            % IFFT chuẩn hoá unitary
        x = [x(end-cfg.Ncp+1:end, :); x];                 % thêm CP
        tx = x(:);
        % --- Kênh ---
        h  = make_channel(cfg);
        rx = filter(h, 1, tx);
        rx = rx + sqrt(sigma2/2) * (randn(size(rx)) + 1j*randn(size(rx)));
        % --- Thu ---
        Y = reshape(rx, N + cfg.Ncp, []);
        Y = Y(cfg.Ncp+1:end, :);                          % bỏ CP
        Y = fft(Y) / sqrt(N);
        Hest = mean(Y(bu, 1:2), 2) ./ cfg.pre;            % ước lượng kênh LS, trung bình 2 preamble
        if cfg.smooth > 1
            Hest = movmean(Hest, cfg.smooth);             % làm trơn theo tần số
        end
        Hf = zeros(N, 1); Hf(bu) = Hest;
        Yd = Y(bd, 3:end); Yp = Y(bp, 3:end);
        Xd = Yd ./ Hf(bd);                                % cân bằng ZF một tap trên mỗi sóng mang
        Xp = Yp ./ Hf(bp);
        phi = angle(sum(Xp .* conj(cfg.pilotVal), 1));    % pha chung (CPE) đo từ pilot
        Xd = Xd .* exp(-1j*phi);
        rbits = qam_demod(Xd(:), cfg.M);
        errs = errs + sum(rbits ~= bits);
        tot  = tot + numel(bits);
        if f == cfg.frames || (errs >= cfg.minErr && f >= 20), break; end
    end
    ber = errs / tot;
    info.rawSym = Yd(:);
    info.eqSym  = Xd(:);
end

function ber = sweep_snr(cfg, snrList)
    ber = zeros(size(snrList));
    for i = 1:numel(snrList)
        cfg.snr = snrList(i);
        ber(i) = run_link(cfg);
    end
end

function h = make_channel(cfg)
% Kênh đa đường: PDP suy giảm theo hàm mũ, các tap Rayleigh độc lập, tổng công suất = 1.
    switch cfg.chan
        case 'awgn'
            h = 1;
        case 'rayleigh'
            p = exp(-(0:cfg.L-1)'/cfg.decay); p = p / sum(p);
            h = sqrt(p/2) .* (randn(cfg.L,1) + 1j*randn(cfg.L,1));
    end
end

function x = gen_ofdm_os(cfg, M, nSym, os)
% Sinh tín hiệu OFDM miền thời gian với oversampling (chưa có CP), dùng cho phổ và PAPR.
    Nos = cfg.Nfft * os; k = log2(M); nd = numel(cfg.dataIdx);
    bits = randi([0 1], k*nd*nSym, 1);
    S = reshape(qam_mod(bits, M), nd, nSym);
    X = zeros(Nos, nSym);
    X(mod(cfg.dataIdx, Nos) + 1, :)  = S;
    X(mod(cfg.pilotIdx, Nos) + 1, :) = repmat(cfg.pilotVal, 1, nSym);
    x = ifft(X) * sqrt(Nos);
end

function papr = ofdm_papr(cfg, M, nSym, clipRatio, os)
    x = gen_ofdm_os(cfg, M, nSym, os);
    if isfinite(clipRatio)
        A = clipRatio * sqrt(mean(abs(x(:)).^2));
        mag = abs(x); idx = mag > A;
        x(idx) = x(idx) ./ mag(idx) * A;              % cắt biên độ, giữ pha
    end
    papr = max(abs(x).^2, [], 1) ./ mean(abs(x).^2, 1);
end

function sym = qam_mod(bits, M)
% Điều chế M-QAM vuông, ánh xạ Gray, công suất trung bình = 1.
    k = log2(M); kb = k/2; L = sqrt(M);
    b  = reshape(bits, k, []);
    w  = 2.^(kb-1:-1:0);
    gI = uint32(w * b(1:kb, :));
    gQ = uint32(w * b(kb+1:k, :));
    iI = double(gray2bin(gI)); iQ = double(gray2bin(gQ));
    sym = ((2*iI - (L-1)) + 1j*(2*iQ - (L-1))) / sqrt(2*(M-1)/3);
    sym = sym(:);
end

function bits = qam_demod(sym, M)
% Giải điều chế quyết định cứng (khoảng cách nhỏ nhất).
    k = log2(M); kb = k/2; L = sqrt(M);
    r  = sym(:).' * sqrt(2*(M-1)/3);
    iI = min(max(round((real(r) + (L-1))/2), 0), L-1);
    iQ = min(max(round((imag(r) + (L-1))/2), 0), L-1);
    gI = bitxor(uint32(iI), bitshift(uint32(iI), -1));
    gQ = bitxor(uint32(iQ), bitshift(uint32(iQ), -1));
    pos = repmat((kb:-1:1).', 1, numel(gI));
    bI = bitget(repmat(gI, kb, 1), pos);
    bQ = bitget(repmat(gQ, kb, 1), pos);
    bits = double([bI; bQ]);
    bits = bits(:);
end

function b = gray2bin(g)
    b = g; s = bitshift(g, -1);
    while any(s(:))
        b = bitxor(b, s);
        s = bitshift(s, -1);
    end
end
