%% MODEM OFDM QUA ÂM THANH (loa -> mic) - dùng cùng nguyên lý với mô phỏng
%  Chỉ cần MATLAB cơ bản (audioplayer, audiorecorder, audiowrite).
%  Lưu file ở dạng UTF-8.
%
%  Tín hiệu OFDM được tạo thẳng ở băng thông âm thanh (thực, quanh fc = 12 kHz)
%  bằng cách đặt sóng mang con vào các bin FFT và lấy đối xứng liên hợp (Hermite),
%  nên không cần bộ trộn tần hay lọc nội suy.
%
%  Cấu trúc khung:
%  [im lặng] [chirp đồng bộ] [guard] [2 ký hiệu preamble] [nSym ký hiệu dữ liệu] [im lặng]
%
%  Máy thu: tương quan với chirp -> tìm điểm bắt đầu khung -> bỏ CP -> FFT
%           -> ước lượng lệch tần (từ 2 preamble giống nhau) -> ước lượng kênh LS
%           -> cân bằng ZF -> bám pha bằng pilot -> giải điều chế -> BER
%
%  CÁC CHẾ ĐỘ (p.mode):
%   'sim'   : kiểm tra code, kênh mô phỏng (đa đường + nhiễu), không cần loa/mic
%   'audio' : phát qua loa và thu bằng mic của cùng một máy tính
%   'tx'    : chỉ phát (máy phát khi dùng 2 thiết bị)
%   'file'  : giải mã file ghi âm p.rxFile (ghi từ máy/điện thoại khác)
%  Cả hai máy phải dùng cùng tham số và cùng p.seed.

clear; close all; clc;

%% ===== THAM SỐ =====
p.mode   = 'audio';
p.rxFile = 'tx.wav';
p.fs     = 48000;                 % tần số lấy mẫu âm thanh
p.N      = 512;                   % FFT -> khoảng cách sóng mang con = fs/N = 93.75 Hz
p.Ncp    = 128;                   % CP = 1/4 N (2.67 ms)
p.fc     = 12000;                 % tần số trung tâm (Hz)
p.M      = 16;                    % 4, 16 hoặc 64
p.nSym   = 100;                   % số ký hiệu dữ liệu
p.used   = [-26:-1, 1:26];
p.pilotIdx = [-21 -7 7 21];
p.pilotVal = [1 1 1 -1].';
p.dataIdx  = setdiff(p.used, p.pilotIdx);
p.scs    = p.fs / p.N;
p.bin0   = round(p.fc / p.scs);   % bin ứng với fc (phải là số nguyên)
p.lead   = 0.3; p.guard = 0.02; p.tail = 0.3;   % giây
p.chirpT = 0.1; p.chirpF = [5000 19000];        % chirp đồng bộ
p.amp    = 0.7;                   % biên độ đỉnh (tránh méo loa)
p.seed   = 2024;                  % hạt giống cho bit và preamble (hai đầu phải giống nhau)
% chỉ dùng cho chế độ 'sim'
p.simSNR   = 25;                  % SNR toàn băng (dB)
p.simSeed  = 2025;                 % seed riêng cho nhiễu; đổi thành 2025 ở lần chạy tiếp theo
p.simDelay = 1500;                % trễ (mẫu)
p.simH     = zeros(61,1); p.simH([1 21 61]) = [1 0.4 0.2];   % đa đường: 0, 0.42 ms, 1.27 ms

fprintf('Băng chiếm: %.2f - %.2f kHz | tốc độ bit ~ %.0f bit/s | CP = %.2f ms\n', ...
    (p.bin0+min(p.used))*p.scs/1e3, (p.bin0+max(p.used))*p.scs/1e3, ...
    numel(p.dataIdx)*log2(p.M)/((p.N+p.Ncp)/p.fs), p.Ncp/p.fs*1e3);

%% ===== PHÁT =====
[tx, ref] = build_frame(p);
audiowrite('tx.wav', tx, p.fs);
fprintf('Thời lượng phát: %.2f s (đã lưu tx.wav)\n', numel(tx)/p.fs);

%% ===== KÊNH / THU =====
switch p.mode
    case 'sim'
        rx = [zeros(p.simDelay,1); filter(p.simH, 1, tx)];
        sp = mean(tx(tx ~= 0).^2);

        % Cố định seed riêng cho nhiễu để có thể tái lập từng lần chạy
        rng(p.simSeed, 'twister');
        noise = sqrt(sp*10^(-p.simSNR/10)) * randn(size(rx));
        rx = rx + noise;

    case 'audio'
        disp('Đang phát và thu... (để loa và mic ở vị trí cố định, giữ yên lặng)');
        rx = audio_loopback(tx, p);
    case 'tx'
        disp('Đang phát...'); playblocking(audioplayer(tx, p.fs)); disp('Xong.'); return;
    case 'file'
        [rx, fsr] = audioread(p.rxFile); rx = rx(:,1);
        if fsr ~= p.fs, error('File thu có fs = %d, cần %d. Hãy ghi hoặc đổi sang %d Hz.', fsr, p.fs, p.fs); end
    otherwise
        error('p.mode không hợp lệ');
end

%% ===== GIẢI MÃ =====
res = decode_frame(rx, p, ref);

fprintf('\n===== KẾT QUẢ =====\n');
fprintf('Chất lượng đồng bộ (đỉnh/trung vị): %.1f\n', res.syncQ);
fprintf('Lệch tần ước lượng: %.2f Hz\n', res.cfoHz);
fprintf('SNR ước lượng từ preamble: %.1f dB\n', res.snrEst);
fprintf('EVM: %.1f dB\n', res.evmdB);
fprintf('BER = %.3e  (%d lỗi / %d bit)\n', res.ber, res.nErr, res.nBit);

%% ===== HÌNH =====
figure('Name','Modem OFDM âm thanh','Position',[60 60 1150 800]);

subplot(2,2,1);
t = ((1:numel(res.corr)) - res.i0) / p.fs * 1e3;
plot(t, res.corr/max(res.corr)); grid on; xlim([-20 20]);
xlabel('Lệch so với đỉnh (ms)'); ylabel('Tương quan chuẩn hoá'); title('Đồng bộ: tương quan với chirp');

subplot(2,2,2);
seg = 2048; seg0 = res.s0; nseg = floor((numel(res.rxFrame))/seg);
w = 0.5 - 0.5*cos(2*pi*(0:seg-1)'/seg); P = zeros(seg,1);
for i = 1:nseg, P = P + abs(fft(res.rxFrame((i-1)*seg+1:i*seg).*w)).^2; end
fk = (0:seg/2-1)'/seg*p.fs/1e3;
plot(fk, 10*log10(P(1:seg/2)/max(P(1:seg/2)))); grid on; ylim([-80 5]);
xlabel('Tần số (kHz)'); ylabel('dB'); title('Phổ tín hiệu thu (đoạn OFDM)');

subplot(2,2,3);
fH = (p.bin0 + p.used)*p.scs/1e3;
plot(fH, 20*log10(abs(res.Hest)), '.-'); grid on;
xlabel('Tần số (kHz)'); ylabel('|H| (dB, tương đối)'); title('Đáp ứng tần số kênh ước lượng (loa + phòng + mic)');

subplot(2,2,4);
plot(real(res.Xd), imag(res.Xd), '.'); axis equal; grid on;
xlabel('I'); ylabel('Q'); title(sprintf('Chòm sao sau cân bằng, BER = %.2e', res.ber));

%% =====================================================================
%%                           HÀM CỤC BỘ
%% =====================================================================

function rx = audio_loopback(tx, p)
    rec = audiorecorder(p.fs, 16, 1);
    ply = audioplayer(tx, p.fs);
    record(rec);
    pause(0.3);
    play(ply);
    pause(numel(tx)/p.fs + 0.5);
    stop(rec);
    rx = getaudiodata(rec);
    rx = rx(:);
end

function [y, ref] = build_frame(p)
    N = p.N; k = log2(p.M); nd = numel(p.dataIdx);
    st = rng; rng(p.seed);
    ref.pre  = 1 - 2*randi([0 1], numel(p.used), 1);          % preamble BPSK đã biết
    ref.bits = randi([0 1], k*nd*p.nSym, 1);                  % bit dữ liệu (cả hai đầu biết để đo BER)
    rng(st);
    S = reshape(qam_mod(ref.bits, p.M), nd, p.nSym);
    nTot = p.nSym + 2;
    bu = p.bin0 + p.used; bd = p.bin0 + p.dataIdx; bp = p.bin0 + p.pilotIdx;
    Xf = zeros(N, nTot);
    Xf(bu+1, 1:2)   = repmat(ref.pre, 1, 2);                  % 2 preamble giống nhau
    Xf(bd+1, 3:end) = S;
    Xf(bp+1, 3:end) = repmat(p.pilotVal, 1, p.nSym);
    Xf(N-bu+1, :)   = conj(Xf(bu+1, :));                      % đối xứng Hermite -> tín hiệu thực
    x = real(ifft(Xf));
    x = [x(end-p.Ncp+1:end, :); x];                           % thêm CP
    ofdm = x(:);
    ofdm = ofdm / max(abs(ofdm)) * p.amp;
    % chirp đồng bộ
    Lc = round(p.chirpT*p.fs); t = (0:Lc-1)'/p.fs;
    c = sin(2*pi*(p.chirpF(1)*t + (p.chirpF(2)-p.chirpF(1))/(2*p.chirpT)*t.^2));
    c = c .* (0.5 - 0.5*cos(2*pi*(0:Lc-1)'/(Lc-1))) * p.amp;
    ref.chirp = c;
    y = [zeros(round(p.lead*p.fs),1); c; zeros(round(p.guard*p.fs),1); ofdm; zeros(round(p.tail*p.fs),1)];
end

function res = decode_frame(rx, p, ref)
    N = p.N; Ns = N + p.Ncp; k = log2(p.M); nTot = p.nSym + 2;
    rx = rx(:);
    % --- đồng bộ ---
    c = xcorr_fft(rx, ref.chirp);
    [pk, i0] = max(c);
    res.corr = c; res.i0 = i0;
    res.syncQ = pk / (median(abs(c)) + eps);
    if res.syncQ < 8
        warning('Không tìm thấy chirp rõ ràng (chất lượng đồng bộ thấp). Tăng âm lượng hoặc giảm tiếng ồn.');
    end
    s0 = i0 + numel(ref.chirp) + round(p.guard*p.fs);    % mẫu đầu tiên của preamble 1 (đầu CP)
    back = 16;                                           % lùi vào trong CP để an toàn về thời gian
    last = s0 + (nTot-1)*Ns + p.Ncp - back + N - 1;
    if last > numel(rx)
        error('Bản ghi quá ngắn hoặc đồng bộ sai (cần %d mẫu, có %d).', last, numel(rx));
    end
    res.s0 = s0;
    res.rxFrame = rx(s0:last);
    % --- FFT từng ký hiệu ---
    bu = p.bin0 + p.used;
    Y = zeros(numel(p.used), nTot);
    for m = 1:nTot
        a = s0 + (m-1)*Ns + p.Ncp - back;
        F = fft(rx(a:a+N-1));
        Y(:, m) = F(bu+1);
    end
    % --- lệch tần từ 2 preamble giống nhau ---
    theta = angle(sum(Y(:,2) .* conj(Y(:,1))));
    res.cfoHz = theta * p.fs / (2*pi*Ns);
    Y = Y .* exp(-1j*theta*(0:nTot-1));
    % --- ước lượng kênh LS ---
    Hest = mean(Y(:,1:2), 2) ./ ref.pre;
    res.Hest = Hest;
    nv   = mean(abs(Y(:,1) - Y(:,2)).^2) / 2;
    sigH = mean(abs(Hest).^2) - nv/2;
    res.snrEst = 10*log10(max(sigH, eps) / nv);
    % --- cân bằng ZF + bám pha bằng pilot ---
    [~, ip] = ismember(p.pilotIdx, p.used);
    [~, id] = ismember(p.dataIdx,  p.used);
    Xd = Y(id, 3:end) ./ Hest(id);
    Xp = Y(ip, 3:end) ./ Hest(ip);
    ang  = angle(Xp .* conj(p.pilotVal));                % 4 x nSym
    A    = [ones(numel(p.pilotIdx),1), p.pilotIdx(:)];
    coef = A \ ang;                                      % pha = a + b*k (b ứng với lệch thời gian còn dư)
    Xd   = Xd .* exp(-1j*(coef(1,:) + p.dataIdx(:)*coef(2,:)));
    res.Xd = Xd(:);
    % --- giải điều chế ---
    rbits = qam_demod(Xd(:), p.M);
    res.nErr = sum(rbits ~= ref.bits); res.nBit = numel(ref.bits);
    res.ber  = res.nErr / res.nBit;
    ideal = qam_mod(ref.bits, p.M);
    res.evmdB = 20*log10(sqrt(mean(abs(Xd(:) - ideal).^2)));
end

function c = xcorr_fft(x, r)
% c(i) = sum_j x(i+j-1) * r(j): tương quan chéo tính bằng FFT.
    n = 2^nextpow2(numel(x) + numel(r));
    c = real(ifft(fft(x, n) .* conj(fft(r, n))));
    c = c(1:numel(x));
end

function sym = qam_mod(bits, M)
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