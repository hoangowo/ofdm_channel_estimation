import numpy as np


class MultipathRayleighChannel:
    """
    Module mô phỏng Kênh truyền Rayleigh Fading đa đường (Multipath Rayleigh Fading Channel)
    kết hợp với Nhiễu Phức AWGN (Additive White Gaussian Noise).
    """
    def __init__(self, pdp):
        """
        Khởi tạo thông số kênh truyền với Hồ sơ Phân bố Năng lượng Trễ - Power Delay Profile (PDP).
        - pdp: Mảng chứa công suất trung bình của các đường truyền (taps).
               Ví dụ: [0.7, 0.2, 0.1] cho kênh 3 đường (3-tap).
        """
        # Chuẩn hóa mảng PDP sao cho tổng năng lượng tất cả các đường truyền đúng bằng 1.0 (0.7 + 0.2 + 0.1 = 1.0)
        self.pdp = np.array(pdp, dtype=float)
        self.pdp = self.pdp / np.sum(self.pdp)
        
        # Số lượng đường truyền (taps) của kênh truyền
        self.num_taps = len(self.pdp)

    def generate_taps(self):
        """
        Sinh các hệ số tap kênh Rayleigh Fading phức ngẫu nhiên: h_i ~ CN(0, P_i).
        Mỗi tap là một số phức Gauss có:
        - Phần Thực (Real) ~ N(0, P_i / 2)
        - Phần Ảo (Imag) ~ N(0, P_i / 2)
        - Biên độ |h_i| tuân theo Phân bố Rayleigh.
        - Công suất trung bình E[|h_i|^2] = P_i.
        """
        # Sinh phần Thực và Ảo độc lập theo phân bố chuẩn Gauss (chuẩn hóa độ lệch chuẩn theo sqrt(P_i / 2))
        real_part = np.random.randn(self.num_taps) * np.sqrt(self.pdp / 2.0)
        imag_part = np.random.randn(self.num_taps) * np.sqrt(self.pdp / 2.0)
        
        # Ghép thành số phức h = Real + j*Imag
        h_taps = real_part + 1j * imag_part
        return h_taps

    def get_frequency_response(self, h_taps, n_fft):
        """
        Tính Đáp ứng Tần số Thực tế (True Channel Frequency Response) H_true của kênh truyền.
        - Chuyển đổi Đáp ứng Xung thời gian h_taps (độ dài L) sang miền Tần số H_true (độ dài N_FFT = 64) 
          bằng phép biến đổi FFT N_FFT điểm: H_true = FFT(h_taps, n_fft).
        - H_true đóng vai trò là "Perfect CSI" (Thông tin kênh truyền hoàn hảo) để làm cơ sở so sánh.
        """
        return np.fft.fft(h_taps, n_fft)

    def apply_channel(self, x_cp, h_taps, snr_db):
        """
        Truyền tín hiệu OFDM miền thời gian x_cp qua kênh Fading đa đường và cộng Nhiễu AWGN.
        Quy trình:
          1. Tích chập tuyến tính tín hiệu x_cp với h_taps: y_conv = x_cp * h.
          2. Cắt độ dài y_conv về đúng độ dài tín hiệu ban đầu.
          3. Tính công suất tín hiệu và phương sai nhiễu sigma_z^2 dựa trên SNR_dB.
          4. Sinh nhiễu phức Gauss AWGN z ~ CN(0, sigma_z^2) và cộng vào tín hiệu.
        """
        # 1. Tích chập tuyến tính miền thời gian giữa tín hiệu phát x_cp và các tap kênh h_taps
        y_conv = np.convolve(x_cp, h_taps, mode='full')
        
        # Cắt lấy đoạn đầu có độ dài bằng đúng x_cp (do CP đã hấp thụ độ trễ đa đường)
        y_cp_noiseless = y_conv[:len(x_cp)]
        
        # 2. Cộng nhiễu AWGN
        # Tính công suất tín hiệu thu chưa nhiễu (Signal Power)
        signal_power = np.mean(np.abs(y_cp_noiseless)**2)
        if signal_power == 0:
            signal_power = 1.0
            
        # Chuyển SNR từ đơn vị dB sang tỷ lệ tuyến tính (Linear SNR)
        # SNR_linear = 10^(SNR_dB / 10)
        snr_linear = 10.0 ** (snr_db / 10.0)
        
        # Tính phương sai công suất nhiễu: noise_variance = Signal_Power / SNR_linear
        noise_variance = signal_power / snr_linear
        
        # Sinh nhiễu phức AWGN z ~ CN(0, noise_variance)
        # Phần thực và ảo chia đôi phương sai: variance_part = noise_variance / 2.0
        noise_real = np.random.randn(len(x_cp)) * np.sqrt(noise_variance / 2.0)
        noise_imag = np.random.randn(len(x_cp)) * np.sqrt(noise_variance / 2.0)
        noise = noise_real + 1j * noise_imag
        
        # Tín hiệu thu hoàn chỉnh bao gồm cả Fading và Nhiễu AWGN
        y_cp = y_cp_noiseless + noise
        return y_cp


class TimeVaryingRayleighChannel:
    """
    Module mô phỏng Kênh Rayleigh Fading Đa đường Biến thiên theo Thời gian
    (Time-Varying Multipath Rayleigh Fading Channel) sử dụng Mô hình Jakes.

    Mô hình Jakes (Sum-of-Sinusoids - SoS):
    - Sinh kênh Rayleigh phức có tương quan thời gian tuân theo hàm Bessel bậc 0:
        R_H(Δn) = J_0(2π * fd_norm * Δn)
    - Trong đó fd_norm = f_d * T_sym là tần số Doppler chuẩn hóa.
    - Mỗi tap kênh h_l[n] được tổng hợp từ num_sinusoids sóng sin độc lập.

    Hỗ trợ 2 chế độ truyền Frame:
    - 'block_aging': Kênh không đổi trong 1 OFDM symbol, biến đổi giữa các symbol (phổ biến).
    - 'sample_varying': Kênh biến thiên từng mẫu thời gian (gây ICI, dùng để khảo sát).
    """

    def __init__(self, pdp, fd_norm, n_fft, n_cp, num_symbols_per_frame,
                 num_sinusoids=8):
        """
        Khởi tạo kênh biến thiên theo thời gian.

        Triển khai theo đúng Program 2.20 (trang 53-54) của giáo trình:
        "MIMO-OFDM Wireless Communications with MATLAB" - Cho, Kim, Yang, Kang (2010).

        Tham số:
        - pdp              : Hồ sơ Phân bố Năng lượng Trễ (list/array), ví dụ [0.7, 0.2, 0.1].
        - fd_norm          : Tần số Doppler chuẩn hóa fd_norm = f_d * T_sym.
                             fd_norm = 0.0 -> kênh tĩnh; fd_norm = 0.1 -> cao tốc.
        - n_fft            : Kích thước FFT (số sóng mang phụ).
        - n_cp             : Độ dài Cyclic Prefix.
        - num_symbols_per_frame : Số ký tự OFDM trong một frame.
        - num_sinusoids    : Số sóng sin N0 trong mô hình Jakes.
                             Mặc định N0=8 ("As suggested by Jakes", Program 2.20 trang 53).
                             Số xấp xỉ thực tế N = 4*N0 + 2 = 34 sóng mang đa đường.
        """
        # Chuẩn hóa PDP để tổng công suất = 1
        self.pdp = np.array(pdp, dtype=float)
        self.pdp = self.pdp / np.sum(self.pdp)

        self.num_taps              = len(self.pdp)      # Số đường truyền (taps)
        self.fd_norm               = fd_norm            # Tần số Doppler chuẩn hóa
        self.n_fft                 = n_fft              # Kích thước FFT
        self.n_cp                  = n_cp               # Độ dài CP
        self.num_symbols_per_frame = num_symbols_per_frame  # Số symbol/frame
        self.num_sinusoids         = num_sinusoids      # N0 = số sóng sin mô hình Jakes

        # Độ dài tín hiệu thời gian của 1 OFDM symbol (gồm cả CP)
        self.symbol_len = n_fft + n_cp

        # =====================================================================
        # Tính các thông số Jakes theo Program 2.20 (trang 53) và Eq. 2.26 (trang 51)
        # =====================================================================
        # N = 4*N0 + 2: Tổng số sóng mang đa đường xấp xỉ (Eq. 2.19, trang 50)
        N_actual = 4 * num_sinusoids + 2

        # Góc tới của sóng phẳng thứ n: theta_n = 2*pi*n / N (Eq. 2.26, trang 51)
        # n = 1, 2, ..., N0
        n_idx = np.arange(1, num_sinusoids + 1)
        self._theta = 2.0 * np.pi * n_idx / N_actual  # (N0,)

        # Pha ban đầu: phi_n = pi*n / (N0+1) (Eq. 2.24, trang 51) - phần tất định
        # Thêm offset ngẫu nhiên riêng cho từng tap để đảm bảo các tap kênh độc lập nhau
        phi_base = np.pi * n_idx / (num_sinusoids + 1)  # (N0,) - tất định
        random_offset = np.random.uniform(0, 2 * np.pi, (self.num_taps, num_sinusoids))
        self._phi = phi_base[None, :] + random_offset   # (num_taps, N0)

    def _generate_tap_at_sample(self, tap_idx, sample_index):
        """
        Sinh hệ số tap thứ tap_idx tại mẫu thời gian rời rạc sample_index.

        Công thức Jakes gốc theo Eq. 2.23a, 2.23b, 2.25 (trang 51) của giáo trình:

            h_I(t) = 2*sum_{n=1}^{N0} cos(phi_n)*cos(w_n*t) + sqrt(2)*cos(w_d*t)
            h_Q(t) = 2*sum_{n=1}^{N0} sin(phi_n)*cos(w_n*t) + 0  [phi_N=0 -> sin(0)=0]
            h(t)   = E0 / sqrt(2*N0 + 1) * (h_I(t) + j*h_Q(t))

        Trong đó:
        - w_n = 2*pi*fd_norm*cos(theta_n): Tần số góc dịch tần Doppler của sóng phẳng thứ n
        - w_d = 2*pi*fd_norm              : Tần số Doppler cực đại (specular component)
        - N0                              : Số sóng sin (num_sinusoids = 8)
        - E0 = sqrt(P_l)                  : Hệ số tỉ lệ để chuẩn hóa công suất tap
        - phi_N = 0                       : Pha ban đầu của thành phần specular (tất định)

        Hàm tự tương quan thời gian của kết quả thỏa mãn Eq. 2.27:
            E[h*(t) h(t+tau)] = P_l * J_0(2*pi*fd_norm*tau)
        """
        p_l   = self.pdp[tap_idx]  # Công suất của tap l từ PDP
        theta = self._theta         # (N0,) - góc tới theo Eq. 2.26
        phi   = self._phi[tap_idx]  # (N0,) - pha ban đầu theo Eq. 2.24 + offset ngẫu nhiên

        # Tần số góc của từng sóng phẳng: w_n = 2*pi*fd_norm*cos(theta_n) (Eq. 2.26)
        w_n = 2.0 * np.pi * self.fd_norm * np.cos(theta)  # (N0,)
        # Tần số Doppler cực đại (cho thành phần specular): w_d = 2*pi*fd_norm
        w_d = 2.0 * np.pi * self.fd_norm

        t = float(sample_index)  # Thời gian rời rạc (mẫu)

        # Thành phần In-phase (Eq. 2.23a), với phi_N=0 -> cos(phi_N)=1
        h_I = 2.0 * np.sum(np.cos(phi) * np.cos(w_n * t)) + np.sqrt(2.0) * np.cos(w_d * t)

        # Thành phần Quadrature (Eq. 2.23b), với phi_N=0 -> sin(phi_N)=0
        h_Q = 2.0 * np.sum(np.sin(phi) * np.cos(w_n * t))

        # Chuẩn hóa năng lượng tap (Eq. 2.25): E0/sqrt(2*N0+1) đảm bảo E[|h|^2] = P_l
        E0    = np.sqrt(p_l)
        h_tap = (E0 / np.sqrt(2 * self.num_sinusoids + 1)) * (h_I + 1j * h_Q)
        return h_tap

    def _generate_h_at_symbol(self, symbol_idx):
        """
        Sinh vectơ tap kênh h[l] tại thời điểm đầu của OFDM symbol thứ symbol_idx.
        Chế độ 'block_aging': kênh được lấy mẫu 1 lần tại thời điểm n = symbol_idx * symbol_len.

        Trả về:
        - h_taps : Mảng phức kích thước (num_taps,)
        """
        # Mẫu thời gian tương ứng với đầu symbol (chỉ số mẫu tuyệt đối)
        n_sample = symbol_idx * self.symbol_len

        h_taps = np.array([
            self._generate_tap_at_sample(l, n_sample)
            for l in range(self.num_taps)
        ])
        return h_taps

    def _apply_channel_one_symbol_block(self, x_cp, h_taps, noise_std):
        """
        Truyền 1 OFDM symbol qua kênh block-fading (kênh tĩnh trong symbol).

        Tham số:
        - x_cp      : Tín hiệu thời gian miền phát, dài symbol_len = N_FFT + N_CP.
        - h_taps    : Vectơ tap kênh phức tại symbol này, dài num_taps.
        - noise_std : Độ lệch chuẩn AWGN (mỗi phần thực/ảo).

        Trả về:
        - y_cp : Tín hiệu thu (đã qua fading + AWGN), dài symbol_len.
        """
        # Tích chập tuyến tính và cắt về đúng độ dài (CP đã hấp thụ ISI)
        y_conv = np.convolve(x_cp, h_taps, mode='full')
        y_noiseless = y_conv[:len(x_cp)]

        # Sinh nhiễu AWGN phức: phần thực và ảo mỗi cái ~ N(0, noise_std^2)
        n_len = len(x_cp)
        noise = (np.random.randn(n_len) + 1j * np.random.randn(n_len)) * noise_std
        return y_noiseless + noise

    def _apply_channel_one_symbol_sample_varying(self, symbol_idx, x_cp, noise_std):
        """
        Truyền 1 OFDM symbol qua kênh biến thiên từng mẫu (sample-varying).
        Kênh thay đổi tại MỖI mẫu thời gian n trong symbol, gây ra nhiễu ICI.

        Tham số:
        - symbol_idx : Chỉ số symbol trong frame (để tính n tuyệt đối).
        - x_cp       : Tín hiệu thời gian miền phát, dài symbol_len.
        - noise_std  : Độ lệch chuẩn AWGN.

        Trả về:
        - y_cp : Tín hiệu thu (mỗi mẫu qua kênh khác nhau + AWGN).
        """
        n_samples = len(x_cp)
        y_noiseless = np.zeros(n_samples, dtype=complex)

        for n in range(n_samples):
            # Mẫu tuyệt đối trong toàn bộ chuỗi thời gian
            n_abs = symbol_idx * self.symbol_len + n

            # Sinh tap kênh tại mẫu n_abs
            h_taps_n = np.array([
                self._generate_tap_at_sample(l, n_abs)
                for l in range(self.num_taps)
            ])

            # Đóng góp của mẫu n vào đầu ra: nhân tích chập với các tap trễ
            for l, h_l in enumerate(h_taps_n):
                if n - l >= 0:
                    y_noiseless[n] += h_l * x_cp[n - l]

        # Cộng nhiễu AWGN
        noise = (np.random.randn(n_samples) + 1j * np.random.randn(n_samples)) * noise_std
        return y_noiseless + noise

    def transmit_frame(self, X_frame_cp, snr_db, mode='block_aging'):
        """
        Truyền toàn bộ một Frame OFDM qua kênh biến thiên thời gian và cộng AWGN.

        Tham số:
        - X_frame_cp : Tín hiệu thời gian của cả frame, kích thước (num_symbols_per_frame, symbol_len).
                       Mỗi hàng là 1 OFDM symbol (đã có CP) sẵn sàng phát.
        - snr_db     : Tỷ lệ tín hiệu trên nhiễu (dB).
        - mode       : Chế độ kênh:
                       'block_aging'    - Kênh tĩnh trong symbol, biến đổi giữa các symbol (mặc định).
                       'sample_varying' - Kênh biến thiên từng mẫu trong symbol (khảo sát ICI).

        Trả về:
        - Y_frame_cp    : Tín hiệu thu của cả frame, kích thước (num_symbols_per_frame, symbol_len).
        - H_true_frame  : Ma trận đáp ứng tần số thực sự tại mỗi symbol,
                          kích thước (N_FFT, num_symbols_per_frame).
                          Cột thứ n là H_true[:, n] = FFT(h_taps[n], N_FFT).
        """
        # Ước tính công suất tín hiệu trung bình (từ toàn bộ frame)
        sig_power = np.mean(np.abs(X_frame_cp) ** 2)
        if sig_power == 0:
            sig_power = 1.0

        # Tính phương sai nhiễu dựa trên SNR
        snr_linear    = 10.0 ** (snr_db / 10.0)
        noise_var     = sig_power / snr_linear
        # Độ lệch chuẩn AWGN (mỗi phần thực/ảo): noise_std = sqrt(noise_var / 2)
        noise_std     = np.sqrt(noise_var / 2.0)

        # Khởi tạo các mảng đầu ra
        Y_frame_cp   = np.zeros_like(X_frame_cp)
        H_true_frame = np.zeros((self.n_fft, self.num_symbols_per_frame), dtype=complex)

        for sym_idx in range(self.num_symbols_per_frame):
            # Sinh tap kênh tại symbol này
            h_taps_sym = self._generate_h_at_symbol(sym_idx)

            # Đáp ứng tần số thực tế H_true[:, sym_idx] = FFT(h_taps, N_FFT)
            H_true_frame[:, sym_idx] = np.fft.fft(h_taps_sym, self.n_fft)

            # Truyền tín hiệu qua kênh theo chế độ đã chọn
            x_cp_sym = X_frame_cp[sym_idx]

            if mode == 'block_aging':
                Y_frame_cp[sym_idx] = self._apply_channel_one_symbol_block(
                    x_cp_sym, h_taps_sym, noise_std
                )
            elif mode == 'sample_varying':
                Y_frame_cp[sym_idx] = self._apply_channel_one_symbol_sample_varying(
                    sym_idx, x_cp_sym, noise_std
                )
            else:
                raise ValueError(
                    f"Chế độ kênh '{mode}' không hợp lệ! Dùng 'block_aging' hoặc 'sample_varying'."
                )

        return Y_frame_cp, H_true_frame
