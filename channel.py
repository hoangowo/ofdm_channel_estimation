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
