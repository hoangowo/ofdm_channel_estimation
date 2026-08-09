import numpy as np

class LSEstimator:
    """
    Module thực hiện các thuật toán xử lý tín hiệu trọng tâm:
    - Ước lượng kênh LS (Least Squares).
    - Ước lượng kênh Simplified MMSE (Wiener Filter).
    - Nội suy tuyến tính bảo toàn biên tuần hoàn (Linear Interpolation with Periodic Extension).
    - Cân bằng kênh miền tần số 1-tap Zero-Forcing (ZF Equalization).
    - Tính toán sai số bình phương trung bình (MSE).
    """
    def __init__(self):
        pass

    def estimate(self, Y_freq, pilot_indices, pilot_value, n_fft):
        """
        THUẬT TOÁN ƯỚC LƯỢNG KÊNH LS (Least Squares):
        - Công thức: H_hat_LS = Y_pilot / X_pilot
        - Sau khi ước lượng tại các subcarrier pilot, thực hiện nội suy tuyến tính (Linear Interpolation) 
          ra toàn bộ N_FFT (64) subcarriers.
        """
        # 1. Trích xuất tín hiệu thu tại các subcarrier pilot
        Y_pilots = Y_freq[pilot_indices]
        
        # Áp dụng công thức LS: Chia giá trị thu được cho giá trị pilot gốc đã biết
        H_ls_pilots = Y_pilots / pilot_value

        # 2. Tạo mảng chỉ số toàn bộ 64 subcarriers [0, 1, 2, ..., 63]
        all_subcarriers = np.arange(n_fft)
        
        # Mở rộng biên tuần hoàn (Periodic boundary extension) để xử lý méo ngoại suy ở 2 đầu subcarrier biên:
        # Nếu vị trí pilot cuối cùng (ví dụ: 60) nhỏ hơn subcarrier cuối cùng (63):
        # Ta giả định vị trí thứ n_fft (64) có giá trị đáp ứng kênh bằng đúng vị trí pilot đầu tiên (vị trí 0)
        # do tính chất tuần hoàn của phép biến đổi FFT.
        if pilot_indices[-1] < n_fft - 1:
            xp = np.append(pilot_indices, n_fft)
            fp_real = np.append(np.real(H_ls_pilots), np.real(H_ls_pilots[0]))
            fp_imag = np.append(np.imag(H_ls_pilots), np.imag(H_ls_pilots[0]))
        else:
            xp = pilot_indices
            fp_real = np.real(H_ls_pilots)
            fp_imag = np.imag(H_ls_pilots)

        # 3. Nội suy tuyến tính riêng rẽ phần Thực (Real) và phần Ảo (Imag) bằng np.interp
        H_est_real = np.interp(all_subcarriers, xp, fp_real)
        H_est_imag = np.interp(all_subcarriers, xp, fp_imag)
        
        # Hợp nhất lại thành đáp ứng kênh phức ước lượng H_est
        H_est = H_est_real + 1j * H_est_imag

        return H_est

    def estimate_mmse(self, Y_freq, pilot_indices, pilot_value, n_fft, snr_db):
        """
        THUẬT TOÁN ƯỚC LƯỢNG KÊNH SIMPLIFIED MMSE:
        - Công thức Wiener thu nhỏ: H_hat_MMSE = [SNR_linear / (SNR_linear + 1)] * H_hat_LS
        - Thuật toán này nhân thêm hệ số co (shrinkage factor) để giảm bớt năng lượng nhiễu AWGN 
          tác động lên các pilot tone, đặc biệt hiệu quả ở vùng SNR thấp.
        """
        # Chuyển đổi SNR từ dB sang tuyến tính
        snr_linear = 10.0 ** (snr_db / 10.0)
        
        # Tính hệ số lọc Wiener (luôn < 1)
        mmse_coeff = snr_linear / (snr_linear + 1.0)
        
        # 1. Thực hiện ước lượng LS tại các pilot trước
        Y_pilots = Y_freq[pilot_indices]
        H_ls_pilots = Y_pilots / pilot_value
        
        # 2. Nhân với hệ số co MMSE để giảm ảnh hưởng của nhiễu nền
        H_mmse_pilots = mmse_coeff * H_ls_pilots
        
        # 3. Tiến hành nội suy tuyến tính tuần hoàn biên giống như thuật toán LS
        all_subcarriers = np.arange(n_fft)
        if pilot_indices[-1] < n_fft - 1:
            xp = np.append(pilot_indices, n_fft)
            fp_real = np.append(np.real(H_mmse_pilots), np.real(H_mmse_pilots[0]))
            fp_imag = np.append(np.imag(H_mmse_pilots), np.imag(H_mmse_pilots[0]))
        else:
            xp = pilot_indices
            fp_real = np.real(H_mmse_pilots)
            fp_imag = np.imag(H_mmse_pilots)

        H_est_real = np.interp(all_subcarriers, xp, fp_real)
        H_est_imag = np.interp(all_subcarriers, xp, fp_imag)
        return H_est_real + 1j * H_est_imag

    def equalize(self, Y_data, H_est_data):
        """
        CÂN BẰNG KÊNH MIỀN TẦN SỐ 1-TAP ZERO-FORCING (ZF):
        - Công thức: X_hat = Y / H_est
        - Chức năng: Chia tín hiệu thu (Y_data) cho đáp ứng kênh ước lượng (H_est_data) 
          để đảo ngược tác động làm méo tín hiệu của kênh truyền.
        """
        # Thêm một hằng số cực nhỏ eps để tránh lỗi chia cho số 0
        eps = 1e-12
        H_safe = np.where(np.abs(H_est_data) == 0, eps, H_est_data)
        
        # Thực hiện chia 1-tap từng phần tử
        return Y_data / H_safe

    def equalize_mmse(self, Y_data, H_est_data, snr_db):
        """
        CÂN BẰNG KÊNH MIỀN TẦN SỐ 1-TAP MMSE (Minimum Mean Square Error):
        - Công thức: X_hat = H_est^* * Y / (|H_est|^2 + 1/SNR_linear)
        - So với ZF Equalization, MMSE Equalization không khuếch đại nhiễu tại các subcarrier 
          bị deep fade (|H(k)| ≈ 0), do đó cải thiện BER đáng kể ở vùng SNR thấp.
        - Khi SNR → ∞, MMSE Equalization tiệm cận ZF Equalization.
        """
        snr_linear = 10.0 ** (snr_db / 10.0)
        
        # Tính liên hợp phức (conjugate) của đáp ứng kênh ước lượng
        H_conj = np.conj(H_est_data)
        
        # Mẫu số MMSE: |H|^2 + 1/SNR (cân bằng giữa đảo ngược kênh và triệt nhiễu)
        denominator = np.abs(H_est_data)**2 + 1.0 / snr_linear
        
        # Bộ lọc MMSE 1-tap: X_hat = H* · Y / (|H|^2 + 1/SNR)
        return H_conj * Y_data / denominator

    def compute_mse(self, H_true, H_est):
        """
        ĐO SAI SỐ BÌNH PHƯƠNG TRUNG BÌNH (MSE):
        - Tính toán sự sai lệch bình phương trung bình giữa đáp ứng tần số kênh thực tế (H_true) 
          và đáp ứng kênh ước lượng (H_est).
        - Công thức: MSE = 1/N * sum( |H_true - H_est|^2 )
        """
        return np.mean(np.abs(H_true - H_est)**2)
