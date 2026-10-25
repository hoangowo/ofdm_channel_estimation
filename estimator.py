import numpy as np


class LSEstimator:
    """
    Module thực hiện các thuật toán xử lý tín hiệu trọng tâm:
    - Ước lượng kênh LS (Least Squares) 1D (baseline - giữ nguyên).
    - Ước lượng kênh Simplified MMSE (Wiener Filter) 1D.
    - Nội suy tuyến tính bảo toàn biên tuần hoàn (Linear Interpolation with Periodic Extension).
    - Cân bằng kênh miền tần số 1-tap Zero-Forcing (ZF Equalization).
    - Tính toán sai số bình phương trung bình (MSE).
    - Ước lượng kênh LS 2D cho frame (MỚI).
    - Ước lượng kênh DFT-Based 2D cho frame (MỚI).
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

    # ===========================================================================
    # PHẦN MỚI: ƯỚC LƯỢNG KÊNH 2D CHO FRAME (FREQUENCY x TIME)
    # ===========================================================================

    def _interp_freq_axis(self, H_pilots_at_subcarriers, pilot_subcarriers, n_fft):
        """
        Hàm phụ: Nội suy tuyến tính theo chiều tần số (frequency axis) với mở rộng biên tuần hoàn.
        Được dùng chung cho cả LS 2D và DFT-based 2D.

        Tham số:
        - H_pilots_at_subcarriers : Mảng phức 1D kích thước (num_pilot_freq,) - giá trị kênh tại pilot.
        - pilot_subcarriers        : Mảng 1D chỉ số subcarrier pilot.
        - n_fft                    : Số sóng mang phụ (FFT size).

        Trả về:
        - H_interp : Mảng phức 1D kích thước (n_fft,) - đã nội suy đầy đủ.
        """
        all_k = np.arange(n_fft)

        # Mở rộng biên tuần hoàn: nếu pilot cuối chưa tới subcarrier cuối, giả định tuần hoàn
        if pilot_subcarriers[-1] < n_fft - 1:
            xp      = np.append(pilot_subcarriers, n_fft)
            fp_real = np.append(np.real(H_pilots_at_subcarriers), np.real(H_pilots_at_subcarriers[0]))
            fp_imag = np.append(np.imag(H_pilots_at_subcarriers), np.imag(H_pilots_at_subcarriers[0]))
        else:
            xp      = pilot_subcarriers
            fp_real = np.real(H_pilots_at_subcarriers)
            fp_imag = np.imag(H_pilots_at_subcarriers)

        H_real = np.interp(all_k, xp, fp_real)
        H_imag = np.interp(all_k, xp, fp_imag)
        return H_real + 1j * H_imag

    def _interp_time_axis(self, H_at_pilot_syms, pilot_sym_indices, num_symbols):
        """
        Hàm phụ: Nội suy tuyến tính theo chiều thời gian (time axis) cho toàn bộ frame.
        Dùng nội suy tuyến tính giữa các symbol pilot, giữ nguyên giá trị ở biên (hold).

        Tham số:
        - H_at_pilot_syms   : Ma trận phức kích thước (n_fft, num_pilot_syms) - H tại các cột pilot.
        - pilot_sym_indices : Mảng 1D chỉ số symbol (cột) là pilot.
        - num_symbols       : Tổng số symbol trong frame.

        Trả về:
        - H_full : Ma trận phức kích thước (n_fft, num_symbols) - đã nội suy đầy đủ chiều thời gian.
        """
        n_fft        = H_at_pilot_syms.shape[0]
        all_n        = np.arange(num_symbols)
        H_full_real  = np.zeros((n_fft, num_symbols))
        H_full_imag  = np.zeros((n_fft, num_symbols))

        for k in range(n_fft):
            # Nội suy từng subcarrier theo chiều thời gian
            # np.interp giữ giá trị đầu/cuối ở biên (không ngoại suy) nhờ tham số left/right mặc định
            H_full_real[k] = np.interp(all_n, pilot_sym_indices,
                                        np.real(H_at_pilot_syms[k]))
            H_full_imag[k] = np.interp(all_n, pilot_sym_indices,
                                        np.imag(H_at_pilot_syms[k]))

        return H_full_real + 1j * H_full_imag

    def estimate_ls_2d(self, Y_frame, pilot_subcarriers, pilot_symbols, pilot_value, n_fft, num_symbols):
        """
        ƯỚC LƯỢNG KÊNH LS 2D CHO TOÀN BỘ FRAME OFDM.

        Thuật toán thực hiện theo 3 bước:
          Bước 1: Ước lượng LS tại các nút pilot (k_p, n_p).
                  H_hat_LS[k_p, n_p] = Y[k_p, n_p] / pilot_value

          Bước 2: Nội suy tuyến tính theo chiều TẦN SỐ với mở rộng biên tuần hoàn,
                  tại mỗi cột pilot_symbol n_p.
                  Kết quả: H_ls_at_pilot_syms[k, n_p] với k = 0..N_FFT-1.

          Bước 3: Nội suy tuyến tính theo chiều THỜI GIAN giữa các cột pilot_symbol,
                  với giá trị biên được giữ cố định (hold) để tránh ngoại suy.
                  Kết quả: H_ls_frame[k, n] với k = 0..N_FFT-1, n = 0..num_symbols-1.

        Tham số:
        - Y_frame          : Ma trận tín hiệu thu 2D, kích thước (N_FFT, num_symbols).
        - pilot_subcarriers: Mảng 1D chỉ số subcarrier pilot (k_p).
        - pilot_symbols    : Mảng 1D chỉ số symbol pilot (n_p).
        - pilot_value      : Giá trị pilot đã biết trước.
        - n_fft            : Kích thước FFT.
        - num_symbols      : Số OFDM symbol trong frame.

        Trả về:
        - H_ls_frame : Ma trận ước lượng kênh 2D, kích thước (N_FFT, num_symbols).
        """
        num_pilot_syms  = len(pilot_symbols)
        # Ma trận lưu giá trị kênh đã nội suy tần số, tại các cột pilot symbol
        H_at_pilot_syms = np.zeros((n_fft, num_pilot_syms), dtype=complex)

        for idx, n_p in enumerate(pilot_symbols):
            # Bước 1: LS tại các pilot subcarrier của symbol n_p
            Y_pilots_at_np = Y_frame[pilot_subcarriers, n_p]
            H_ls_at_pilots  = Y_pilots_at_np / pilot_value

            # Bước 2: Nội suy theo chiều tần số (periodic extension)
            H_at_pilot_syms[:, idx] = self._interp_freq_axis(
                H_ls_at_pilots, pilot_subcarriers, n_fft
            )

        # Bước 3: Nội suy theo chiều thời gian giữa các cột pilot symbol
        H_ls_frame = self._interp_time_axis(H_at_pilot_syms, pilot_symbols, num_symbols)
        return H_ls_frame

    def estimate_dft_2d(self, Y_frame, pilot_subcarriers, pilot_symbols, pilot_value,
                        n_fft, num_symbols, window_len=None):
        """
        ƯỚC LƯỢNG KÊNH DFT-BASED 2D CHO TOÀN BỘ FRAME OFDM.

        Thuật toán kết hợp nội suy LS 2D với bộ lọc DFT miền trễ (Delay-Domain Filter)
        để triệt nhiễu AWGN ngoài vùng trễ của kênh:

          Bước 1: Thực hiện LS 2D như trên để thu H_ls_frame (N_FFT x num_symbols).

          Bước 2: Với mỗi OFDM symbol n (sau khi đã có H_ls toàn frame):
                  a. Biến đổi IFFT vào miền trễ: h_est[l] = IFFT(H_ls[:, n])
                  b. Cửa sổ hóa (Windowing) trong miền trễ:
                     - Giữ lại các mẫu l = 0..window_len-1 (năng lượng kênh thực tế).
                     - Đặt các mẫu l = window_len..N_FFT-1 = 0 (triệt nhiễu AWGN ngoài vùng trễ).
                  c. Biến đổi FFT ngược để thu đáp ứng kênh tần số đã lọc: H_dft[:, n] = FFT(h_windowed)

        Ưu điểm so với LS 2D:
        - Triệt nhiễu AWGN hiệu quả bằng cách chỉ giữ năng lượng trong vùng trễ kênh.
        - MSE giảm rõ rệt ở vùng SNR thấp và trung bình.
        - Hiệu quả nhất khi window_len = len(CHANNEL_PDP) (đúng với độ trễ kênh thực tế).

        Tham số:
        - Y_frame          : Ma trận tín hiệu thu 2D, kích thước (N_FFT, num_symbols).
        - pilot_subcarriers: Mảng 1D chỉ số subcarrier pilot.
        - pilot_symbols    : Mảng 1D chỉ số symbol pilot.
        - pilot_value      : Giá trị pilot đã biết.
        - n_fft            : Kích thước FFT.
        - num_symbols      : Số OFDM symbol trong frame.
        - window_len       : Độ dài cửa sổ trễ (số mẫu được giữ lại).
                             Nếu None, mặc định dùng n_cp (Cyclic Prefix length) / 2 hoặc
                             len(channel_pdp). Có thể truyền N_CP hoặc len(CHANNEL_PDP).
                             Thông thường: window_len = 16 (N_CP) hoặc window_len = 3 (len(PDP)).

        Trả về:
        - H_dft_frame : Ma trận ước lượng kênh DFT 2D, kích thước (N_FFT, num_symbols).
        """
        # Bước 1: Khởi tạo ước lượng LS 2D đầy đủ frame
        H_ls_frame = self.estimate_ls_2d(
            Y_frame, pilot_subcarriers, pilot_symbols, pilot_value, n_fft, num_symbols
        )

        # Xác định window_len mặc định nếu không truyền vào
        if window_len is None:
            # Mặc định lấy N_FFT // 4 (bảo toàn an toàn) nhưng thực tế nên truyền vào
            window_len = max(1, n_fft // 4)

        # Bước 2: Áp dụng bộ lọc DFT miền trễ cho từng OFDM symbol
        H_dft_frame = np.zeros((n_fft, num_symbols), dtype=complex)

        for n in range(num_symbols):
            # a. Biến đổi IFFT: từ miền tần số sang miền trễ (delay domain)
            #    h_est[l] = IFFT(H_ls[:, n]) - đây là đáp ứng xung ước lượng
            h_est = np.fft.ifft(H_ls_frame[:, n], axis=0)

            # b. Cửa sổ hóa trong miền trễ:
            #    - Tạo mặt nạ cửa sổ: 1 cho l = 0..window_len-1, 0 cho phần còn lại
            #    - Các mẫu bị đặt 0 là nhiễu AWGN "tràn vào" ngoài vùng trễ kênh thực tế
            h_windowed = np.zeros(n_fft, dtype=complex)
            h_windowed[:window_len] = h_est[:window_len]  # Giữ phần kênh

            # c. Biến đổi FFT ngược: quay lại miền tần số đã được lọc nhiễu
            #    H_dft[:, n] = FFT(h_windowed) - ước lượng kênh sau khi triệt nhiễu DFT
            H_dft_frame[:, n] = np.fft.fft(h_windowed, axis=0)

        return H_dft_frame
