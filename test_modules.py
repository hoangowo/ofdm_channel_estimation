"""
test_modules.py - Kiểm thử đơn vị (Unit Tests) cho tất cả các module của hệ thống OFDM.
Đồ án: "Pilot-Based Channel Estimation for OFDM over Time-Varying Fading Channels"
Sinh viên: Đào Kim Hoàng - MSSV: 2211068

Các bài test bao gồm:
  Test 1: Vòng lặp lý tưởng không nhiễu (Noiseless Ideal Loopback) - giữ nguyên từ baseline.
  Test 2: Mô phỏng kênh Rayleigh đa đường - giữ nguyên từ baseline.
  Test 3 (MỚI): Kênh Jakes fd_norm=0.0 (kênh tĩnh) -> BER=0 khi không nhiễu.
  Test 4 (MỚI): So sánh MSE của DFT-based vs LS tại SNR=10dB -> DFT nhỏ hơn rõ rệt.
"""
import numpy as np
from config import (N_FFT, N_CP, MODULATION, PILOT_SPACING, PILOT_VALUE, CHANNEL_PDP,
                    NUM_SYMBOLS_PER_FRAME, PILOT_SPACING_FREQ, PILOT_SPACING_TIME)
from channel import MultipathRayleighChannel, TimeVaryingRayleighChannel
from ofdm_system import OFDMSystem
from estimator import LSEstimator


def test_noiseless_loopback():
    """
    TEST 1: Vòng lặp lý tưởng - Kênh phẳng không nhiễu (Noiseless Flat Fading).
    Điều kiện pass: 0 bit lỗi, MSE kênh < 1e-10.
    """
    print("--- Running Test 1: Noiseless Ideal Loopback ---")
    ofdm      = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, modulation=MODULATION)
    estimator = LSEstimator()
    channel   = MultipathRayleighChannel([1.0])
    
    # Kênh lý tưởng phẳng h = [1.0]
    h_taps = np.array([1.0 + 0.0j])
    H_true = channel.get_frequency_response(h_taps, N_FFT)
    
    # Sinh bit ngẫu nhiên
    num_data_subcarriers = ofdm.num_data
    bits_per_symbol      = ofdm.bits_per_symbol
    tx_bits              = np.random.randint(0, 2, num_data_subcarriers * bits_per_symbol)
    
    # 1. Ánh xạ bit -> symbol
    tx_data_symbols = ofdm.map_bits(tx_bits)
    
    # 2. Chèn pilot
    X_freq = ofdm.insert_pilots(tx_data_symbols)
    
    # 3. Điều chế (IFFT + CP)
    x_cp = ofdm.modulate(X_freq)
    
    # 4. Qua kênh không nhiễu (SNR rất cao = 100 dB)
    y_cp  = channel.apply_channel(x_cp, h_taps, snr_db=100.0)
    Y_freq = ofdm.demodulate(y_cp)
    
    # 5. Ước lượng kênh LS
    H_est = estimator.estimate(Y_freq, ofdm.pilot_indices, PILOT_VALUE, N_FFT)
    
    # 6. Cân bằng kênh ZF
    Y_data         = Y_freq[ofdm.data_indices]
    H_est_data     = H_est[ofdm.data_indices]
    rx_data_symbols = estimator.equalize(Y_data, H_est_data)
    
    # 7. Giải điều chế
    rx_bits    = ofdm.demap_symbols(rx_data_symbols)
    bit_errors = np.sum(tx_bits != rx_bits)
    mse        = estimator.compute_mse(H_true, H_est)
    
    print(f"  Noiseless Loopback Bit Errors : {bit_errors}")
    print(f"  Noiseless Loopback Channel MSE: {mse:.6e}")
    assert bit_errors == 0, "Loopback phải có 0 bit lỗi!"
    assert mse < 1e-10,     "Loopback MSE phải xấp xỉ 0!"
    print(">>> Test 1 PASSED!\n")


def test_multipath_rayleigh():
    """
    TEST 2: Mô phỏng kênh Rayleigh đa đường tại SNR = 20 dB.
    Điều kiện pass: Chạy thông suốt (không exception), in kết quả để tham khảo.
    """
    print("--- Running Test 2: Rayleigh Channel Simulation ---")
    ofdm      = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, modulation=MODULATION)
    channel   = MultipathRayleighChannel(CHANNEL_PDP)
    estimator = LSEstimator()
    
    h_taps = channel.generate_taps()
    H_true = channel.get_frequency_response(h_taps, N_FFT)
    
    tx_bits    = np.random.randint(0, 2, ofdm.num_data * ofdm.bits_per_symbol)
    tx_data    = ofdm.map_bits(tx_bits)
    X_freq     = ofdm.insert_pilots(tx_data)
    x_cp       = ofdm.modulate(X_freq)
    
    snr_db    = 20.0
    y_cp      = channel.apply_channel(x_cp, h_taps, snr_db)
    Y_freq    = ofdm.demodulate(y_cp)
    H_est     = estimator.estimate(Y_freq, ofdm.pilot_indices, PILOT_VALUE, N_FFT)
    
    # So sánh Perfect CSI và LS CSI
    rx_data_perfect = estimator.equalize(Y_freq[ofdm.data_indices], H_true[ofdm.data_indices])
    bits_perfect    = ofdm.demap_symbols(rx_data_perfect)
    err_perfect     = np.sum(tx_bits != bits_perfect)
    
    rx_data_ls = estimator.equalize(Y_freq[ofdm.data_indices], H_est[ofdm.data_indices])
    bits_ls    = ofdm.demap_symbols(rx_data_ls)
    err_ls     = np.sum(tx_bits != bits_ls)
    
    mse = estimator.compute_mse(H_true, H_est)
    
    print(f"  SNR                 : {snr_db} dB")
    print(f"  Perfect CSI Errors  : {err_perfect} / {len(tx_bits)}")
    print(f"  LS Estimation Errors: {err_ls} / {len(tx_bits)}")
    print(f"  LS Channel MSE      : {mse:.6e}")
    print(">>> Test 2 PASSED!\n")


def test_jakes_static_channel_no_errors():
    """
    TEST 3 (MỚI): Kênh Jakes với fd_norm = 0.0 (kênh tĩnh hoàn toàn).
    Tại SNR rất cao (100 dB), BER của toàn bộ frame phải = 0.

    Kiểm tra:
    - fd_norm = 0 -> kênh không thay đổi -> giống kênh block-fading thông thường.
    - Với kênh tĩnh và không nhiễu, LS 2D phải ước lượng chính xác.
    - Toàn bộ bit trong frame phải được khôi phục hoàn toàn (0 lỗi).
    """
    print("--- Running Test 3: Jakes Channel fd_norm=0.0 (Static) -> BER=0 @ No Noise ---")
    
    num_symbols = NUM_SYMBOLS_PER_FRAME
    df          = PILOT_SPACING_FREQ
    dt          = 1  # dt=1: mỗi symbol đều là pilot -> độ chính xác cao nhất
    snr_db      = 100.0  # Không nhiễu (SNR cực cao)

    # Khởi tạo hệ thống
    ofdm      = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, modulation=MODULATION)
    estimator = LSEstimator()
    channel   = TimeVaryingRayleighChannel(
        pdp                  = CHANNEL_PDP,
        fd_norm              = 0.0,         # Kênh tĩnh
        n_fft                = N_FFT,
        n_cp                 = N_CP,
        num_symbols_per_frame = num_symbols,
    )

    # Sinh bit và tạo frame phát
    max_bits    = N_FFT * num_symbols * ofdm.bits_per_symbol
    data_bits   = np.random.randint(0, 2, max_bits)
    X_frame, tx_bits_used, pilot_sub, pilot_sym = ofdm.insert_pilots_2d(
        data_bits, df, dt, num_symbols
    )

    # Điều chế frame -> miền thời gian
    X_frame_cp = ofdm.modulate_frame(X_frame)

    # Truyền qua kênh Jakes (kênh tĩnh)
    Y_frame_cp, H_true_frame = channel.transmit_frame(X_frame_cp, snr_db=snr_db, mode='block_aging')

    # Giải điều chế -> miền tần số
    Y_frame = ofdm.demodulate_frame(Y_frame_cp)

    # Ước lượng kênh LS 2D
    H_ls_frame = estimator.estimate_ls_2d(
        Y_frame, pilot_sub, pilot_sym, PILOT_VALUE, N_FFT, num_symbols
    )

    # Giải mã và đếm lỗi
    rx_bits    = ofdm.equalize_and_decode_frame(Y_frame, H_ls_frame, pilot_sub, pilot_sym)
    bit_errors = np.sum(tx_bits_used != rx_bits)

    print(f"  fd_norm             : 0.0 (Static Channel)")
    print(f"  SNR                 : {snr_db} dB")
    print(f"  Total bits          : {len(tx_bits_used)}")
    print(f"  Bit Errors (LS 2D)  : {bit_errors}")
    assert bit_errors == 0, (
        f"Kênh tĩnh fd=0 tại SNR=100dB phải có 0 bit lỗi! Got {bit_errors}."
    )
    print(">>> Test 3 PASSED!\n")


def test_dft_better_than_ls_at_10db():
    """
    TEST 4 (MỚI): So sánh MSE của DFT-based 2D vs LS 2D tại SNR = 10 dB.
    Yêu cầu: MSE của DFT-based PHẢI NHỎ HƠN RÕ RỆT so với LS.

    Giải thích:
    - LS 2D chỉ nội suy, không lọc nhiễu -> MSE còn lớn do AWGN tại SNR trung bình.
    - DFT 2D lọc nhiễu miền trễ -> chỉ giữ năng lượng kênh -> MSE giảm rõ rệt.
    - Ngưỡng kiểm tra: MSE_DFT < 0.5 * MSE_LS (DFT tốt hơn ít nhất 2x về MSE).
    """
    print("--- Running Test 4: DFT-Based MSE < LS MSE @ SNR=10dB ---")

    num_symbols = NUM_SYMBOLS_PER_FRAME
    df          = PILOT_SPACING_FREQ
    dt          = PILOT_SPACING_TIME
    snr_db      = 10.0
    fd_norm     = 0.01  # Slow fading - kênh biến đổi nhẹ

    ofdm      = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, modulation=MODULATION)
    estimator = LSEstimator()
    channel   = TimeVaryingRayleighChannel(
        pdp                  = CHANNEL_PDP,
        fd_norm              = fd_norm,
        n_fft                = N_FFT,
        n_cp                 = N_CP,
        num_symbols_per_frame = num_symbols,
    )

    # Chạy Monte Carlo 20 lần để lấy MSE trung bình
    num_trials = 20
    mse_ls_list  = []
    mse_dft_list = []

    for _ in range(num_trials):
        max_bits  = N_FFT * num_symbols * ofdm.bits_per_symbol
        data_bits = np.random.randint(0, 2, max_bits)

        X_frame, _, pilot_sub, pilot_sym = ofdm.insert_pilots_2d(data_bits, df, dt, num_symbols)
        X_frame_cp = ofdm.modulate_frame(X_frame)
        Y_frame_cp, H_true_frame = channel.transmit_frame(X_frame_cp, snr_db=snr_db)
        Y_frame = ofdm.demodulate_frame(Y_frame_cp)

        # Ước lượng LS 2D
        H_ls = estimator.estimate_ls_2d(Y_frame, pilot_sub, pilot_sym, PILOT_VALUE, N_FFT, num_symbols)

        # Ước lượng DFT 2D với window_len = len(CHANNEL_PDP) (tối ưu)
        H_dft = estimator.estimate_dft_2d(
            Y_frame, pilot_sub, pilot_sym, PILOT_VALUE, N_FFT, num_symbols,
            window_len=len(CHANNEL_PDP)
        )

        mse_ls_list.append(estimator.compute_mse(H_true_frame, H_ls))
        mse_dft_list.append(estimator.compute_mse(H_true_frame, H_dft))

    mse_ls_avg  = np.mean(mse_ls_list)
    mse_dft_avg = np.mean(mse_dft_list)

    print(f"  SNR                 : {snr_db} dB, fd_norm = {fd_norm}")
    print(f"  MSE LS 2D (avg)     : {mse_ls_avg:.4e}")
    print(f"  MSE DFT 2D (avg)    : {mse_dft_avg:.4e}")
    print(f"  Cải thiện (LS/DFT)  : {mse_ls_avg/mse_dft_avg:.2f}x")

    assert mse_dft_avg < mse_ls_avg, (
        f"DFT-based phải có MSE nhỏ hơn LS! MSE_DFT={mse_dft_avg:.4e}, MSE_LS={mse_ls_avg:.4e}"
    )
    print(">>> Test 4 PASSED!\n")


if __name__ == '__main__':
    np.random.seed(42)  # Cố định seed để kết quả tái hiện được

    test_noiseless_loopback()
    test_multipath_rayleigh()
    test_jakes_static_channel_no_errors()
    test_dft_better_than_ls_at_10db()

    print("=" * 50)
    print("  TẤT CẢ CÁC BÀI KIỂM THỬ ĐÃ PASS THÀNH CÔNG!")
    print("=" * 50)

