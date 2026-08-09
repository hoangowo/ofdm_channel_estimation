import os
import time
import matplotlib
# Sử dụng backend 'Agg' để vẽ đồ thị tĩnh và lưu ảnh PNG không cần môi trường hiển thị đồ họa (headless)
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
import numpy as np

# Nhập các thông số cấu hình và các khối lớp mô phỏng
from config import (
    N_FFT, N_CP, MODULATION, PILOT_SPACING, PILOT_VALUE,
    SNR_DB_RANGE, NUM_OFDM_SYMBOLS, CHANNEL_PDP
)
from channel import MultipathRayleighChannel
from ofdm_system import OFDMSystem
from estimator import LSEstimator

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_simulation(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, 
                   modulation=MODULATION, snr_range=SNR_DB_RANGE, 
                   num_symbols=NUM_OFDM_SYMBOLS, pdp=CHANNEL_PDP):
    """
    HÀM CHẠY MÔ PHỎNG OFDM (CHỈ SỬ DỤNG THUẬT TOÁN LS):
    Tính toán tỷ lệ BER của Perfect CSI vs LS CSI, và MSE kênh truyền LS theo từng mức SNR.
    """
    # Khởi tạo các module
    ofdm = OFDMSystem(n_fft=n_fft, n_cp=n_cp, pilot_spacing=pilot_spacing, 
                      modulation=modulation, pilot_value=PILOT_VALUE)
    channel = MultipathRayleighChannel(pdp)
    estimator = LSEstimator()
    
    # Số bit truyền trên một symbol OFDM
    total_data_bits_per_symbol = ofdm.num_data * ofdm.bits_per_symbol
    
    ber_perfect = []
    ber_ls = []
    mse_ls = []
    
    # Vòng lặp quét SNR
    for snr_db in snr_range:
        total_bit_errors_perfect = 0
        total_bit_errors_ls = 0
        total_bits = 0
        accumulated_mse = 0.0
        
        # Vòng lặp Monte Carlo sinh các symbol OFDM
        for sym_idx in range(num_symbols):
            # 1. Khối Phát (Transmitter): sinh bit -> điều chế QPSK/16QAM -> chèn pilot -> IFFT + CP
            tx_bits = np.random.randint(0, 2, total_data_bits_per_symbol)
            tx_data_symbols = ofdm.map_bits(tx_bits)
            X_freq = ofdm.insert_pilots(tx_data_symbols)
            x_cp = ofdm.modulate(X_freq)
            
            # 2. Khối Kênh truyền: Đi qua kênh Rayleigh đa đường + cộng nhiễu AWGN
            h_taps = channel.generate_taps()
            H_true = channel.get_frequency_response(h_taps, n_fft)
            y_cp = channel.apply_channel(x_cp, h_taps, snr_db)
            
            # 3. Khối Nhận: Tháo bỏ CP + Thực hiện FFT đưa về miền tần số
            Y_freq = ofdm.demodulate(y_cp)
            
            # 4. Khối Ước lượng LS: chia Y_pilot / X_pilot và nội suy tuyến tính
            H_est = estimator.estimate(Y_freq, ofdm.pilot_indices, PILOT_VALUE, n_fft)
            accumulated_mse += estimator.compute_mse(H_true, H_est)
            
            # 5a. Khối thu Perfect CSI: Cân bằng kênh lý tưởng ZF và đếm lỗi bit
            rx_data_perfect = estimator.equalize(Y_freq[ofdm.data_indices], H_true[ofdm.data_indices])
            rx_bits_perfect = ofdm.demap_symbols(rx_data_perfect)
            total_bit_errors_perfect += np.sum(tx_bits != rx_bits_perfect)
            
            # 5b. Khối thu LS CSI: Cân bằng kênh ZF bằng H_est và đếm lỗi bit
            rx_data_ls = estimator.equalize(Y_freq[ofdm.data_indices], H_est[ofdm.data_indices])
            rx_bits_ls = ofdm.demap_symbols(rx_data_ls)
            total_bit_errors_ls += np.sum(tx_bits != rx_bits_ls)
            
            total_bits += total_data_bits_per_symbol
            
        # Tính kết quả trung bình cho mức SNR này
        avg_ber_perfect = total_bit_errors_perfect / total_bits
        avg_ber_ls = total_bit_errors_ls / total_bits
        avg_mse_ls = accumulated_mse / num_symbols
        
        ber_perfect.append(avg_ber_perfect)
        ber_ls.append(avg_ber_ls)
        mse_ls.append(avg_mse_ls)
        
        print(f"SNR = {snr_db:2d} dB | BER Perfect CSI: {avg_ber_perfect:.6f} | BER LS Est: {avg_ber_ls:.6f} | MSE LS: {avg_mse_ls:.6e}")
        
    return {
        'ber_perfect': np.array(ber_perfect),
        'ber_ls': np.array(ber_ls),
        'mse_ls': np.array(mse_ls)
    }

def plot_main_results(snr_range, results, output_dir=SCRIPT_DIR):
    """
    Vẽ đồ thị BER vs SNR và MSE vs SNR cho thuật toán LS.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Đồ thị BER vs SNR
    plt.figure(figsize=(8, 6))
    plt.semilogy(snr_range, results['ber_perfect'], 'b-o', label='Perfect CSI (Lower Bound)', linewidth=2, markersize=6)
    plt.semilogy(snr_range, results['ber_ls'], 'r-s', label=f'LS Estimation (Df={PILOT_SPACING})', linewidth=2, markersize=6)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.xlabel('SNR (dB)', fontsize=12)
    plt.ylabel('Bit Error Rate (BER)', fontsize=12)
    plt.title('OFDM System Performance: BER vs SNR under Rayleigh Fading', fontsize=13)
    plt.legend(fontsize=11)
    plt.tight_layout()
    ber_path = os.path.join(output_dir, 'ber_vs_snr.png')
    plt.savefig(ber_path, dpi=300)
    plt.close()
    print(f"Saved figure: {ber_path}")
    
    # 2. Đồ thị MSE vs SNR
    plt.figure(figsize=(8, 6))
    plt.semilogy(snr_range, results['mse_ls'], 'g-^', label=f'LS Channel Estimation MSE (Df={PILOT_SPACING})', linewidth=2, markersize=6)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.xlabel('SNR (dB)', fontsize=12)
    plt.ylabel('Mean Squared Error (MSE)', fontsize=12)
    plt.title('LS Channel Estimation Performance: MSE vs SNR', fontsize=13)
    plt.legend(fontsize=11)
    plt.tight_layout()
    mse_path = os.path.join(output_dir, 'mse_vs_snr.png')
    plt.savefig(mse_path, dpi=300)
    plt.close()
    print(f"Saved figure: {mse_path}")

def run_pilot_spacing_tradeoff(output_dir=SCRIPT_DIR):
    """
    Nhiệm vụ 7: Đánh giá sự đánh đổi khoảng cách pilot D_f = 2, 4, 8.
    """
    print("\n=======================================================")
    print("Task 7: Running Pilot Spacing Trade-off Experiments (Df = 2, 4, 8)")
    print("=======================================================")
    
    spacings = [2, 4, 8]
    spacing_results = {}
    
    # Sử dụng 1000 symbols cho thử nghiệm trade-off
    tradeoff_symbols = 1000
    
    for df in spacings:
        print(f"\n--- Simulating Pilot Spacing D_f = {df} ---")
        res = run_simulation(pilot_spacing=df, num_symbols=tradeoff_symbols)
        spacing_results[df] = res
        
    # Tạo đồ thị so sánh 2 nhánh: BER và MSE
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    colors = {2: 'purple', 4: 'red', 8: 'orange'}
    markers = {2: 'o', 4: 's', 8: '^'}
    
    # BER plot
    ax1.semilogy(SNR_DB_RANGE, spacing_results[4]['ber_perfect'], 'k--', label='Perfect CSI', linewidth=2)
    for df in spacings:
        overhead = (N_FFT // df) / N_FFT * 100
        ax1.semilogy(SNR_DB_RANGE, spacing_results[df]['ber_ls'], color=colors[df], marker=markers[df],
                     label=f'LS D_f={df} (Overhead: {overhead:.1f}%)', linewidth=2, markersize=5)
    ax1.grid(True, which='both', linestyle='--', alpha=0.7)
    ax1.set_xlabel('SNR (dB)', fontsize=12)
    ax1.set_ylabel('Bit Error Rate (BER)', fontsize=12)
    ax1.set_title('BER Comparison for Different Pilot Spacings (D_f)', fontsize=12)
    ax1.legend(fontsize=10)
    
    # MSE plot
    for df in spacings:
        ax2.semilogy(SNR_DB_RANGE, spacing_results[df]['mse_ls'], color=colors[df], marker=markers[df],
                     label=f'MSE D_f={df}', linewidth=2, markersize=5)
    ax2.grid(True, which='both', linestyle='--', alpha=0.7)
    ax2.set_xlabel('SNR (dB)', fontsize=12)
    ax2.set_ylabel('Mean Squared Error (MSE)', fontsize=12)
    ax2.set_title('MSE Comparison for Different Pilot Spacings (D_f)', fontsize=12)
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    comp_path = os.path.join(output_dir, 'pilot_spacing_comparison.png')
    plt.savefig(comp_path, dpi=300)
    plt.close()
    print(f"\nSaved pilot spacing trade-off comparison figure: {comp_path}")

def main():
    print("=======================================================")
    print("OFDM Channel Estimation Simulation (LS Algorithm)")
    print(f"Parameters: N_FFT={N_FFT}, N_CP={N_CP}, Mod={MODULATION}, D_f={PILOT_SPACING}, Symbols={NUM_OFDM_SYMBOLS}")
    print("=======================================================\n")
    
    start_time = time.time()
    
    # 1. Chạy Mô phỏng Chính (Mặc định D_f = 4)
    print("--- Running Main Simulation (D_f = 4) ---")
    results = run_simulation()
    plot_main_results(SNR_DB_RANGE, results)
    
    # 2. Khảo sát đánh đổi khoảng cách pilot (D_f = 2, 4, 8)
    run_pilot_spacing_tradeoff()
    
    elapsed = time.time() - start_time
    print(f"\n=======================================================")
    print(f"Simulation completed successfully in {elapsed:.2f} seconds!")
    print("=======================================================")

if __name__ == '__main__':
    main()
