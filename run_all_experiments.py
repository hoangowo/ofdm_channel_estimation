import os
import time
import csv
import matplotlib
# Cấu hình matplotlib sử dụng backend 'Agg' để vẽ biểu đồ và lưu file PNG ở chế độ không giao diện (headless)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Nhập cấu hình và các khối module đã xây dựng
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
                   num_symbols=NUM_OFDM_SYMBOLS, pdp=CHANNEL_PDP, include_mmse=True):
    """
    HÀM THỰC THI MÔ PHỎNG MONTE CARLO CHO HỆ THỐNG OFDM:
    Chạy thử nghiệm hệ thống qua toàn bộ dải SNR (0 đến 24 dB) để tính toán tỷ lệ BER và sai số MSE.
    """
    # Khởi tạo các đối tượng hệ thống, kênh truyền và bộ ước lượng kênh
    ofdm = OFDMSystem(n_fft=n_fft, n_cp=n_cp, pilot_spacing=pilot_spacing, 
                      modulation=modulation, pilot_value=PILOT_VALUE)
    channel = MultipathRayleighChannel(pdp)
    estimator = LSEstimator()
    
    # Tính số lượng bit dữ liệu truyền trong 1 OFDM symbol
    total_data_bits_per_symbol = ofdm.num_data * ofdm.bits_per_symbol
    
    # Các mảng tích lũy kết quả để phục vụ vẽ đồ thị
    ber_perfect = []
    ber_ls = []
    ber_mmse = []
    mse_ls = []
    mse_mmse = []
    
    # Vòng lặp quét ngoài: chạy qua từng mức SNR trong dải (dB)
    for snr_db in snr_range:
        total_bit_errors_perfect = 0
        total_bit_errors_ls = 0
        total_bit_errors_mmse = 0
        total_bits = 0
        accumulated_mse_ls = 0.0
        accumulated_mse_mmse = 0.0
        
        # Vòng lặp quét trong: Chạy mô phỏng qua hàng nghìn symbol OFDM độc lập (Monte Carlo)
        for sym_idx in range(num_symbols):
            # 1. KHỐI PHÁT (TRANSMITTER)
            # Sinh chuỗi bit ngẫu nhiên -> Ánh xạ QPSK/16QAM -> Chèn Pilot -> Điều chế IFFT + CP
            tx_bits = np.random.randint(0, 2, total_data_bits_per_symbol)
            tx_data_symbols = ofdm.map_bits(tx_bits)
            X_freq = ofdm.insert_pilots(tx_data_symbols)
            x_cp = ofdm.modulate(X_freq)
            
            # 2. KHỐI KÊNH TRUYỀN (CHANNEL)
            # Sinh kênh Rayleigh ngẫu nhiên -> Truyền tín hiệu -> Cộng nhiễu AWGN
            h_taps = channel.generate_taps()
            H_true = channel.get_frequency_response(h_taps, n_fft)
            y_cp = channel.apply_channel(x_cp, h_taps, snr_db)
            
            # 3. KHỐI THU (RECEIVER)
            # Giải điều chế FFT bỏ CP
            Y_freq = ofdm.demodulate(y_cp)
            
            # 4. KHỐI ƯỚC LƯỢNG KÊNH (CHANNEL ESTIMATION)
            # Ước lượng LS
            H_ls = estimator.estimate(Y_freq, ofdm.pilot_indices, PILOT_VALUE, n_fft)
            accumulated_mse_ls += estimator.compute_mse(H_true, H_ls)
            
            # Ước lượng MMSE (nếu cấu hình yêu cầu)
            if include_mmse:
                H_mmse = estimator.estimate_mmse(Y_freq, ofdm.pilot_indices, PILOT_VALUE, n_fft, snr_db)
                accumulated_mse_mmse += estimator.compute_mse(H_true, H_mmse)
            
            # 5. KHỐI CÂN BẰNG & GIẢI ĐIỀU CHẾ (EQUALIZATION & DEMAPPING)
            # 5a. Nhánh Perfect CSI (CSI hoàn hảo làm ngưỡng dưới tối ưu)
            rx_data_perfect = estimator.equalize(Y_freq[ofdm.data_indices], H_true[ofdm.data_indices])
            rx_bits_perfect = ofdm.demap_symbols(rx_data_perfect)
            total_bit_errors_perfect += np.sum(tx_bits != rx_bits_perfect)
            
            # 5b. Nhánh LS CSI
            rx_data_ls = estimator.equalize(Y_freq[ofdm.data_indices], H_ls[ofdm.data_indices])
            rx_bits_ls = ofdm.demap_symbols(rx_data_ls)
            total_bit_errors_ls += np.sum(tx_bits != rx_bits_ls)
            
            # 5c. Nhánh MMSE CSI (sử dụng MMSE Equalization thay vì ZF để tránh khuếch đại nhiễu)
            if include_mmse:
                rx_data_mmse = estimator.equalize_mmse(Y_freq[ofdm.data_indices], H_mmse[ofdm.data_indices], snr_db)
                rx_bits_mmse = ofdm.demap_symbols(rx_data_mmse)
                total_bit_errors_mmse += np.sum(tx_bits != rx_bits_mmse)
            
            total_bits += total_data_bits_per_symbol
            
        # Tính tỷ lệ BER trung bình và sai số kênh MSE trung bình của mức SNR hiện tại
        avg_ber_perfect = total_bit_errors_perfect / total_bits
        avg_ber_ls = total_bit_errors_ls / total_bits
        avg_mse_ls = accumulated_mse_ls / num_symbols
        
        ber_perfect.append(avg_ber_perfect)
        ber_ls.append(avg_ber_ls)
        mse_ls.append(avg_mse_ls)
        
        if include_mmse:
            avg_ber_mmse = total_bit_errors_mmse / total_bits
            avg_mse_mmse = accumulated_mse_mmse / num_symbols
            ber_mmse.append(avg_ber_mmse)
            mse_mmse.append(avg_mse_mmse)
            print(f"SNR = {snr_db:2d} dB | BER Perf: {avg_ber_perfect:.6f} | BER LS: {avg_ber_ls:.6f} | BER MMSE: {avg_ber_mmse:.6f} | MSE LS: {avg_mse_ls:.6e}")
        else:
            print(f"SNR = {snr_db:2d} dB | BER Perf: {avg_ber_perfect:.6f} | BER LS: {avg_ber_ls:.6f} | MSE LS: {avg_mse_ls:.6e}")
            
    return {
        'ber_perfect': np.array(ber_perfect),
        'ber_ls': np.array(ber_ls),
        'ber_mmse': np.array(ber_mmse) if include_mmse else None,
        'mse_ls': np.array(mse_ls),
        'mse_mmse': np.array(mse_mmse) if include_mmse else None,
    }

def export_csv_results(snr_range, main_res, tradeoff_res, csv_path):
    """
    XUẤT FILE KẾT QUẢ CSV:
    Ghi lại toàn bộ bảng số liệu BER/MSE theo SNR ra tệp simulation_results.csv.
    """
    with open(csv_path, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['SNR_dB', 'BER_Perfect_CSI', 'BER_LS_Df4', 'MSE_LS_Df4', 
                         'BER_MMSE_Df4', 'MSE_MMSE_Df4', 'BER_LS_Df2', 'MSE_LS_Df2', 
                         'BER_LS_Df8', 'MSE_LS_Df8'])
        
        for idx, snr in enumerate(snr_range):
            writer.writerow([
                snr,
                f"{main_res['ber_perfect'][idx]:.6e}",
                f"{main_res['ber_ls'][idx]:.6e}",
                f"{main_res['mse_ls'][idx]:.6e}",
                f"{main_res['ber_mmse'][idx]:.6e}",
                f"{main_res['mse_mmse'][idx]:.6e}",
                f"{tradeoff_res[2]['ber_ls'][idx]:.6e}",
                f"{tradeoff_res[2]['mse_ls'][idx]:.6e}",
                f"{tradeoff_res[8]['ber_ls'][idx]:.6e}",
                f"{tradeoff_res[8]['mse_ls'][idx]:.6e}",
            ])
    print(f"Successfully exported simulation results CSV: {csv_path}")

def generate_constellation_and_fit_plots(output_dir):
    """
    VẼ ĐỒ THỊ CHÒM SAO IQ VÀ SNAPSHOT ĐÁP ỨNG KÊNH:
    """
    ofdm = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, modulation=MODULATION)
    channel = MultipathRayleighChannel(CHANNEL_PDP)
    estimator = LSEstimator()
    
    snr_db = 15.0
    h_taps = channel.generate_taps()
    H_true = channel.get_frequency_response(h_taps, N_FFT)
    
    num_syms = 5
    tx_bits = np.random.randint(0, 2, ofdm.num_data * ofdm.bits_per_symbol * num_syms)
    tx_data_symbols = ofdm.map_bits(tx_bits)
    
    rx_raw_symbols_all = []
    rx_eq_symbols_all = []
    
    # Thu thập điểm chòm sao trước và sau khi cân bằng
    for i in range(num_syms):
        sym_chunk = tx_data_symbols[i * ofdm.num_data : (i + 1) * ofdm.num_data]
        X_freq = ofdm.insert_pilots(sym_chunk)
        x_cp = ofdm.modulate(X_freq)
        y_cp = channel.apply_channel(x_cp, h_taps, snr_db)
        Y_freq = ofdm.demodulate(y_cp)
        
        H_est_ls = estimator.estimate(Y_freq, ofdm.pilot_indices, PILOT_VALUE, N_FFT)
        rx_data_eq = estimator.equalize(Y_freq[ofdm.data_indices], H_est_ls[ofdm.data_indices])
        
        rx_raw_symbols_all.extend(Y_freq[ofdm.data_indices])
        rx_eq_symbols_all.extend(rx_data_eq)
        
    rx_raw_symbols_all = np.array(rx_raw_symbols_all)
    rx_eq_symbols_all = np.array(rx_eq_symbols_all)
    
    # 1. Vẽ Biểu đồ Chòm sao IQ (Constellation)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))
    
    # Trực quan hóa chòm sao bị nhiễu và méo (trước cân bằng)
    ax1.scatter(np.real(rx_raw_symbols_all), np.imag(rx_raw_symbols_all), color='crimson', alpha=0.6, label='Received (Before Equalization)')
    ax1.scatter(np.real(tx_data_symbols[:len(rx_raw_symbols_all)]), np.imag(tx_data_symbols[:len(rx_raw_symbols_all)]), 
                color='black', marker='+', s=100, linewidth=2, label='Transmitted Ideal')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.set_xlabel('In-Phase (I)', fontsize=12)
    ax1.set_ylabel('Quadrature (Q)', fontsize=12)
    ax1.set_title(f'Constellation Before Equalization (SNR={snr_db} dB)', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.axis('equal')
    
    # Trực quan hóa chòm sao co cụm lại (sau cân bằng ZF)
    ax2.scatter(np.real(rx_eq_symbols_all), np.imag(rx_eq_symbols_all), color='royalblue', alpha=0.7, label='Received (After ZF Equalization)')
    ax2.scatter(np.real(tx_data_symbols[:len(rx_eq_symbols_all)]), np.imag(tx_data_symbols[:len(rx_eq_symbols_all)]), 
                color='black', marker='+', s=100, linewidth=2, label='Transmitted Ideal')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.set_xlabel('In-Phase (I)', fontsize=12)
    ax2.set_ylabel('Quadrature (Q)', fontsize=12)
    ax2.set_title(f'Constellation After ZF Equalization (SNR={snr_db} dB)', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.axis('equal')
    
    plt.tight_layout()
    const_path = os.path.join(output_dir, 'constellation_comparison.png')
    plt.savefig(const_path, dpi=300)
    plt.close()
    print(f"Saved figure: {const_path}")

    # 2. Vẽ đồ thị so sánh đáp ứng tần số kênh truyền |H(k)|
    snr_snap = 10.0
    h_taps_snap = channel.generate_taps()
    H_true_snap = channel.get_frequency_response(h_taps_snap, N_FFT)
    tx_bits_snap = np.random.randint(0, 2, ofdm.num_data * ofdm.bits_per_symbol)
    X_freq_snap = ofdm.insert_pilots(ofdm.map_bits(tx_bits_snap))
    x_cp_snap = ofdm.modulate(X_freq_snap)
    y_cp_snap = channel.apply_channel(x_cp_snap, h_taps_snap, snr_snap)
    Y_freq_snap = ofdm.demodulate(y_cp_snap)
    
    H_ls_snap = estimator.estimate(Y_freq_snap, ofdm.pilot_indices, PILOT_VALUE, N_FFT)
    H_mmse_snap = estimator.estimate_mmse(Y_freq_snap, ofdm.pilot_indices, PILOT_VALUE, N_FFT, snr_snap)
    
    subcarriers = np.arange(N_FFT)
    plt.figure(figsize=(10, 5))
    plt.plot(subcarriers, np.abs(H_true_snap), 'k-', label='True Channel |H_true(k)|', linewidth=2.5)
    plt.plot(subcarriers, np.abs(H_ls_snap), 'r--', label='LS Estimate |H_LS(k)|', linewidth=1.8)
    plt.plot(subcarriers, np.abs(H_mmse_snap), 'g-.', label='MMSE Estimate |H_MMSE(k)|', linewidth=1.8)
    plt.scatter(ofdm.pilot_indices, np.abs(Y_freq_snap[ofdm.pilot_indices]/PILOT_VALUE), color='red', s=50, zorder=5, label='Pilot Measurements')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlabel('Subcarrier Index (k)', fontsize=12)
    plt.ylabel('Magnitude |H(k)|', fontsize=12)
    plt.title(f'Channel Frequency Response Estimation Snapshot (SNR={snr_snap} dB)', fontsize=13)
    plt.legend(fontsize=11)
    plt.tight_layout()
    fit_path = os.path.join(output_dir, 'channel_estimation_fit.png')
    plt.savefig(fit_path, dpi=300)
    plt.close()
    print(f"Saved figure: {fit_path}")

def generate_all_plots(snr_range, main_res, tradeoff_res, output_dir):
    """
    TỔNG HỢP VẼ CẢ 5 ĐỒ THỊ KẾT QUẢ MÔ PHỎNG PHỤC VỤ BÁO CÁO:
    """
    # 1. Đồ thị BER vs SNR (Perfect CSI vs LS vs MMSE)
    plt.figure(figsize=(8, 6))
    plt.semilogy(snr_range, main_res['ber_perfect'], 'b-o', label='Perfect CSI (Lower Bound)', linewidth=2, markersize=6)
    plt.semilogy(snr_range, main_res['ber_ls'], 'r-s', label=f'LS Estimation (Df={PILOT_SPACING})', linewidth=2, markersize=6)
    plt.semilogy(snr_range, main_res['ber_mmse'], 'g-^', label=f'MMSE Estimation (Df={PILOT_SPACING})', linewidth=2, markersize=6)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.xlabel('SNR (dB)', fontsize=12)
    plt.ylabel('Bit Error Rate (BER)', fontsize=12)
    plt.title('OFDM System Performance: BER vs SNR under Rayleigh Fading', fontsize=13)
    plt.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ber_vs_snr.png'), dpi=300)
    plt.close()

    # 2. Đồ thị MSE vs SNR (LS vs MMSE)
    plt.figure(figsize=(8, 6))
    plt.semilogy(snr_range, main_res['mse_ls'], 'r-s', label=f'LS Channel Estimation MSE (Df={PILOT_SPACING})', linewidth=2, markersize=6)
    plt.semilogy(snr_range, main_res['mse_mmse'], 'g-^', label=f'MMSE Channel Estimation MSE (Df={PILOT_SPACING})', linewidth=2, markersize=6)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.xlabel('SNR (dB)', fontsize=12)
    plt.ylabel('Mean Squared Error (MSE)', fontsize=12)
    plt.title('Channel Estimation Performance: MSE vs SNR', fontsize=13)
    plt.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mse_vs_snr.png'), dpi=300)
    plt.close()

    # 3. Đồ thị Phân tích Đánh đổi Pilot spacing (D_f = 2, 4, 8)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    colors = {2: 'purple', 4: 'red', 8: 'orange'}
    markers = {2: 'o', 4: 's', 8: '^'}
    
    ax1.semilogy(snr_range, tradeoff_res[4]['ber_perfect'], 'k--', label='Perfect CSI', linewidth=2)
    for df in [2, 4, 8]:
        overhead = (N_FFT // df) / N_FFT * 100
        ax1.semilogy(snr_range, tradeoff_res[df]['ber_ls'], color=colors[df], marker=markers[df],
                     label=f'LS D_f={df} (Overhead: {overhead:.1f}%)', linewidth=2, markersize=5)
    ax1.grid(True, which='both', linestyle='--', alpha=0.7)
    ax1.set_xlabel('SNR (dB)', fontsize=12)
    ax1.set_ylabel('Bit Error Rate (BER)', fontsize=12)
    ax1.set_title('BER Comparison for Different Pilot Spacings (D_f)', fontsize=12)
    ax1.legend(fontsize=10)
    
    for df in [2, 4, 8]:
        ax2.semilogy(snr_range, tradeoff_res[df]['mse_ls'], color=colors[df], marker=markers[df],
                     label=f'MSE D_f={df}', linewidth=2, markersize=5)
    ax2.grid(True, which='both', linestyle='--', alpha=0.7)
    ax2.set_xlabel('SNR (dB)', fontsize=12)
    ax2.set_ylabel('Mean Squared Error (MSE)', fontsize=12)
    ax2.set_title('MSE Comparison for Different Pilot Spacings (D_f)', fontsize=12)
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pilot_spacing_comparison.png'), dpi=300)
    plt.close()

    # 4. Tạo nốt chòm sao & Fit snapshot
    generate_constellation_and_fit_plots(output_dir)
    
    print("All 5 figures successfully generated and saved!")

def main():
    print("=======================================================")
    print("RUN ALL EXPERIMENTS: OFDM LS & MMSE Channel Estimation")
    print("=======================================================\n")
    
    start_time = time.time()
    
    # Đặt seed ngẫu nhiên cố định để đảm bảo kết quả tái tạo được giữa các lần chạy
    np.random.seed(42)
    
    # Thí nghiệm 1: Chạy mô phỏng chính cho cả LS & MMSE với D_f = 4 (2000 Symbols)
    print("--- [Experiment 1] Main OFDM Simulation (Df = 4, LS & MMSE) ---")
    main_results = run_simulation(num_symbols=2000, include_mmse=True)
    
    # Thí nghiệm 2: Chạy đánh đổi khoảng cách pilot (D_f = 2, 8) so sánh hiệu năng
    print("\n--- [Experiment 2] Pilot Spacing Trade-off (Df = 2, 4, 8) ---")
    tradeoff_results = {4: main_results}
    for df in [2, 8]:
        print(f"\n--- Simulating Pilot Spacing D_f = {df} ---")
        tradeoff_results[df] = run_simulation(pilot_spacing=df, num_symbols=1000, include_mmse=False)
        
    # Ghi xuất tệp dữ liệu CSV & vẽ toàn bộ biểu đồ báo cáo
    csv_path = os.path.join(SCRIPT_DIR, 'simulation_results.csv')
    export_csv_results(SNR_DB_RANGE, main_results, tradeoff_results, csv_path)
    generate_all_plots(SNR_DB_RANGE, main_results, tradeoff_results, SCRIPT_DIR)
    
    elapsed = time.time() - start_time
    print(f"\n=======================================================")
    print(f"All experiments completed successfully in {elapsed:.2f} seconds!")
    print("=======================================================")

if __name__ == '__main__':
    main()
