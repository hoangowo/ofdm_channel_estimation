"""
run_thesis_experiments.py - Script tự động chạy toàn bộ thí nghiệm đồ án.
Đồ án: "Pilot-Based Channel Estimation for OFDM over Time-Varying Fading Channels"
Sinh viên: Đào Kim Hoàng - MSSV: 2211068

Thí nghiệm 1: Baseline Check - So sánh Perfect CSI, LS, DFT-based @ kênh tĩnh/chậm.
Thí nghiệm 2: Ảnh hưởng Mobility/Doppler - Quét fd_norm từ 0 đến 0.10.
Thí nghiệm 3: Khảo sát Lưới Pilot 2D (Df, Dt) Trade-off.
Thí nghiệm 4: Độ phức tạp & Bảng Khuyến nghị Pilot.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Backend không cần GUI (chạy headless)
import matplotlib.pyplot as plt
import csv
import time
import warnings
warnings.filterwarnings('ignore')

from config import (N_FFT, N_CP, CHANNEL_PDP, PILOT_VALUE, SNR_DB_RANGE,
                    NUM_SYMBOLS_PER_FRAME, NUM_FRAMES,
                    PILOT_SPACING_FREQ, PILOT_SPACING_TIME,
                    FD_NORM_LIST, DF_LIST, DT_LIST)
from channel import TimeVaryingRayleighChannel
from ofdm_system import OFDMSystem
from estimator import LSEstimator

# ===========================================================================
# HÀM TIỆN ÍCH CHUNG
# ===========================================================================

def safe_semilogy(ax, x, y, *args, **kwargs):
    """
    [FIX 4b] Vẽ semilogy an toàn: BER=0 được thay bằng NaN để không hiển thị
    trên đồ thị log (thay vì giả mạo thành 1e-6).
    Ghi chú: BER=0 có nghĩa là không quan sát được lỗi nào, KHÔNG phải BER thực sự bằng 0.
    """
    y_safe = np.where(np.array(y) == 0, np.nan, y)
    return ax.semilogy(x, y_safe, *args, **kwargs)


# ===========================================================================

def simulate_frame(ofdm, estimator, channel, snr_db, df, dt, num_symbols,
                   window_len_ncp, window_len_pdp):
    """
    Mô phỏng 1 frame OFDM qua kênh, trả về MSE và BER cho 4 phương pháp:
    Perfect CSI, LS 2D, DFT-based (window=N_CP), DFT-based (window=len(PDP)).

    Tham số:
    - ofdm           : Instance OFDMSystem.
    - estimator      : Instance LSEstimator.
    - channel        : Instance TimeVaryingRayleighChannel.
    - snr_db         : SNR dB.
    - df, dt         : Cấu hình lưới pilot.
    - num_symbols    : Số symbol/frame.
    - window_len_ncp : Độ dài cửa sổ DFT = N_CP.
    - window_len_pdp : Độ dài cửa sổ DFT = len(PDP).

    Trả về:
    - dict gồm các MSE và số bit lỗi/tổng bit cho từng phương pháp.
    """
    # Sinh bit và tạo frame phát
    max_bits  = N_FFT * num_symbols * ofdm.bits_per_symbol
    data_bits = np.random.randint(0, 2, max_bits)
    X_frame, tx_bits, pilot_sub, pilot_sym = ofdm.insert_pilots_2d(data_bits, df, dt, num_symbols)

    # Điều chế + Truyền qua kênh
    X_frame_cp            = ofdm.modulate_frame(X_frame)
    Y_frame_cp, H_true    = channel.transmit_frame(X_frame_cp, snr_db=snr_db, mode='block_aging')
    Y_frame               = ofdm.demodulate_frame(Y_frame_cp)

    n_bits = len(tx_bits)

    # Perfect CSI (cân bằng với H_true)
    rx_perfect  = ofdm.equalize_and_decode_frame(Y_frame, H_true, pilot_sub, pilot_sym)
    err_perfect = int(np.sum(tx_bits != rx_perfect))

    # LS 2D
    H_ls       = estimator.estimate_ls_2d(Y_frame, pilot_sub, pilot_sym, PILOT_VALUE, N_FFT, num_symbols)
    rx_ls      = ofdm.equalize_and_decode_frame(Y_frame, H_ls, pilot_sub, pilot_sym)
    err_ls     = int(np.sum(tx_bits != rx_ls))
    mse_ls     = estimator.compute_mse(H_true, H_ls)

    # DFT-based (window = N_CP)
    H_dft_ncp  = estimator.estimate_dft_2d(Y_frame, pilot_sub, pilot_sym, PILOT_VALUE,
                                            N_FFT, num_symbols, window_len=window_len_ncp)
    rx_dft_ncp = ofdm.equalize_and_decode_frame(Y_frame, H_dft_ncp, pilot_sub, pilot_sym)
    err_dft_ncp= int(np.sum(tx_bits != rx_dft_ncp))
    mse_dft_ncp= estimator.compute_mse(H_true, H_dft_ncp)

    # DFT-based (window = len(PDP))
    H_dft_pdp  = estimator.estimate_dft_2d(Y_frame, pilot_sub, pilot_sym, PILOT_VALUE,
                                            N_FFT, num_symbols, window_len=window_len_pdp)
    rx_dft_pdp = ofdm.equalize_and_decode_frame(Y_frame, H_dft_pdp, pilot_sub, pilot_sym)
    err_dft_pdp= int(np.sum(tx_bits != rx_dft_pdp))
    mse_dft_pdp= estimator.compute_mse(H_true, H_dft_pdp)

    return {
        'n_bits'       : n_bits,
        'err_perfect'  : err_perfect,
        'err_ls'       : err_ls,
        'err_dft_ncp'  : err_dft_ncp,
        'err_dft_pdp'  : err_dft_pdp,
        'mse_ls'       : mse_ls,
        'mse_dft_ncp'  : mse_dft_ncp,
        'mse_dft_pdp'  : mse_dft_pdp,
    }


def run_snr_sweep(fd_norm, df, dt, modulation, num_frames, num_symbols,
                  snr_range, window_len_ncp, window_len_pdp):
    """
    Quét toàn bộ dải SNR, chạy Monte Carlo num_frames frame mỗi điểm SNR.
    Trả về dict các mảng MSE và BER theo SNR.
    """
    ofdm      = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING_FREQ,
                           modulation=modulation, pilot_value=PILOT_VALUE)
    estimator = LSEstimator()
    # [FIX 3] Channel được tạo mới cho mỗi frame trong vòng lặp Monte Carlo bên dưới
    # để đảm bảo tính độc lập thống kê giữa các frame

    mse_ls_arr     = np.zeros(len(snr_range))
    mse_dft_ncp_arr= np.zeros(len(snr_range))
    mse_dft_pdp_arr= np.zeros(len(snr_range))
    ber_perfect_arr= np.zeros(len(snr_range))
    ber_ls_arr     = np.zeros(len(snr_range))
    ber_dft_ncp_arr= np.zeros(len(snr_range))
    ber_dft_pdp_arr= np.zeros(len(snr_range))

    for si, snr_db in enumerate(snr_range):
        total_bits = {'n':0,'ep':0,'els':0,'edcp':0,'edpd':0,
                      'mls':0,'mdcp':0,'mdpd':0}
        for _ in range(num_frames):
            # [FIX 3] Tạo channel mới cho mỗi frame Monte Carlo
            # để đảm bảo các frame có channel realization độc lập
            # (tránh lặp lại cùng channel trajectory do reuse phase offsets)
            channel = TimeVaryingRayleighChannel(
                pdp=CHANNEL_PDP, fd_norm=fd_norm, n_fft=N_FFT, n_cp=N_CP,
                num_symbols_per_frame=num_symbols
            )
            r = simulate_frame(ofdm, estimator, channel, snr_db, df, dt, num_symbols,
                               window_len_ncp, window_len_pdp)
            total_bits['n']    += r['n_bits']
            total_bits['ep']   += r['err_perfect']
            total_bits['els']  += r['err_ls']
            total_bits['edcp'] += r['err_dft_ncp']
            total_bits['edpd'] += r['err_dft_pdp']
            total_bits['mls']  += r['mse_ls']
            total_bits['mdcp'] += r['mse_dft_ncp']
            total_bits['mdpd'] += r['mse_dft_pdp']

        n = total_bits['n']
        mse_ls_arr[si]      = total_bits['mls'] / num_frames
        mse_dft_ncp_arr[si] = total_bits['mdcp'] / num_frames
        mse_dft_pdp_arr[si] = total_bits['mdpd'] / num_frames
        # [FIX 4] Lưu BER thực tế, KHÔNG thay BER=0 bằng 1e-6
        # Ghi nhận số bit lỗi thực tế để báo cáo trung thực
        ber_perfect_arr[si] = total_bits['ep'] / n
        ber_ls_arr[si]      = total_bits['els'] / n
        ber_dft_ncp_arr[si] = total_bits['edcp'] / n
        ber_dft_pdp_arr[si] = total_bits['edpd'] / n

    return {
        'snr'        : snr_range,
        'mse_ls'     : mse_ls_arr,
        'mse_dft_ncp': mse_dft_ncp_arr,
        'mse_dft_pdp': mse_dft_pdp_arr,
        'ber_perfect': ber_perfect_arr,
        'ber_ls'     : ber_ls_arr,
        'ber_dft_ncp': ber_dft_ncp_arr,
        'ber_dft_pdp': ber_dft_pdp_arr,
    }


# ===========================================================================
# THÍ NGHIỆM 1: BASELINE CHECK
# ===========================================================================

def experiment1_baseline():
    """
    THÍ NGHIỆM 1: So sánh kênh tĩnh và slow fading giữa các phương pháp ước lượng.
    - Kịch bản: fd_norm = 0.0 (Static) và fd_norm = 0.005 (Pedestrian chậm).
    - Cấu hình: Dt=1, Df=4 (Pilot dày nhất theo thời gian).
    - Điều chế: QPSK và 16QAM.
    - Đầu ra: fig1_mse_vs_snr.png, fig2_ber_vs_snr.png.
    """
    print("\n" + "=" * 60)
    print("THÍ NGHIỆM 1: BASELINE CHECK - MSE & BER vs SNR")
    print("=" * 60)
    t0 = time.time()

    snr_range      = SNR_DB_RANGE
    num_symbols    = NUM_SYMBOLS_PER_FRAME
    num_frames     = NUM_FRAMES
    df, dt         = PILOT_SPACING_FREQ, 1  # Dt=1: mọi symbol đều có pilot
    window_ncp     = N_CP
    window_pdp     = len(CHANNEL_PDP)

    results_fd0_qpsk  = run_snr_sweep(0.0,   df, dt, 'QPSK',  num_frames, num_symbols,
                                      snr_range, window_ncp, window_pdp)
    results_fd5_qpsk  = run_snr_sweep(0.005, df, dt, 'QPSK',  num_frames, num_symbols,
                                      snr_range, window_ncp, window_pdp)
    results_fd0_16qam = run_snr_sweep(0.0,   df, dt, '16QAM', num_frames, num_symbols,
                                      snr_range, window_ncp, window_pdp)
    results_fd5_16qam = run_snr_sweep(0.005, df, dt, '16QAM', num_frames, num_symbols,
                                      snr_range, window_ncp, window_pdp)

    # --- VẼ ĐỒ THỊ 1: MSE vs SNR ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Thí nghiệm 1: MSE vs SNR - Baseline Check', fontsize=14, fontweight='bold')

    for ax, results, mod_name in [
        (axes[0], results_fd0_qpsk,  'QPSK  (fd=0.0)'),
        (axes[1], results_fd0_16qam, '16QAM (fd=0.0)'),
    ]:
        ax.semilogy(snr_range, results['mse_ls'],      'b-o', label='LS 2D', markersize=5)
        ax.semilogy(snr_range, results['mse_dft_ncp'], 'r-s', label=f'DFT (L={N_CP})', markersize=5)
        ax.semilogy(snr_range, results['mse_dft_pdp'], 'g-^', label=f'DFT (L={window_pdp})', markersize=5)
        ax.set_xlabel('SNR (dB)')
        ax.set_ylabel('MSE')
        ax.set_title(mod_name)
        ax.legend(fontsize=9)
        ax.grid(True, which='both', alpha=0.4)
        ax.set_xlim([snr_range[0], snr_range[-1]])

    plt.tight_layout()
    plt.savefig('fig1_mse_vs_snr.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  [OK] fig1_mse_vs_snr.png da luu.")

    # --- VẼ ĐỒ THỊ 2: BER vs SNR ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Thí nghiệm 1: BER vs SNR - Baseline Check', fontsize=14, fontweight='bold')

    for ax, (r0, r5, mod_name) in [
        (axes[0], (results_fd0_qpsk,  results_fd5_qpsk,  'QPSK')),
        (axes[1], (results_fd0_16qam, results_fd5_16qam, '16QAM')),
    ]:
        safe_semilogy(ax, snr_range, r0['ber_perfect'], 'k--',  label='Perfect CSI (fd=0)', linewidth=2)
        safe_semilogy(ax, snr_range, r0['ber_ls'],      'b-o',  label='LS 2D (fd=0)', markersize=5)
        safe_semilogy(ax, snr_range, r0['ber_dft_pdp'], 'g-^',  label=f'DFT L={window_pdp} (fd=0)', markersize=5)
        safe_semilogy(ax, snr_range, r5['ber_dft_pdp'], 'g--^', label=f'DFT L={window_pdp} (fd=0.005)', markersize=5, alpha=0.7)
        ax.set_xlabel('SNR (dB)')
        ax.set_ylabel('BER')
        ax.set_title(mod_name)
        ax.legend(fontsize=8)
        ax.grid(True, which='both', alpha=0.4)
        ax.set_xlim([snr_range[0], snr_range[-1]])
        ax.set_ylim([1e-4, 1])

    plt.tight_layout()
    plt.savefig('fig2_ber_vs_snr.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  [OK] fig2_ber_vs_snr.png da luu.")

    elapsed = time.time() - t0
    print(f"  Thoi gian thi nghiem 1: {elapsed:.1f}s")

    # Xuất CSV
    with open('thesis_results_snr.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['modulation', 'fd_norm', 'snr_db',
                         'mse_ls', 'mse_dft_ncp', 'mse_dft_pdp',
                         'ber_perfect', 'ber_ls', 'ber_dft_ncp', 'ber_dft_pdp'])
        for mod, fd, res in [
            ('QPSK',  0.0,   results_fd0_qpsk),
            ('QPSK',  0.005, results_fd5_qpsk),
            ('16QAM', 0.0,   results_fd0_16qam),
            ('16QAM', 0.005, results_fd5_16qam),
        ]:
            for i, snr in enumerate(snr_range):
                writer.writerow([mod, fd, snr,
                                 res['mse_ls'][i], res['mse_dft_ncp'][i], res['mse_dft_pdp'][i],
                                 res['ber_perfect'][i], res['ber_ls'][i],
                                 res['ber_dft_ncp'][i], res['ber_dft_pdp'][i]])
    print("  [OK] thesis_results_snr.csv da luu.")
    return results_fd0_qpsk


# ===========================================================================
# THÍ NGHIỆM 2: ẢNH HƯỞNG CỦA MOBILITY / DOPPLER
# ===========================================================================

def experiment2_mobility():
    """
    THÍ NGHIỆM 2: Quét fd_norm từ 0.0 đến 0.10 tại SNR cố định.
    - Cấu hình: Df=4, Dt=2.
    - SNR = 10 dB và 20 dB.
    - Đầu ra: fig3_mobility_impact_mse_ber.png, thesis_results_doppler.csv.
    """
    print("\n" + "=" * 60)
    print("THÍ NGHIỆM 2: ANH HUONG MOBILITY / DOPPLER")
    print("=" * 60)
    t0 = time.time()

    fd_list    = FD_NORM_LIST
    snr_list   = [10.0, 20.0]  # 2 mức SNR cố định
    df, dt     = PILOT_SPACING_FREQ, PILOT_SPACING_TIME
    num_symbols= NUM_SYMBOLS_PER_FRAME
    num_frames = NUM_FRAMES
    window_ncp = N_CP
    window_pdp = len(CHANNEL_PDP)
    snr_range  = SNR_DB_RANGE

    results_doppler = {}  # key: (snr_db, fd_norm)

    for snr_db in snr_list:
        for fd_norm in fd_list:
            print(f"  Dang chay: SNR={snr_db} dB, fd_norm={fd_norm}...")
            res = run_snr_sweep(fd_norm, df, dt, 'QPSK', num_frames, num_symbols,
                                np.array([snr_db]), window_ncp, window_pdp)
            results_doppler[(snr_db, fd_norm)] = {
                'mse_ls'     : res['mse_ls'][0],
                'mse_dft_ncp': res['mse_dft_ncp'][0],
                'mse_dft_pdp': res['mse_dft_pdp'][0],
                'ber_ls'     : res['ber_ls'][0],
                'ber_dft_ncp': res['ber_dft_ncp'][0],
                'ber_dft_pdp': res['ber_dft_pdp'][0],
                'ber_perfect': res['ber_perfect'][0],
            }

    # --- VẼ ĐỒ THỊ 3 ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Thi nghiem 2: Anh huong Doppler len Hieu nang He thong',
                 fontsize=13, fontweight='bold')

    fd_arr = np.array(fd_list)

    for col_idx, snr_db in enumerate(snr_list):
        mse_ls      = np.array([results_doppler[(snr_db, f)]['mse_ls']      for f in fd_list])
        mse_dft_ncp = np.array([results_doppler[(snr_db, f)]['mse_dft_ncp'] for f in fd_list])
        mse_dft_pdp = np.array([results_doppler[(snr_db, f)]['mse_dft_pdp'] for f in fd_list])
        ber_ls      = np.array([results_doppler[(snr_db, f)]['ber_ls']      for f in fd_list])
        ber_dft_ncp = np.array([results_doppler[(snr_db, f)]['ber_dft_ncp'] for f in fd_list])
        ber_dft_pdp = np.array([results_doppler[(snr_db, f)]['ber_dft_pdp'] for f in fd_list])
        ber_perfect = np.array([results_doppler[(snr_db, f)]['ber_perfect'] for f in fd_list])

        # MSE vs fd_norm
        ax = axes[0, col_idx]
        ax.semilogy(fd_arr, mse_ls,      'b-o', label='LS 2D',  markersize=6)
        ax.semilogy(fd_arr, mse_dft_ncp, 'r-s', label=f'DFT (L={N_CP})', markersize=6)
        ax.semilogy(fd_arr, mse_dft_pdp, 'g-^', label=f'DFT (L={window_pdp})', markersize=6)
        ax.set_title(f'MSE vs Doppler (SNR={snr_db} dB)')
        ax.set_xlabel('Normalized Doppler fd_norm')
        ax.set_ylabel('MSE')
        ax.legend(fontsize=9)
        ax.grid(True, which='both', alpha=0.4)

        # BER vs fd_norm
        ax = axes[1, col_idx]
        safe_semilogy(ax, fd_arr, ber_perfect, 'k--',  label='Perfect CSI', linewidth=2)
        safe_semilogy(ax, fd_arr, ber_ls,      'b-o',  label='LS 2D', markersize=6)
        safe_semilogy(ax, fd_arr, ber_dft_ncp, 'r-s',  label=f'DFT (L={N_CP})', markersize=6)
        safe_semilogy(ax, fd_arr, ber_dft_pdp, 'g-^',  label=f'DFT (L={window_pdp})', markersize=6)
        ax.set_title(f'BER vs Doppler (SNR={snr_db} dB)')
        ax.set_xlabel('Normalized Doppler fd_norm')
        ax.set_ylabel('BER')
        ax.legend(fontsize=9)
        ax.grid(True, which='both', alpha=0.4)
        ax.set_ylim([1e-4, 1.0])

    plt.tight_layout()
    plt.savefig('fig3_mobility_impact_mse_ber.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  [OK] fig3_mobility_impact_mse_ber.png da luu.")

    elapsed = time.time() - t0
    print(f"  Thoi gian thi nghiem 2: {elapsed:.1f}s")

    # Xuất CSV
    with open('thesis_results_doppler.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['snr_db', 'fd_norm',
                         'mse_ls', 'mse_dft_ncp', 'mse_dft_pdp',
                         'ber_perfect', 'ber_ls', 'ber_dft_ncp', 'ber_dft_pdp'])
        for snr_db in snr_list:
            for fd_norm in fd_list:
                r = results_doppler[(snr_db, fd_norm)]
                writer.writerow([snr_db, fd_norm,
                                 r['mse_ls'], r['mse_dft_ncp'], r['mse_dft_pdp'],
                                 r['ber_perfect'], r['ber_ls'], r['ber_dft_ncp'], r['ber_dft_pdp']])
    print("  [OK] thesis_results_doppler.csv da luu.")
    return results_doppler


# ===========================================================================
# THÍ NGHIỆM 3: KHẢO SÁT LƯỚI PILOT 2D (Df, Dt) TRADE-OFF
# ===========================================================================

def experiment3_pilot_grid():
    """
    THÍ NGHIỆM 3: Khảo sát các cấu hình lưới pilot (Df, Dt) ở Doppler thấp và cao.
    - Cấu hình khảo sát: (2,1),(4,1),(4,2),(4,4),(8,4).
    - fd_norm = 0.01 (thấp) và fd_norm = 0.05 (cao).
    - Đầu ra: fig4_pilot_spacing_2d_tradeoff.png, fig5_spectral_efficiency_vs_snr.png.
    """
    print("\n" + "=" * 60)
    print("THÍ NGHIỆM 3: KHAO SAT LUOI PILOT 2D (Df, Dt) TRADE-OFF")
    print("=" * 60)
    t0 = time.time()

    configs    = [(2,1), (4,1), (4,2), (4,4), (8,4)]  # (Df, Dt)
    fd_scenarios = [0.01, 0.05]
    snr_ref    = 15.0   # SNR tham chiếu cho biểu đồ tradeoff
    snr_range  = SNR_DB_RANGE
    num_symbols= NUM_SYMBOLS_PER_FRAME
    num_frames = NUM_FRAMES // 2  # Ít frame hơn để tiết kiệm thời gian
    window_pdp = len(CHANNEL_PDP)

    results_pilot = {}   # key: (df, dt, fd_norm)
    se_results    = {}   # key: (df, dt, fd_norm)

    ofdm_qpsk = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING_FREQ,
                            modulation='QPSK', pilot_value=PILOT_VALUE)

    for fd_norm in fd_scenarios:
        for (df, dt) in configs:
            print(f"  (Df={df}, Dt={dt}, fd={fd_norm})...")
            res = run_snr_sweep(fd_norm, df, dt, 'QPSK', num_frames, num_symbols,
                                snr_range, N_CP, window_pdp)
            results_pilot[(df, dt, fd_norm)] = res

            # Tính Pilot Overhead và SE (dùng get_pilot_overhead trả về 3 giá trị)
            ovhd_pct, _, _ = ofdm_qpsk.get_pilot_overhead(df, dt, num_symbols)
            ber_arr = res['ber_dft_pdp']
            se_arr  = np.array([ofdm_qpsk.get_effective_spectral_efficiency(
                                    df, dt, num_symbols, ber) for ber in ber_arr])
            se_results[(df, dt, fd_norm)] = se_arr

    # --- VẼ ĐỒ THỊ 4: BER tại SNR=15dB vs cấu hình pilot ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Thi nghiem 3: Pilot Grid 2D (Df, Dt) Trade-off', fontsize=13, fontweight='bold')

    colors  = ['b', 'r', 'g', 'm', 'orange']
    markers = ['o', 's', '^', 'D', 'v']

    for col, fd_norm in enumerate(fd_scenarios):
        ax = axes[col]
        for ci, (df, dt) in enumerate(configs):
            res = results_pilot[(df, dt, fd_norm)]
            snr_arr = res['snr']
            safe_semilogy(ax, snr_arr, res['ber_dft_pdp'],
                        color=colors[ci], marker=markers[ci],
                        label=f'(Df={df}, Dt={dt})', markersize=5)
        ax.set_title(f'BER vs SNR, fd_norm={fd_norm}')
        ax.set_xlabel('SNR (dB)')
        ax.set_ylabel('BER')
        ax.legend(fontsize=9)
        ax.grid(True, which='both', alpha=0.4)
        ax.set_ylim([1e-4, 1.0])
        ax.set_xlim([snr_range[0], snr_range[-1]])

    plt.tight_layout()
    plt.savefig('fig4_pilot_spacing_2d_tradeoff.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  [OK] fig4_pilot_spacing_2d_tradeoff.png da luu.")

    # --- VẼ ĐỒ THỊ 5: Spectral Efficiency vs SNR ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Thi nghiem 3: Effective Spectral Efficiency vs SNR', fontsize=13, fontweight='bold')

    for col, fd_norm in enumerate(fd_scenarios):
        ax = axes[col]
        for ci, (df, dt) in enumerate(configs):
            se_arr  = se_results[(df, dt, fd_norm)]
            ovhd_pct, _, _ = ofdm_qpsk.get_pilot_overhead(df, dt, num_symbols)
            ax.plot(snr_range, se_arr,
                    color=colors[ci], marker=markers[ci],
                    label=f'(Df={df},Dt={dt}) OH={ovhd_pct:.1f}%', markersize=5)
        ax.set_title(f'Spectral Efficiency, fd_norm={fd_norm}')
        ax.set_xlabel('SNR (dB)')
        ax.set_ylabel('SE (bps/Hz)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.4)
        ax.set_xlim([snr_range[0], snr_range[-1]])
        ax.set_ylim([0, ofdm_qpsk.bits_per_symbol])

    plt.tight_layout()
    plt.savefig('fig5_spectral_efficiency_vs_snr.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  [OK] fig5_spectral_efficiency_vs_snr.png da luu.")

    elapsed = time.time() - t0
    print(f"  Thoi gian thi nghiem 3: {elapsed:.1f}s")
    return results_pilot, se_results


# ===========================================================================
# THÍ NGHIỆM 4: ĐỘ PHỨC TẠP & BẢNG KHUYẾN NGHỊ
# ===========================================================================

def experiment4_complexity_and_recommendations():
    """
    THÍ NGHIỆM 4: Đo thời gian chạy LS vs DFT-based (ms/frame).
    In và xuất Bảng Khuyến nghị Mật độ Pilot theo vùng tốc độ di chuyển.
    Đầu ra: thesis_pilot_recommendations.csv.
    """
    print("\n" + "=" * 60)
    print("THÍ NGHIỆM 4: DO PHUC TAP & BANG KHUYEN NGHI PILOT")
    print("=" * 60)

    num_symbols = NUM_SYMBOLS_PER_FRAME
    df, dt      = PILOT_SPACING_FREQ, PILOT_SPACING_TIME
    window_pdp  = len(CHANNEL_PDP)
    snr_db      = 20.0
    fd_norm     = 0.01
    n_timing    = 50  # Số lần đo để lấy trung bình

    ofdm      = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING_FREQ,
                           modulation='QPSK', pilot_value=PILOT_VALUE)
    estimator = LSEstimator()
    channel   = TimeVaryingRayleighChannel(
        pdp=CHANNEL_PDP, fd_norm=fd_norm, n_fft=N_FFT, n_cp=N_CP,
        num_symbols_per_frame=num_symbols
    )

    print(f"\n  Do thoi gian chay ({n_timing} frames moi phuong phap)...")

    # Sinh frame mẫu để đo thời gian
    max_bits  = N_FFT * num_symbols * ofdm.bits_per_symbol
    data_bits = np.random.randint(0, 2, max_bits)
    X_frame, _, pilot_sub, pilot_sym = ofdm.insert_pilots_2d(data_bits, df, dt, num_symbols)
    X_frame_cp = ofdm.modulate_frame(X_frame)
    Y_frame_cp, _ = channel.transmit_frame(X_frame_cp, snr_db=snr_db)
    Y_frame = ofdm.demodulate_frame(Y_frame_cp)

    # Đo LS 2D
    t_ls = time.perf_counter()
    for _ in range(n_timing):
        estimator.estimate_ls_2d(Y_frame, pilot_sub, pilot_sym, PILOT_VALUE, N_FFT, num_symbols)
    t_ls = (time.perf_counter() - t_ls) / n_timing * 1000  # ms/frame

    # Đo DFT 2D
    t_dft = time.perf_counter()
    for _ in range(n_timing):
        estimator.estimate_dft_2d(Y_frame, pilot_sub, pilot_sym, PILOT_VALUE,
                                   N_FFT, num_symbols, window_len=window_pdp)
    t_dft = (time.perf_counter() - t_dft) / n_timing * 1000  # ms/frame

    print(f"\n  --- Ket qua Do phuc tap (ms/frame) ---")
    print(f"  LS 2D    : {t_ls:.3f} ms/frame")
    print(f"  DFT 2D   : {t_dft:.3f} ms/frame")
    print(f"  Ty le    : DFT overhead = +{t_dft - t_ls:.3f} ms ({(t_dft/t_ls - 1)*100:.1f}%)")

    # --- BẢNG KHUYẾN NGHỊ MẬT ĐỘ PILOT ---
    # Phân loại theo fd_norm -> vùng tốc độ di chuyển:
    # fd_norm = f_d * T_sym, T_sym = (N_FFT + N_CP) / fs
    # Ví dụ: fs=10 MHz, T_sym=80/10M = 8 µs
    # fd_norm=0.01 -> f_d=1250 Hz -> v~km/h (tùy băng tần)

    recommendations = [
        # (Kịch bản, fd_norm điển hình, Df_rec, Dt_rec, Ghi chú)
        ('Static Channel',          '0.000',       4,  7, 'Kenh tinh: co the dung Dt lon'),
        ('Pedestrian (di bo)',       '0.005~0.01',  4,  4, 'Di bo: Dt=4 du tot, Df=4 an toan'),
        ('Urban Vehicular (do thi)', '0.01~0.05',   4,  2, 'Do thi: Dt=2 dam bao theo doi Doppler'),
        ('High-Speed (cao toc)',     '0.05~0.10',   4,  1, 'Cao toc: can Dt=1 (moi sym la pilot)'),
    ]

    print("\n" + "=" * 75)
    print("  BANG KHUYEN NGHI MAT DO PILOT (Df, Dt) THEO VUNG TOC DO DI CHUYEN")
    print("=" * 75)
    print(f"  {'Kich ban':<28} {'fd_norm':<14} {'Df_rec':>6} {'Dt_rec':>6}  {'Ghi chu'}")
    print("-" * 75)
    for (scenario, fd_range, df_r, dt_r, note) in recommendations:
        print(f"  {scenario:<28} {fd_range:<14} {df_r:>6} {dt_r:>6}  {note}")
    print("=" * 75)

    # Xuất CSV
    with open('thesis_pilot_recommendations.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['kich_ban', 'fd_norm_range', 'Df_rec', 'Dt_rec',
                         'pilot_overhead_pct', 'ghi_chu'])
        for (scenario, fd_range, df_r, dt_r, note) in recommendations:
            ovhd, _, _ = ofdm.get_pilot_overhead(df_r, dt_r, num_symbols)
            writer.writerow([scenario, fd_range, df_r, dt_r, f'{ovhd:.1f}', note])

    # Thêm thông tin thời gian chạy vào CSV
    with open('thesis_pilot_recommendations.csv', 'a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([])
        writer.writerow(['Complexity (ms/frame)', '', '', '', '', ''])
        writer.writerow(['LS 2D', f'{t_ls:.3f}', '', '', '', ''])
        writer.writerow(['DFT 2D', f'{t_dft:.3f}', '', '', '', ''])

    print("  [OK] thesis_pilot_recommendations.csv da luu.")
    return t_ls, t_dft


# ===========================================================================
# CHẠY TẤT CẢ THÍ NGHIỆM
# ===========================================================================

if __name__ == '__main__':
    np.random.seed(2025)  # Seed cố định để kết quả có thể tái hiện
    total_start = time.time()

    print("=" * 60)
    print("  BAT DAU CHAY TOAN BO THI NGHIEM DO AN")
    print(f"  Cau hinh: N_FFT={N_FFT}, N_CP={N_CP}, NUM_FRAMES={NUM_FRAMES}")
    print(f"  SNR_DB_RANGE: {SNR_DB_RANGE[0]} - {SNR_DB_RANGE[-1]} dB (buoc {SNR_DB_RANGE[1]-SNR_DB_RANGE[0]} dB)")
    print("=" * 60)

    # Bước 1-4: Chạy các thí nghiệm
    experiment1_baseline()
    experiment2_mobility()
    experiment3_pilot_grid()
    experiment4_complexity_and_recommendations()

    total_elapsed = time.time() - total_start
    print("\n" + "=" * 60)
    print(f"  HOAN THANH TOAN BO THI NGHIEM!")
    print(f"  Tong thoi gian: {total_elapsed:.1f}s ({total_elapsed/60:.1f} phut)")
    print("  Cac file dau ra:")
    print("    - fig1_mse_vs_snr.png")
    print("    - fig2_ber_vs_snr.png")
    print("    - fig3_mobility_impact_mse_ber.png")
    print("    - fig4_pilot_spacing_2d_tradeoff.png")
    print("    - fig5_spectral_efficiency_vs_snr.png")
    print("    - thesis_results_snr.csv")
    print("    - thesis_results_doppler.csv")
    print("    - thesis_pilot_recommendations.csv")
    print("=" * 60)
