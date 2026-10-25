"""
config.py - Tập trung toàn bộ tham số cấu hình cho hệ thống ước lượng kênh OFDM.
Đồ án: "Pilot-Based Channel Estimation for OFDM over Time-Varying Fading Channels"
Sinh viên: Đào Kim Hoàng - MSSV: 2211068
"""
import numpy as np

# ===========================================================================
# I. THAM SỐ HỆ THỐNG OFDM CƠ BẢN (GIỮ NGUYÊN ĐỂ LÀM BASELINE)
# ===========================================================================
N_FFT = 64               # Tổng số sóng mang phụ (FFT size), tiêu chuẩn thường dùng 64, 128, 256...
N_CP  = 16               # Độ dài Cyclic Prefix (số mẫu), N_CP >= L_channel - 1 để tránh ISI

MODULATION  = 'QPSK'     # Sơ đồ điều chế mặc định ('QPSK' hoặc '16QAM')
PILOT_VALUE = 1.0 + 0.0j # Giá trị ký tự pilot chuẩn (BPSK: biên độ 1, pha 0)

# Khoảng cách sóng mang pilot mặc định (Df) - dùng cho chế độ 1D (baseline cũ)
PILOT_SPACING = 4

# Dải SNR mô phỏng: từ 0 dB đến 24 dB, bước 2 dB (tổng 13 điểm)
SNR_DB_RANGE = np.arange(0, 26, 2)

# Số ký tự OFDM cũ (giữ để tương thích ngược với baseline)
NUM_OFDM_SYMBOLS = 2000

# Hồ sơ phân bố năng lượng trễ (Power Delay Profile - PDP) chuẩn hóa
# Kênh 3-tap: tap 0 mạnh nhất (0.7), tap 1 (0.2), tap 2 (0.1)
CHANNEL_PDP = np.array([0.7, 0.2, 0.1])

# ===========================================================================
# II. THAM SỐ KHUNG TRUYỀN (FRAME) - MỚI THÊM CHO ĐỒ ÁN
# ===========================================================================
# Số ký tự OFDM trong một frame (tương tự cấu trúc 1 slot LTE/NR: 14 symbol/slot)
NUM_SYMBOLS_PER_FRAME = 14

# Số frame mô phỏng Monte Carlo tại mỗi mức SNR
# (200 frame x 14 symbol = 2800 symbol/SNR, đủ tốt cho ước lượng BER tin cậy)
NUM_FRAMES = 200

# ===========================================================================
# III. THAM SỐ LƯỚI PILOT 2D (FREQUENCY x TIME) - MỚI THÊM
# ===========================================================================
# --- Chiều tần số (Frequency axis) ---
# Khoảng cách sóng mang pilot mặc định Df = 4 (tức cứ 4 subcarrier đặt 1 pilot)
# Điều kiện Nyquist tần số: Df <= N_FFT / (2 * (L_channel - 1)) = 64/(2*2) = 16 -> Df=4 rất an toàn
PILOT_SPACING_FREQ = 4

# Danh sách Df để khảo sát trade-off: dày (Df=2) - chuẩn (Df=4) - thưa (Df=8)
DF_LIST = [2, 4, 8]

# --- Chiều thời gian (Time axis) ---
# Khoảng cách symbol pilot mặc định Dt = 2 (cứ 2 OFDM symbol đặt 1 hàng pilot)
# Điều kiện Nyquist thời gian: Dt <= 1/(2*fd_norm_max) ~ 1/(2*0.1)=5 -> Dt=2 an toàn với fd=0.05
PILOT_SPACING_TIME = 2

# Danh sách Dt để khảo sát: liên tục (Dt=1) - chuẩn (Dt=2) - trung bình (Dt=4) - thưa (Dt=7)
DT_LIST = [1, 2, 4, 7]

# ===========================================================================
# IV. THAM SỐ DOPPLER - ĐẶC TRƯNG MỨC ĐỘ DI ĐỘNG
# ===========================================================================
# Tần số Doppler chuẩn hóa: fd_norm = f_d * T_sym
# Trong đó: f_d = v * f_c / c (Hz), T_sym = (N_FFT + N_CP) / f_sampling
# Phân loại theo kịch bản di động:
#   fd_norm = 0.000 -> Kênh tĩnh (Static channel, kiểm tra baseline)
#   fd_norm = 0.010 -> Đi bộ nhanh (~5 km/h ở 2.4 GHz)
#   fd_norm = 0.030 -> Đô thị xe máy chậm (~15 km/h)
#   fd_norm = 0.050 -> Đô thị xe máy nhanh (~25 km/h)
#   fd_norm = 0.100 -> Cao tốc (~50+ km/h, Doppler mạnh)
FD_NORM_LIST = [0.0, 0.01, 0.03, 0.05, 0.10]
