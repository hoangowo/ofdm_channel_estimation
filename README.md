# Mô Phỏng Hệ Thống OFDM và Ước Lượng Kênh Truyền Least Squares (LS) & MMSE dựa trên Pilot

Dự án mô phỏng hệ thống **OFDM (Orthogonal Frequency Division Multiplexing)** và thuật toán **Ước lượng Kênh truyền LS (Least Squares)** và **Simplified MMSE** với cấu trúc Pilot dạng Comb-type trong kênh fading Rayleigh đa đường theo chuẩn giáo trình *MIMO-OFDM Wireless Communications with MATLAB* (Yong Soo Cho et al.).

---

## 📌 Hướng dẫn Nộp bài & Yêu cầu của Giảng viên

Dự án đáp ứng đầy đủ các yêu cầu nộp bài của giảng viên:
- **Repo GitHub / Fork**: Tổ chức mã nguồn chuẩn mô đun hoá.
- **Hướng dẫn README**: Chi tiết cách cài đặt, cách chạy, nguồn dữ liệu và sinh lại kết quả.
- **Script chạy tự động**: Tích hợp `run_all_experiments.py` tự động chạy lại toàn bộ thí nghiệm và sinh 5 hình đồ thị + bảng kết quả CSV.
- **Báo cáo & Số liệu**: Đầy đủ hình ảnh đồ thị PNG, biểu đồ chòm sao IQ, bảng CSV kết quả và báo cáo kỹ thuật.

---

## 🎲 Nguồn Dữ Liệu (Data Sources & Generation)

Dự án sử dụng **Dữ liệu tổng hợp ngẫu nhiên (Monte Carlo Simulation)** — không yêu cầu bất kỳ tập dữ liệu bên ngoài nào.

- **Bit dữ liệu đầu vào**: Chuỗi ngẫu nhiên phân phối Bernoulli(0.5) được sinh tự động qua `np.random.randint`.
- **Kênh truyền Rayleigh**: Đáp ứng xung kênh $h_i \sim \mathcal{CN}(0, P_i)$ với Power Delay Profile (PDP) 3-tap chuẩn hóa `[0.7, 0.2, 0.1]`.
- **Nhiễu AWGN**: Nhiễu Phức Gauss $\mathcal{CN}(0, \sigma_z^2)$ với công suất nhiễu được tính trực tiếp từ dải SNR (0 dB đến 24 dB).

---

## 🛠️ Yêu Cầu Cài Đặt (Installation)

### 1. Yêu cầu môi trường
- Python 3.8+ (đã kiểm thử trên Python 3.14).

### 2. Cài đặt thư viện phụ thuộc
Chạy lệnh sau trong terminal:
```bash
pip install -r requirements.txt
```

Nội dung `requirements.txt`:
```
numpy>=1.21.0
matplotlib>=3.4.0
```

---

## 🚀 Cách Chạy & Tái Tạo Kết Quả (Reproduction Guide)

### 1. Chạy toàn bộ thí nghiệm và sinh hình/bảng kết quả (Khuyên dùng)
Để chạy toàn bộ mô phỏng Monte Carlo, xuất bảng dữ liệu CSV và tự động vẽ tất cả 5 đồ thị:
```bash
python run_all_experiments.py
```

Khi chạy xong, hệ thống sẽ tự động tạo các tệp kết quả:
- `simulation_results.csv`: Bảng tổng hợp số liệu BER và MSE chi tiết cho Perfect CSI, LS và MMSE.
- `ber_vs_snr.png`: Đồ thị so sánh BER của Perfect CSI vs LS Estimation vs MMSE Estimation.
- `mse_vs_snr.png`: Đồ thị MSE của LS Estimation vs MMSE Estimation theo SNR.
- `pilot_spacing_comparison.png`: Đồ thị phân tích đánh đổi mật độ Pilot ($D_f = 2, 4, 8$).
- `constellation_comparison.png`: Biểu đồ chòm sao IQ trước và sau cân bằng 1-tap ZF.
- `channel_estimation_fit.png`: Snapshot độ khớp đáp ứng tần số kênh $|H(k)|$ thực tế vs LS vs MMSE.

### 2. Chạy kiểm thử đơn vị (Unit Tests)
Để kiểm tra tính vuông góc OFDM và độ chính xác của IFFT/FFT trong điều kiện kênh lý tưởng (Noiseless Loopback):
```bash
python test_modules.py
```

---

## 📂 Cấu Trúc Thư Mục Dự Án (Project Structure)

```
ofdm_channel_estimation/
├── config.py                 # Hằng số cấu hình hệ thống (N_FFT, N_CP, D_f, SNR, PDP)
├── channel.py                # Module sinh kênh Rayleigh đa đường & nhiễu AWGN
├── ofdm_system.py            # Module phát/thu OFDM, điều chế QPSK/16QAM, IFFT/FFT, CP
├── estimator.py              # Module ước lượng kênh LS, MMSE, nội suy tuyến tính & cân bằng 1-tap
├── main.py                   # Script chạy mô phỏng chính
├── run_all_experiments.py    # ⭐ Script chạy toàn bộ thí nghiệm, xuất CSV và 5 đồ thị
├── test_modules.py           # Unit tests kiểm tra tính chính xác của các khối
├── requirements.txt          # Khai báo các thư viện phụ thuộc
├── simulation_results.csv    # Bảng số liệu kết quả mô phỏng xuất tự động
├── ber_vs_snr.png            # Biểu đồ BER vs SNR (Perfect vs LS vs MMSE)
├── mse_vs_snr.png            # Biểu đồ MSE vs SNR (LS vs MMSE)
├── pilot_spacing_comparison.png # Biểu đồ so sánh khoảng cách pilot (D_f = 2, 4, 8)
├── constellation_comparison.png # Biểu đồ chòm sao IQ trước và sau cân bằng ZF
└── channel_estimation_fit.png   # Snapshot độ khớp đáp ứng tần số kênh |H(k)|
```

---

## 📊 Kết Quả Mô Phỏng & Ý Nghĩa Kỹ Thuật

### 1. Bảng số liệu mô phỏng đại diện (trích từ `simulation_results.csv`)

| SNR (dB) | BER Perfect CSI | BER LS ($D_f=4$) | BER MMSE ($D_f=4$) | MSE LS ($D_f=4$) | MSE MMSE ($D_f=4$) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 0 | 0.188062 | 0.261224 | 0.261224 | 6.890e-01 | 4.237e-01 |
| 4 | 0.096365 | 0.146292 | 0.146292 | 2.700e-01 | 2.212e-01 |
| 8 | 0.039776 | 0.063708 | 0.063708 | 1.106e-01 | 1.023e-01 |
| 12 | 0.015396 | 0.025063 | 0.025063 | 4.270e-02 | 4.202e-02 |
| 16 | 0.005682 | 0.009422 | 0.009422 | 1.748e-02 | 1.755e-02 |
| 20 | 0.002271 | 0.003875 | 0.003875 | 7.074e-03 | 7.168e-03 |
| 24 | 0.000870 | 0.001479 | 0.001479 | 3.107e-03 | 3.152e-03 |

### 2. Phân tích Ý nghĩa Kỹ thuật
- **So sánh MSE: LS vs MMSE**: MMSE cải thiện MSE kênh vượt trội hơn LS ở vùng SNR thấp (0–10 dB) nhờ hệ số co Wiener $SNR/(SNR+1)$ làm giảm tác động của nhiễu AWGN trên pilot. Ở vùng SNR cao, MMSE tiệm cận LS.
- **So sánh BER: LS vs MMSE (QPSK)**: Đối với điều chế QPSK với bộ cân bằng 1-tap, BER của MMSE trùng khớp hoàn toàn với LS. Nguyên nhân: $\hat{H}_{MMSE} = c \cdot \hat{H}_{LS}$ ($c > 0$) nên pha tín hiệu sau cân bằng $Y/\hat{H}$ được bảo toàn 100%, mà QPSK quyết định bit chỉ dựa trên dấu phần Thực/Ảo (góc phần tư).
- **Suy giảm hiệu năng 2–3 dB**: Thuật toán LS Estimation bị suy giảm khoảng 2 dB đến 3 dB so với đường cơ sở Perfect CSI do ảnh hưởng của nhiễu AWGN tác động lên các vị trí pilot và lỗi nội suy tuyến tính giữa các pilot subcarrier.
- **Đặc tính MSE**: Đường cong MSE giảm tuyến tính trên thang log-log khi SNR tăng.
- **Đánh đổi Mật độ Pilot ($D_f$)**:
  - $D_f = 2$ (Pilot dày): Đạt BER và MSE tiệm cận Perfect CSI nhưng tổn hao 50% dung lượng dành cho pilot overhead.
  - $D_f = 4$ (Pilot chuẩn): Điểm cân bằng hợp lý giữa độ chính xác ước lượng và tốc độ truyền dữ liệu (25% overhead).
  - $D_f = 8$ (Pilot thưa): Tối ưu tốc độ dữ liệu (chỉ 12.5% overhead) nhưng lỗi nội suy tăng cao ở vùng SNR cao.
