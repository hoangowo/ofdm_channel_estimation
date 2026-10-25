"""
ofdm_system.py - Module xử lý điều chế/giải điều chế OFDM cấp độ symbol và frame.
Bao gồm ánh xạ bit QPSK/16-QAM, chèn pilot, IFFT/FFT, và xử lý Pilot Grid 2D.
"""
import numpy as np


class OFDMSystem:
    """
    Module quản lý toàn bộ quá trình Điều chế, Giải điều chế OFDM, 
    Ánh xạ bit (Mapping/Demapping QPSK & 16-QAM), và Chèn tín hiệu Pilot.
    Hỗ trợ cả xử lý theo từng symbol (1D - baseline) và theo Frame 2D (mới).
    """
    def __init__(self, n_fft=64, n_cp=16, pilot_spacing=4, modulation='QPSK', pilot_value=1.0+0.0j):
        """
        Khởi tạo thông số cho khối OFDM System.
        - n_fft: Tổng số sóng mang phụ (subcarriers), mặc định 64.
        - n_cp: Độ dài Cyclic Prefix (số mẫu), mặc định 16.
        - pilot_spacing: Khoảng cách giữa các subcarrier pilot (Df), mặc định 4.
        - modulation: Sơ đồ điều chế ('QPSK' hoặc '16QAM').
        - pilot_value: Giá trị ký tự pilot biết trước (BPSK: 1.0 + 0.0j).
        """
        self.n_fft         = n_fft
        self.n_cp          = n_cp
        self.pilot_spacing = pilot_spacing
        self.modulation    = modulation.upper()
        self.pilot_value   = pilot_value

        # --- Định nghĩa vị trí các subcarrier ---
        # pilot_indices: Mảng chứa vị trí các subcarrier pilot (ví dụ: [0, 4, 8, 12, ..., 60])
        self.pilot_indices = np.arange(0, self.n_fft, self.pilot_spacing)
        
        # data_indices: Mảng chứa các vị trí subcarrier còn lại dành cho dữ liệu (Data)
        self.data_indices = np.array([i for i in range(self.n_fft) if i not in self.pilot_indices])
        
        self.num_pilots = len(self.pilot_indices) # Số lượng pilot trong 1 OFDM symbol
        self.num_data   = len(self.data_indices)  # Số lượng subcarrier data trong 1 OFDM symbol

        # Khởi tạo bảng chòm sao điều chế (Constellation Map)
        self._init_constellation()

    def _init_constellation(self):
        """
        Khởi tạo từ điển ánh xạ chòm sao Gray-code và chuẩn hóa công suất trung bình Es = 1.
        """
        if self.modulation == 'QPSK':
            self.bits_per_symbol = 2 # QPSK mã hóa 2 bit trên 1 symbol
            
            # Hệ số chuẩn hóa công suất cho QPSK (Es = 1) -> chia cho sqrt(2)
            scale = 1.0 / np.sqrt(2.0)
            
            # Bảng ánh xạ Gray code QPSK: Key là tuple bit (b0, b1), Value là số phức (I + jQ)
            self.constellation = {
                (0, 0): (1.0 + 1.0j) * scale,   # Góc phần tư thứ I
                (0, 1): (-1.0 + 1.0j) * scale,  # Góc phần tư thứ II
                (1, 1): (-1.0 - 1.0j) * scale,  # Góc phần tư thứ III
                (1, 0): (1.0 - 1.0j) * scale    # Góc phần tư thứ IV
            }
        elif self.modulation == '16QAM':
            self.bits_per_symbol = 4 # 16-QAM mã hóa 4 bit trên 1 symbol
            
            # Hệ số chuẩn hóa công suất cho 16-QAM (Es = 1) -> chia cho sqrt(10)
            scale = 1.0 / np.sqrt(10.0)
            
            # Ánh xạ 1D Gray-code cho các biên độ [-3, -1, +1, +3]
            mapping_1d = {
                (0, 0): -3.0,
                (0, 1): -1.0,
                (1, 1): 1.0,
                (1, 0): 3.0
            }
            # Kết hợp trục Thực (I) và Trục Ảo (Q) để tạo chòm sao 16-QAM 2D
            self.constellation = {}
            for b_re, val_re in mapping_1d.items():
                for b_im, val_im in mapping_1d.items():
                    bit_tuple = b_re + b_im # Ghép thành tuple 4 bit (b0, b1, b2, b3)
                    self.constellation[bit_tuple] = (val_re + 1j * val_im) * scale
        else:
            raise ValueError(f"Sơ đồ điều chế '{self.modulation}' không hỗ trợ! Dùng 'QPSK' hoặc '16QAM'.")

        # Chuyển bảng chòm sao thành danh sách cặp (bit_tuple, constellation_point) để phục vụ giải điều chế nhanh
        self.constellation_list = list(self.constellation.items())

    def map_bits(self, bits):
        """
        CHỨC NĂNG PHÍA PHÁT (TX): Ánh xạ mảng bit nhị phân (0, 1) sang các symbol phức.
        - Input: bits (mảng 1D các số 0 và 1)
        - Output: Mảng m các số phức symbol điều chế (QPSK/16QAM)
        """
        # Tính số lượng symbol thu được = Tổng số bit / số bit mỗi symbol
        num_symbols = len(bits) // self.bits_per_symbol
        symbols = np.zeros(num_symbols, dtype=complex)
        
        for i in range(num_symbols):
            # Cắt lấy nhóm bit tương ứng cho 1 symbol (2 bit với QPSK, 4 bit với 16QAM)
            bit_chunk = tuple(bits[i * self.bits_per_symbol : (i + 1) * self.bits_per_symbol])
            
            # Tra bảng chòm sao để lấy giá trị số phức I + jQ tương ứng
            symbols[i] = self.constellation[bit_chunk]
            
        return symbols

    def demap_symbols(self, symbols):
        """
        CHỨC NĂNG PHÍA THU (RX): Giải điều chế theo khoảng cách Euclidean tối thiểu (Minimum Euclidean Distance).
        - So sánh từng symbol nhận được với tất cả các điểm chòm sao lý tưởng để tìm điểm gần nhất.
        - Input: mảng số phức symbols nhận được sau khi cân bằng kênh.
        - Output: mảng 1D các bit nhị phân khôi phục (0 và 1).
        """
        recovered_bits = []
        for sym in symbols:
            min_dist = float('inf')
            best_bits = None
            
            # Duyệt qua tất cả các điểm chòm sao chuẩn trong bộ nhớ
            for bit_tuple, const_point in self.constellation_list:
                # Tính khoảng cách Euclidean bình phương: |y - c|^2
                dist = np.abs(sym - const_point)**2
                if dist < min_dist:
                    min_dist = dist
                    best_bits = bit_tuple
                    
            # Thêm nhóm bit tìm được vào luồng bit khôi phục
            recovered_bits.extend(best_bits)
            
        return np.array(recovered_bits, dtype=int)

    def insert_pilots(self, data_symbols):
        """
        CHỨC NĂNG PHÍA PHÁT (TX): Ghép ký tự Pilot và ký tự Data vào khung tần số N_FFT subcarrier (Comb-type).
        - Gán giá trị pilot_value (1.0+0.0j) vào các vị trí pilot_indices.
        - Gán các ký tự dữ liệu (data_symbols) vào các vị trí data_indices.
        - Output: Vector miền tần số X[k] có độ dài N_FFT (64).
        """
        X = np.zeros(self.n_fft, dtype=complex)
        X[self.pilot_indices] = self.pilot_value
        X[self.data_indices]  = data_symbols
        return X

    def modulate(self, X_freq):
        """
        CHỨC NĂNG PHÍA PHÁT (TX): Biến đổi tín hiệu miền tần số sang miền thời gian (IFFT) & Chèn Cyclic Prefix (CP).
        - Bước 1: Thực hiện IFFT chuẩn hóa năng lượng: x[n] = 1/sqrt(N_FFT) * IFFT(X[k])
        - Bước 2: Chèn CP bằng cách lấy N_CP mẫu cuối dán lên đầu chuỗi tín hiệu thời gian.
        - Output: Mảng tín hiệu thời gian x_cp có độ dài N_FFT + N_CP (64 + 16 = 80 mẫu).
        """
        # 1. Normalized IFFT: Biến đổi từ miền tần số (subcarrier) sang miền thời gian
        # np.fft.ifft mặc định chia cho N_FFT, nên nhân sqrt(N_FFT) để chuẩn hóa năng lượng
        x_time = np.fft.ifft(X_freq) * np.sqrt(self.n_fft)
        
        # 2. Thêm Cyclic Prefix (CP): Lấy N_CP mẫu cuối cùng dán vào trước tín hiệu
        cp    = x_time[-self.n_cp:]
        x_cp  = np.concatenate([cp, x_time])
        return x_cp

    def demodulate(self, y_cp):
        """
        CHỨC NĂNG PHÍA THU (RX): Bỏ Cyclic Prefix (CP) & Biến đổi tín hiệu miền thời gian sang miền tần số (FFT).
        - Bước 1: Cắt bỏ N_CP mẫu đầu tiên (Cyclic Prefix).
        - Bước 2: Thực hiện FFT chuẩn hóa năng lượng: Y[k] = 1/sqrt(N_FFT) * FFT(y[n])
        - Output: Vector miền tần số Y[k] nhận được có độ dài N_FFT (64).
        """
        # 1. Bỏ Cyclic Prefix (CP): Chỉ giữ lại N_FFT mẫu tín hiệu phía sau
        y_time  = y_cp[self.n_cp : self.n_cp + self.n_fft]
        
        # 2. Normalized FFT: Biến đổi từ miền thời gian trở lại miền tần số (subcarriers)
        Y_freq  = np.fft.fft(y_time) / np.sqrt(self.n_fft)
        return Y_freq

    # ===========================================================================
    # PHẦN MỚI: XỬ LÝ CẤP ĐỘ FRAME 2D (FREQUENCY x TIME)
    # ===========================================================================

    def get_pilot_grid(self, df, dt, num_symbols):
        """
        Tính toán và trả về vị trí lưới Pilot 2D trong frame (df x dt).

        Tham số:
        - df          : Khoảng cách tần số giữa các pilot (Pilot Spacing Frequency).
                        Ví dụ df=4: pilot ở subcarrier 0, 4, 8, ..., 60.
        - dt          : Khoảng cách thời gian giữa các hàng pilot (Pilot Spacing Time).
                        Ví dụ dt=2: pilot ở symbol 0, 2, 4, ..., tức là cứ 2 symbol đặt 1 hàng pilot.
        - num_symbols : Số ký tự OFDM trong frame.

        Trả về:
        - pilot_subcarrier_indices : Mảng 1D chứa các chỉ số subcarrier pilot theo chiều tần số.
        - pilot_symbol_indices     : Mảng 1D chứa các chỉ số symbol (cột thời gian) chứa pilot.

        Lưu ý về biên thời gian:
        - Symbol 0 luôn là pilot (đặt ở đầu frame để khởi tạo).
        - Nếu symbol cuối cùng (num_symbols - 1) không nằm trong danh sách pilot,
          ta thêm nó vào để tránh ngoại suy ra ngoài biên (extrapolation gây sai số lớn).
        """
        # Vị trí subcarrier pilot theo chiều tần số (comb-type)
        pilot_subcarrier_indices = np.arange(0, self.n_fft, df)

        # Vị trí symbol pilot theo chiều thời gian: bắt đầu từ 0, bước dt
        pilot_symbol_indices = np.arange(0, num_symbols, dt)

        # Đảm bảo symbol cuối frame cũng là pilot để tránh ngoại suy về phía phải
        last_sym = num_symbols - 1
        if pilot_symbol_indices[-1] != last_sym:
            pilot_symbol_indices = np.append(pilot_symbol_indices, last_sym)

        return pilot_subcarrier_indices, pilot_symbol_indices

    def insert_pilots_2d(self, data_bits, df, dt, num_symbols):
        """
        CHỨC NĂNG PHÍA PHÁT (TX): Tạo Frame tần số 2D kích thước (N_FFT, num_symbols)
        bằng cách chèn pilot theo lưới 2D (df, dt) và điền data vào các vị trí còn lại.

        Tham số:
        - data_bits   : Mảng bit dữ liệu 1D (sẽ được tự động cắt vừa đủ).
        - df          : Khoảng cách pilot theo tần số.
        - dt          : Khoảng cách pilot theo thời gian.
        - num_symbols : Số OFDM symbol trong frame.

        Trả về:
        - X_frame     : Ma trận miền tần số frame 2D, kích thước (N_FFT, num_symbols).
                        Dtype complex. Chứa pilot và data.
        - tx_bits_used: Mảng bit thực sự được dùng (đã cắt vừa đủ để ghép vào frame).
        - pilot_sub   : Vị trí subcarrier pilot (để truyền cho estimator).
        - pilot_sym   : Vị trí symbol pilot (để truyền cho estimator).
        """
        # Lấy vị trí lưới pilot 2D
        pilot_sub, pilot_sym = self.get_pilot_grid(df, dt, num_symbols)

        # Xác định tập hợp các ô (k, n) là pilot trong frame
        pilot_set = set()
        for n in pilot_sym:
            for k in pilot_sub:
                pilot_set.add((k, n))

        # Số subcarrier data = tổng subcarrier - số pilot trong TỪNG SYMBOL
        # Pilot symbol: chứa cả pilot và data (data ở các subcarrier không phải pilot)
        # Non-pilot symbol: toàn bộ N_FFT subcarrier đều là data
        num_data_cells = self.n_fft * num_symbols - len(pilot_set)

        # Tính số bit cần thiết và ánh xạ sang symbol
        num_data_symbols = num_data_cells
        num_bits_needed  = num_data_symbols * self.bits_per_symbol

        # Cắt bit đầu vào cho vừa đủ
        tx_bits_used = data_bits[:num_bits_needed]
        data_symbols = self.map_bits(tx_bits_used)

        # Điền vào frame 2D
        X_frame     = np.zeros((self.n_fft, num_symbols), dtype=complex)
        data_idx    = 0

        for n in range(num_symbols):
            for k in range(self.n_fft):
                if (k, n) in pilot_set:
                    # Ô pilot: gán giá trị pilot chuẩn
                    X_frame[k, n] = self.pilot_value
                else:
                    # Ô data: lấy symbol tiếp theo từ luồng data
                    X_frame[k, n] = data_symbols[data_idx]
                    data_idx += 1

        return X_frame, tx_bits_used, pilot_sub, pilot_sym

    def modulate_frame(self, X_frame):
        """
        CHỨC NĂNG PHÍA PHÁT (TX): Biến đổi toàn bộ frame tần số 2D sang frame thời gian.
        Áp dụng IFFT + CP cho từng OFDM symbol trong frame.

        Tham số:
        - X_frame : Ma trận tần số 2D kích thước (N_FFT, num_symbols).

        Trả về:
        - X_frame_cp : Ma trận thời gian 2D kích thước (num_symbols, symbol_len).
                       Mỗi hàng là 1 OFDM symbol đã có CP, sẵn sàng phát.
        """
        num_symbols  = X_frame.shape[1]
        symbol_len   = self.n_fft + self.n_cp
        X_frame_cp   = np.zeros((num_symbols, symbol_len), dtype=complex)

        for n in range(num_symbols):
            # Lấy vector tần số của symbol n và điều chế (IFFT + CP)
            X_freq_n       = X_frame[:, n]
            X_frame_cp[n]  = self.modulate(X_freq_n)

        return X_frame_cp

    def demodulate_frame(self, Y_frame_cp):
        """
        CHỨC NĂNG PHÍA THU (RX): Biến đổi frame thời gian 2D về frame tần số 2D.
        Áp dụng bỏ CP + FFT cho từng OFDM symbol trong frame.

        Tham số:
        - Y_frame_cp : Ma trận thời gian 2D kích thước (num_symbols, symbol_len).

        Trả về:
        - Y_frame : Ma trận tần số 2D kích thước (N_FFT, num_symbols).
        """
        num_symbols = Y_frame_cp.shape[0]
        Y_frame     = np.zeros((self.n_fft, num_symbols), dtype=complex)

        for n in range(num_symbols):
            # Giải điều chế từng symbol: bỏ CP + FFT
            Y_frame[:, n] = self.demodulate(Y_frame_cp[n])

        return Y_frame

    def equalize_and_decode_frame(self, Y_frame, H_est_frame, pilot_sub, pilot_sym):
        """
        CHỨC NĂNG PHÍA THU (RX): Cân bằng kênh ZF và giải điều chế toàn bộ frame 2D.

        Tham số:
        - Y_frame       : Ma trận tần số thu 2D kích thước (N_FFT, num_symbols).
        - H_est_frame   : Ma trận ước lượng kênh 2D kích thước (N_FFT, num_symbols).
        - pilot_sub     : Vị trí subcarrier pilot (để loại bỏ khỏi luồng dữ liệu).
        - pilot_sym     : Vị trí symbol pilot (để loại bỏ khỏi luồng dữ liệu).

        Trả về:
        - rx_bits : Mảng bit giải điều chế 1D.
        """
        num_symbols = Y_frame.shape[1]
        pilot_set   = set()
        for n in pilot_sym:
            for k in pilot_sub:
                pilot_set.add((k, n))

        # Cân bằng ZF: X_hat[k,n] = Y[k,n] / H_est[k,n]
        eps   = 1e-12
        H_safe = np.where(np.abs(H_est_frame) < eps, eps, H_est_frame)
        X_hat  = Y_frame / H_safe

        # Trích xuất các symbol data theo đúng thứ tự chèn vào frame
        data_symbols_rx = []
        for n in range(num_symbols):
            for k in range(self.n_fft):
                if (k, n) not in pilot_set:
                    data_symbols_rx.append(X_hat[k, n])

        data_symbols_rx = np.array(data_symbols_rx)
        rx_bits = self.demap_symbols(data_symbols_rx)
        return rx_bits

    def get_pilot_overhead(self, df, dt, num_symbols):
        """
        Tính tỷ lệ phần trăm tài nguyên trong frame dùng làm pilot (Pilot Overhead).

        Công thức:
            Pilot Overhead (%) = (Số ô pilot) / (Tổng số ô trong frame) * 100

        Tham số:
        - df          : Khoảng cách tần số pilot.
        - dt          : Khoảng cách thời gian pilot.
        - num_symbols : Số symbol trong frame.

        Trả về:
        - overhead_percent : Tỷ lệ pilot overhead (%).
        - num_pilot_cells  : Số ô pilot trong frame.
        - total_cells      : Tổng số ô tài nguyên trong frame.
        """
        pilot_sub, pilot_sym = self.get_pilot_grid(df, dt, num_symbols)
        num_pilot_cells   = len(pilot_sub) * len(pilot_sym)
        total_cells       = self.n_fft * num_symbols
        overhead_percent  = num_pilot_cells / total_cells * 100.0
        return overhead_percent, num_pilot_cells, total_cells

    def get_effective_spectral_efficiency(self, df, dt, num_symbols, ber):
        """
        Tính Hiệu suất Phổ Hiệu dụng (Effective Spectral Efficiency - SE) tính bằng bps/Hz.

        Công thức:
            SE = bits_per_symbol * (1 - Pilot_Overhead) * (N_FFT / (N_FFT + N_CP)) * (1 - BER)

        Trong đó:
        - bits_per_symbol          : Số bit/symbol của điều chế (2 với QPSK, 4 với 16QAM).
        - (1 - Pilot_Overhead)     : Hiệu suất sau khi trừ phần dùng làm pilot.
        - N_FFT / (N_FFT + N_CP)   : Hiệu suất sau khi trừ Cyclic Prefix overhead.
        - (1 - BER)                : Tỷ lệ bit được truyền chính xác.

        Tham số:
        - df, dt      : Cấu hình lưới pilot.
        - num_symbols : Số symbol trong frame.
        - ber         : Tỷ lệ bit lỗi (Bit Error Rate).

        Trả về:
        - se : Hiệu suất phổ hiệu dụng (bps/Hz).
        """
        overhead_pct, _, _ = self.get_pilot_overhead(df, dt, num_symbols)
        # Hiệu suất thực sau pilot overhead (tỷ lệ data)
        data_efficiency  = 1.0 - overhead_pct / 100.0
        # Hiệu suất sau CP overhead
        cp_efficiency    = self.n_fft / (self.n_fft + self.n_cp)
        # Tính SE hiệu dụng
        se = self.bits_per_symbol * data_efficiency * cp_efficiency * (1.0 - ber)
        return se
