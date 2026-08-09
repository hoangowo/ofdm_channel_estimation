import numpy as np

class OFDMSystem:
    """
    Module quản lý toàn bộ quá trình Điều chế, Giải điều chế OFDM, 
    Ánh xạ bit (Mapping/Demapping QPSK & 16-QAM), và Chèn tín hiệu Pilot.
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
        self.n_fft = n_fft
        self.n_cp = n_cp
        self.pilot_spacing = pilot_spacing
        self.modulation = modulation.upper()
        self.pilot_value = pilot_value

        # --- Định nghĩa vị trí các subcarrier ---
        # pilot_indices: Mảng chứa vị trí các subcarrier pilot (ví dụ: [0, 4, 8, 12, ..., 60])
        self.pilot_indices = np.arange(0, self.n_fft, self.pilot_spacing)
        
        # data_indices: Mảng chứa các vị trí subcarrier còn lại dành cho dữ liệu (Data)
        self.data_indices = np.array([i for i in range(self.n_fft) if i not in self.pilot_indices])
        
        self.num_pilots = len(self.pilot_indices) # Số lượng pilot trong 1 OFDM symbol
        self.num_data = len(self.data_indices)     # Số lượng subcarrier data trong 1 OFDM symbol

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
        X[self.data_indices] = data_symbols
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
        cp = x_time[-self.n_cp:]
        x_cp = np.concatenate([cp, x_time])
        return x_cp

    def demodulate(self, y_cp):
        """
        CHỨC NĂNG PHÍA THU (RX): Bỏ Cyclic Prefix (CP) & Biến đổi tín hiệu miền thời gian sang miền tần số (FFT).
        - Bước 1: Cắt bỏ N_CP mẫu đầu tiên (Cyclic Prefix).
        - Bước 2: Thực hiện FFT chuẩn hóa năng lượng: Y[k] = 1/sqrt(N_FFT) * FFT(y[n])
        - Output: Vector miền tần số Y[k] nhận được có độ dài N_FFT (64).
        """
        # 1. Bỏ Cyclic Prefix (CP): Chỉ giữ lại N_FFT mẫu tín hiệu phía sau
        y_time = y_cp[self.n_cp : self.n_cp + self.n_fft]
        
        # 2. Normalized FFT: Biến đổi từ miền thời gian trở lại miền tần số (subcarriers)
        Y_freq = np.fft.fft(y_time) / np.sqrt(self.n_fft)
        return Y_freq
