# SỔ TAY VỠ LÒNG: HIỂU BẢN CHẤT VÀ HỌC THUỘC ĐỒ ÁN OFDM
### Đề tài: Ước lượng kênh truyền cho hệ thống OFDM trong môi trường Fading biến thiên thời gian

---

> **DÀNH CHO BẠN (NGƯỜI MỚI BẮT ĐẦU):**
> Tài liệu này được viết theo cách **bình dân nhất, dễ hiểu nhất**, dùng các hình ảnh đời thường để bạn không bị ngợp bởi toán học hay thuật ngữ kỹ thuật. Bạn chỉ cần đọc từ trên xuống dưới là sẽ **hiểu bản chất** và **thuộc lòng** để tự tin trả lời thầy cô.

---

## PHẦN 0: ĐOẠN VĂN HỌC THUỘC LÒNG (KHI THẦY HỎI TỔNG QUAN)

> **Nếu thầy hỏi:** *"Em hãy tóm tắt ngắn gọn đề tài của em làm cái gì và giải quyết bài toán gì?"*
> 
> **Học thuộc và trả lời đúng 4 câu này:**
> 1. *"Dạ thưa thầy, khi truyền tín hiệu vô tuyến (như 4G/5G/Wi-Fi), sóng phát ra bị nhà cửa, xe cộ cản trở tạo ra nhiều đường dội lại (đa đường), và người dùng di chuyển tạo ra hiệu ứng Doppler làm tín hiệu bị méo mó liên tục theo thời gian."*
> 2. *"Để máy thu giải mã được dữ liệu mà không bị lỗi, máy thu bắt buộc phải biết kênh truyền đã làm méo tín hiệu như thế nào — bài toán này gọi là **Ước lượng kênh truyền (Channel Estimation)**."*
> 3. *"Trong đồ án này, em dùng phương pháp **chèn các ký tự mẫu đã biết trước (gọi là Pilot)** vào khung dữ liệu 2D (Tần số x Thời gian)."*
> 4. *"Sau đó, em cài đặt và so sánh 2 thuật toán: thuật toán cơ bản **LS (Least Squares)** và thuật toán cải tiến **DFT-Based** để lọc nhiễu miền trễ, giúp giảm lỗi BER và tăng chất lượng truyền dẫn khi người dùng di chuyển ở các tốc độ khác nhau."*

---

## PHẦN 1: BẢN CHẤT CỦA CÁC KHÁI NIỆM (VÍ DỤ ĐỜI THƯỜNG)

### 1. Kênh đa đường (Multipath) là gì?
- **Ví dụ thực tế:** Hãy tưởng tượng bạn đứng ở đầu một con hẻm hẹp và hét to một câu sang đầu bên kia. Người ở đầu kia sẽ nghe thấy gì?
  - Tiếng nói trực tiếp đến trước.
  - Sau đó là tiếng vọng đập vào tường nhà bên trái dội lại.
  - Tiếp theo là tiếng vọng đập vào tường nhà bên phải dội lại trễ hơn một chút.
- **Trong viễn thông:** Sóng điện từ cũng vậy. Nó đập vào các tòa nhà tạo ra nhiều "đường" (gọi là các **Taps**). 
  - Đường 0: Đi thẳng (mạnh nhất, đến sớm nhất, mang 70% năng lượng).
  - Đường 1: Dội 1 lần (yếu hơn, đến trễ hơn, mang 20% năng lượng).
  - Đường 2: Dội xa hơn (yếu nhất, đến trễ nhất, mang 10% năng lượng).
  - Đây chính là cấu hình `PDP = [0.7, 0.2, 0.1]` trong bài của bạn!

### 2. OFDM là gì? Tại sao không phát một luồng thật nhanh cho xong?
- **Vấn đề:** Nếu bạn đọc một câu thật nhanh trong con hẻm vang tiếng, chữ nọ sẽ dội đè lên chữ kia $\rightarrow$ người nghe không hiểu gì cả (đây gọi là **Nhiễu xuyên ký tự - ISI**).
- **Cách giải quyết của OFDM:** Thay vì 1 người nói thật nhanh, ta chia câu nói đó cho **64 người cùng nói một lúc**, nhưng mỗi người nói **thật chậm** một từ trên một "tông giọng" (tần số con - subcarrier) riêng biệt.
  - Vì mỗi từ nói rất chậm nên tiếng vang của từ trước không thể đè lên từ sau nữa!

### 3. Cyclic Prefix (CP) là gì?
- **Hình dung:** Khi 64 người nói xong một câu, họ phải nghỉ lấy hơi một nhịp ngắn trước khi nói câu tiếp theo để tiếng vang của câu cũ tan hết trong không khí.
- **CP chính là khoảng nghỉ đó!** Người ta lấy một đoạn ngắn ở cuối câu nói dán lên đầu câu.
  - Trong bài của bạn: $N_{FFT} = 64$ mẫu, $N_{CP} = 16$ mẫu. Tổng cộng 1 symbol dài $64 + 16 = 80$ mẫu.
  - Khi máy thu nhận được âm thanh, nó chỉ việc cắt bỏ 16 mẫu CP đi $\rightarrow$ sạch bóng tiếng vang, không còn nhiễu ISI nữa!

### 4. Hiệu ứng Doppler là gì?
- **Ví dụ thực tế:** Khi một chiếc xe cứu thương hú còi chạy vụt qua bạn:
  - Khi xe chạy lại gần: Tiếng còi nghe the thé (tần số bị đẩy lên cao).
  - Khi xe chạy ra xa: Tiếng còi nghe trầm hẳn xuống (tần số bị kéo tụt xuống).
- **Trong viễn thông:** Khi điện thoại của bạn di chuyển (đi bộ, đi xe máy, ngồi tàu hỏa), tần số sóng vô tuyến bị dịch đi liên tục. Thiết bị chạy càng nhanh thì tín hiệu biến đổi càng chóng mặt (Doppler lớn).
  - `fd_norm = 0`: Đứng yên.
  - `fd_norm = 0.01`: Đi bộ.
  - `fd_norm = 0.05`: Xe máy chạy.
  - `fd_norm = 0.10`: Cao tốc/Tàu hỏa.

### 5. Pilot là gì? Tại sao phải có Pilot?
- **Hình dung:** Bạn gửi cho bạn mình một bức ảnh bị chỉnh màu (bị ngả vàng, mờ đi). Làm sao bạn của bạn biết ảnh gốc màu gì để chỉnh lại?
  - Rất đơn giản: Ở góc bức ảnh, bạn dán một **miếng dán màu trắng tinh khiết chuẩn 100%**.
  - Bạn của bạn mở ảnh ra, nhìn vào miếng dán: *"À, miếng dán trắng này đang bị ám vàng 20% và mờ đi 10%"*.
  - Từ đó, người bạn đảo ngược lại: giảm vàng 20%, tăng nét 10% cho toàn bộ bức ảnh $\rightarrow$ bức ảnh phục hồi hoàn hảo!
- **Miếng dán chuẩn đó chính là PILOT!**
  - Cả máy phát và máy thu đều biết trước: *"Tại vị trí này, tao sẽ phát giá trị chuẩn là $1.0 + 0j$ nhé"*.
  - Máy thu đo vị trí đó thấy nhận được số khác $\rightarrow$ suy ra kênh truyền đã bẻ cong tín hiệu thế nào!

---

## PHẦN 1.5: VÍ DỤ MINH HỌA TRỰC QUAN: HÀNH TRÌNH BIẾN ĐỔI CỦA 2 BIT DỮ LIỆU [0, 1] QUA TỪNG GIAI ĐOẠN

Để hiểu rõ dữ liệu biến đổi như thế nào qua từng khâu trong code, hãy theo dõi hành trình của đúng **2 bit dữ liệu `[0, 1]`** từ máy phát, bay qua không gian bị méo mó, đến lúc máy thu khôi phục lại nguyên vẹn:

```
[MÁY PHÁT - TX]
  Bước 1: Chuỗi bit đầu vào: [0, 1]
  Bước 2: Điều chế QPSK -> Biến thành toạ độ phức: (-0.707 + 0.707j)
  Bước 3: Ghép Pilot -> Xếp Data cạnh Pilot chuẩn (1.0 + 0j) vào ma trận 64 sóng mang
  Bước 4: IFFT (64 điểm) -> Chuyển thành chuỗi 64 rung động thời gian x[n]
  Bước 5: Ghép Cyclic Prefix (CP) -> Lấy 16 mẫu đuôi dán lên đầu thành 80 mẫu
        │
════════╪════════════ KHÔNG KHÍ (MÔI TRƯỜNG TRUYỀN SÓNG) ════════════
        │  1. Bị dội tường (Multipath): 3 đường [0.7, 0.2, 0.1] đè lên nhau
        │  2. Xe chạy (Doppler): Tín hiệu bị xoay góc pha
        │  3. Nhiễu môi trường: Cộng thêm nhiễu trắng AWGN
        ▼
[MÁY THU - RX]
  Bước 6: Cắt bỏ CP -> Vứt bỏ 16 mẫu đầu (nơi hứng trọn tiếng vang) còn 64 mẫu
  Bước 7: FFT (64 điểm) -> Đưa về miền tần số Y[k] (tín hiệu bị méo xẹo)
  Bước 8: Ước lượng kênh & Cân bằng:
          - Đo Pilot xem bị méo bao nhiêu -> Đoán được kênh truyền H
          - Lấy tín hiệu thu Y chia cho H để "nắn thẳng" lại dữ liệu
  Bước 9: Giải mã chòm sao (Demapping) -> Khôi phục chính xác 100% cặp bit [0, 1]!
```

---

### BƯỚC ĐI CỤ THỂ VỚI CON SỐ THỰC TẾ:

#### 🔹 Bước 1: Có 2 bit dữ liệu cần gửi
- Điện thoại của bạn muốn gửi đi 2 bit: **`[0, 1]`**.
- Ăng-ten không thể phát trực tiếp số 0 hay 1 ra không gian, mà phải biến thành **biên độ và góc pha của sóng điện từ**.

#### 🔹 Bước 2: Điều chế QPSK (Ánh xạ chòm sao IQ)
QPSK gom cứ 2 bit thành 1 tọa độ phức:
- `[0, 0]` $\rightarrow$ Góc phần tư I: $(+0.707 + 0.707j)$
- **`[0, 1]` $\rightarrow$ Góc phần tư II: $(-0.707 + 0.707j)$** *(chính là 2 bit của bạn!)*
- `[1, 1]` $\rightarrow$ Góc phần tư III: $(-0.707 - 0.707j)$
- `[1, 0]` $\rightarrow$ Góc phần tư IV: $(+0.707 - 0.707j)$
> **Kết quả:** 2 bit `[0, 1]` đã biến thành con số phức: **$X_{data} = -0.707 + 0.707j$**.

#### 🔹 Bước 3: Chèn Pilot (Cắm mốc chuẩn)
Hệ thống chia băng thông thành 64 sóng mang con ($k = 0, 1, 2 \dots 63$):
- Sóng mang $k = 0$ (Pilot): Gán giá trị chuẩn **$X[0] = 1.0 + 0j$** (Cả 2 bên đều biết trước).
- Sóng mang $k = 1$ (Data của bạn): Gán giá trị **$X[1] = -0.707 + 0.707j$**.
- Sóng mang $k = 2, 3 \dots$: Điền các bit dữ liệu khác.
- Sóng mang $k = 4$ (Pilot tiếp theo): Lại cắm mốc **$X[4] = 1.0 + 0j$**.

#### 🔹 Bước 4: Biến đổi IFFT (64 điểm)
- Biến 64 con số phức ở miền tần số thành 64 mẫu tín hiệu thời gian rung rinh $x[0], x[1] \dots x[63]$.

#### 🔹 Bước 5: Ghép Cyclic Prefix (CP dài 16 mẫu)
- Lấy 16 mẫu cuối ($x[48] \dots x[63]$) dán lên đầu chuỗi.
- Tổng cộng chuỗi dài $16 + 64 = 80$ mẫu. Ăng-ten phát chuỗi 80 mẫu này ra không khí.

#### 🔹 Bước 6: Bay qua không gian (Bị kênh truyền làm méo)
- Sóng dội qua các tòa nhà và thiết bị di chuyển khiến tín hiệu bị suy giảm và xoay góc pha.
- Giả sử kênh truyền tác động một hệ số làm méo là $H \approx 0.5 + 0.3j$, cộng thêm nhiễu trắng môi trường.

#### 🔹 Bước 7: Máy thu nhận tín hiệu và chạy FFT
- Máy thu nhận được 80 mẫu sóng.
- **Cắt bỏ CP:** Vứt bỏ 16 mẫu đầu tiên (nơi dính toàn bộ tiếng vang dội tường của ký tự trước).
- **FFT (64 điểm):** Đưa 64 mẫu còn lại trở về miền tần số $Y[k]$.
- Lúc này máy thu kiểm tra giá trị nhận được:
  - **Tại sóng mang 0 (Pilot):** Ban đầu phát $1.0$, nhưng thu được:
    $$Y[0] \approx 0.51 + 0.29j \quad (\text{bị méo đi!})$$
  - **Tại sóng mang 1 (Data của bạn):** Ban đầu phát $-0.707 + 0.707j$, nhưng thu được:
    $$Y[1] \approx -0.57 + 0.15j \quad (\text{méo mó hoàn toàn, nếu đọc ngay sẽ sai bét!})$$

#### 🔹 Bước 8: Ước lượng kênh & Cân bằng Zero-Forcing (Phép màu ở đây!)
1. **Đo độ méo tại Pilot (Ước lượng LS):**
   $$\hat{H}[0] = \frac{Y[0]}{X_{pilot}} = \frac{0.51 + 0.29j}{1.0} = 0.51 + 0.29j$$
   Máy thu phát hiện: *"Kênh truyền đang làm méo tín hiệu theo hệ số $0.51 + 0.29j$!"*
2. **Nội suy sang sóng mang 1:**
   Vì sóng mang 1 nằm ngay sát sóng mang 0, nên kênh tại sóng mang 1 cũng bị méo khoảng $\hat{H}[1] \approx 0.51 + 0.29j$.
   *(Nếu dùng thêm bộ lọc **DFT-Based**, máy thu sẽ gạt bỏ bớt nhiễu ngoài vùng trễ, giúp $\hat{H}[1]$ chính xác hơn nữa!)*
3. **Cân bằng kênh Zero-Forcing (Nắn thẳng lại dữ liệu):**
   Máy thu lấy dữ liệu thu được $Y[1]$ chia ngược lại cho kênh $\hat{H}[1]$:
   $$\hat{X}[1] = \frac{Y[1]}{\hat{H}[1]} = \frac{-0.57 + 0.15j}{0.51 + 0.29j} \approx -0.71 + 0.69j$$

> 🎯 **Hãy nhìn con số kỳ diệu này:**  
> Con số sau khi nắn là **$-0.71 + 0.69j$**, nó đã quay trở về **cực kỳ sát với giá trị ban đầu** là **$-0.707 + 0.707j$**!

#### 🔹 Bước 9: Giải mã chòm sao (Demapping) ra lại bit
- Máy thu chấm điểm $-0.71 + 0.69j$ lên mặt phẳng:
  - Điểm này nằm trọn trong **Góc phần tư thứ II** (Trục thực âm, trục ảo dương).
  - Tra bảng chòm sao: Góc phần tư thứ II chính là cặp bit **`[0, 1]`**!
- **Kết luận:** Máy thu khôi phục chính xác 100% cặp bit `[0, 1]` ban đầu mà không có bất kỳ lỗi nào!

---

## PHẦN 2: THẦY BẢO SỬA 4 LỖI GÌ? (NÓI BẰNG TIẾNG VIỆT DỄ HIỂU)

Khi thầy hỏi: *"Đợt trước thầy nhận xét 4 lỗi, em đã hiểu và sửa như thế nào?"* $\rightarrow$ Hãy trả lời rành mạch từng lỗi như sau:

---

### LỖI 1: Nhầm đơn vị thời gian trong mô hình Jakes (Chia cho 80)
- **Bản chất lỗi:** 
  - Tần số Doppler chuẩn hóa $f_{d,norm}$ được tính theo đơn vị **1 ký tự OFDM** ($T_{sym}$).
  - Nhưng 1 ký tự OFDM trong code gồm có 64 mẫu dữ liệu + 16 mẫu CP = **80 mẫu thời gian**.
  - Trước đây, trong code lấy thẳng chỉ số mẫu ($0, 1, 2, 3 \dots$) làm thời gian $t$. Nghĩa là cứ qua 1 mẫu thời gian thì code lại tưởng là đã trôi qua 1 ký tự OFDM!
  - Kết quả: **Thời gian bị chạy nhanh gấp 80 lần!** Xe đạp biến thành tên lửa siêu thanh!
- **Cách em sửa:**
  ```python
  # Trong file channel.py:
  # Lấy chỉ số mẫu chia cho 80 (symbol_len) để quy đổi đúng về đơn vị ký tự OFDM
  t = float(sample_index) / float(self.symbol_len)
  ```
- **Kiểm chứng:** Em đo hàm tự tương quan thời gian của kênh, đường cong đo được khớp hoàn toàn 100% với hàm toán học Bessel $J_0(2\pi f_d \Delta t)$ của giáo trình.

---

### LỖI 2: Tính sai công suất kênh (Hệ số $\sqrt{4N_0 + 1}$)
- **Bản chất lỗi:**
  - Tổng công suất của 3 đường kênh truyền phải bằng đúng $1.0$ ($0.7 + 0.2 + 0.1 = 1.0$) thì việc đo đạc tỷ lệ tín hiệu trên nhiễu (SNR) mới chính xác.
  - Công thức mô hình Jakes tổng hợp kênh từ $N_0 = 8$ sóng sin. Phần thực có $2N_0 + 1$ thành phần, phần ảo có $2N_0$ thành phần $\rightarrow$ Tổng công suất thô là $(2N_0 + 1) + 2N_0 = 4N_0 + 1$.
  - Trước đây, code chia nhầm cho $\sqrt{2N_0 + 1}$, làm cho công suất kênh bị thổi phồng lên gần gấp rưỡi (sai lệch 47%)!
- **Cách em sửa:**
  ```python
  # Trong file channel.py:
  # Sửa mẫu số thành đúng công thức giải tích 4*N0 + 1
  scale = np.sqrt(p_l / (4.0 * self.num_sinusoids + 1.0))
  h_tap = scale * (h_I + 1j * h_Q)
  ```
- **Kiểm chứng:** Sau khi sửa, em chạy 5000 lần thử, đo công suất 3 đường lần lượt ra đúng $0.7$, $0.2$, $0.1$, tổng công suất bằng $1.0$ và biên độ tuân theo đúng phân bố Rayleigh.

---

### LỖI 3: Dùng lại một kênh cho 200 frame trong mô phỏng Monte Carlo
- **Bản chất lỗi:**
  - Mô phỏng Monte Carlo là phương pháp thử nghiệm nhiều lần (200 frame) trong các điều kiện ngẫu nhiên khác nhau để lấy kết quả trung bình đáng tin cậy.
  - Trước đây, code tạo ra 1 đối tượng kênh ở ngoài vòng lặp, rồi dùng đi dùng lại kênh đó cho cả 200 frame.
  - Việc này giống như thử độ bền của xe ô tô trên 200 đoạn đường, nhưng đoạn đường nào cũng có đúng 1 ổ gà ở đúng 1 vị trí! Kết quả bị thiên lệch (bias), mất tính ngẫu nhiên.
- **Cách em sửa:**
  ```python
  # Trong file run_thesis_experiments.py:
  # Chuyển lệnh tạo channel vào TRONG vòng lặp từng frame:
  for _ in range(num_frames):
      channel = TimeVaryingRayleighChannel(...) # Mỗi frame là 1 kênh ngẫu nhiên mới!
      r = simulate_frame(...)
  ```
- **Kiểm chứng:** 200 frame bây giờ trải qua 200 trạng thái kênh độc lập khác nhau, số liệu thống kê chuẩn xác tuyệt đối.

---

### LỖI 4: Tự bịa số liệu BER = $10^{-6}$ khi không có lỗi
- **Bản chất lỗi:**
  - Ở mức SNR cao (ví dụ 24 dB), tín hiệu rất rõ nên trong 200 frame không có bit nào bị lỗi ($\text{Số bit lỗi} = 0 \rightarrow \text{BER} = 0$).
  - Đồ thị dùng thang đo logarit thì không vẽ được số 0. Code cũ dùng lệnh `max(ber, 1e-6)` để ép số 0 thành $10^{-6}$ (một phần triệu) để vẽ cho đẹp.
  - Thầy phê bình: Mô phỏng chỉ có khoảng vài chục nghìn bit, không thể có chuyện đo được độ chính xác 1 phần triệu! Làm vậy là "ngụy tạo số liệu".
- **Cách em sửa:**
  - Giữ nguyên số bit lỗi thực tế chia cho tổng số bit: `ber = errors / total_bits`.
  - Viết hàm vẽ đồ thị `safe_semilogy`: Điểm nào có $\text{BER} = 0$ thì đổi thành `NaN` (bỏ trống điểm đó trên đồ thị, không vẽ sàn giả).

---

## PHẦN 3: 2 THUẬT TOÁN ƯỚC LƯỢNG KÊNH HOẠT ĐỘNG THẾ NÀO?

Trong đồ án có 2 thuật toán chính mà bạn cần phân biệt rành rọt:

### 1. Thuật toán LS (Least Squares - Bình phương tối thiểu)
- **Cách làm:** Cực kỳ đơn giản.
  - Bước 1: Tại những ô có cắm cờ Pilot, lấy giá trị thu được chia cho giá trị Pilot gốc: 
    $$\hat{H}_{LS} = \frac{Y_{pilot}}{X_{pilot}}$$
  - Bước 2: Tại những ô dữ liệu không có Pilot, dùng phương pháp **Nội suy (nối các điểm lại với nhau theo đường thẳng)** để đoán giá trị kênh.
- **Ưu điểm:** Chạy siêu nhanh (chỉ mất ~0.9 ms cho 1 frame).
- **Nhược điểm:** Kém cỏi trước nhiễu. Nếu tín hiệu bị nhiễu AWGN đè lên thì thuật toán LS nhận luôn cả nhiễu vào kênh ước lượng, dẫn đến giải mã bit bị sai nhiều.

### 2. Thuật toán DFT-Based (Bộ lọc khử nhiễu miền trễ)
- **Bản chất (Mẹo thông minh):**
  - Tín hiệu kênh truyền thật chỉ xuất hiện trongL,,,,,,,,,,,,     một khoảng trễ rất ngắn ở đầu (tối đa là 16 mẫu CP).
  - Nhưng nhiễu trắng AWGN thì lan tràn khắp nơi trong toàn bộ 64 mẫu.
- **Cách làm 3 bước:**
  1. Lấy kết quả từ LS, dùng phép biến đổi **IFFT** để chuyển về miền thời gian trễ.
  2. Dùng một **"cái kéo" (Cửa sổ Window)**: Giữ lại các mẫu đầu tiên (nơi có sóng thật), còn tất cả các mẫu phía sau (chỉ toàn là rác và nhiễu) thì gạt phăng đi, gán bằng $0$!
  3. Dùng phép biến đổi **FFT** đưa ngược trở lại miền tần số.
- **Kết quả:** Nhiễu bị triệt tiêu phần lớn! Đường cong MSE giảm sâu, BER giảm rõ rệt, chất lượng hệ thống vượt trội so với LS mà thời gian chạy chỉ nhỉnh hơn một chút (~1.3 ms so với 0.9 ms).

---

## PHẦN 4: HƯỚNG DẪN ĐỌC TỪNG FILE CODE (DÀNH CHO NGƯỜI CHƯA BIẾT CODE)

Thư mục code của bạn gồm 5 file chính:

### 1. `config.py` — Bảng thiết lập hằng số
- Chứa các con số cố định của hệ thống:
  - `N_FFT = 64`: 64 sóng mang con.
  - `N_CP = 16`: Tiền tố lặp dài 16 mẫu.
  - `CHANNEL_PDP = [0.7, 0.2, 0.1]`: Kênh 3 đường với công suất 70%, 20%, 10%.
  - `NUM_SYMBOLS_PER_FRAME = 14`: 1 frame có 14 symbol.
  - `SNR_DB_RANGE = np.arange(0, 26, 2)`: Quét SNR từ 0 đến 24 dB, mỗi bước 2 dB.

### 2. `ofdm_system.py` — Cỗ máy phát và thu tín hiệu
- **Bên phát (TX):**
  - `map_bits`: Biến chuỗi bit $0, 1$ thành các toạ độ phức (QPSK: mỗi 2 bit thành 1 điểm).
  - `insert_pilots_2d`: Tạo bảng caro $64 \times 14$, cắm các điểm Pilot chuẩn vào các vị trí quy định, chỗ trống còn lại thì điền dữ liệu người dùng.
  - `modulate_frame`: Chạy IFFT và dán đầu CP 16 mẫu cho từng symbol.
- **Bên thu (RX):**
  - `demodulate_frame`: Cắt bỏ 16 mẫu CP rồi chạy FFT đưa về tần số.
  - `equalize_and_decode_frame`: Lấy tín hiệu thu chia cho kênh ước lượng ($\hat{X} = Y / \hat{H}$) để nắn thẳng lại tín hiệu (Zero-Forcing), sau đó dò tìm xem điểm đó gần số $0$ hay số $1$ nhất để trả về chuỗi bit gốc.

### 3. `channel.py` — Môi trường không khí giả lập
- Lớp `TimeVaryingRayleighChannel`:
  - Dùng các hàm $\cos, \sin$ của mô hình Jakes để tạo ra độ trôi pha Doppler.
  - Làm tín hiệu bị biến dạng theo thời gian thực khi xe di chuyển.
  - Cộng thêm tiếng ồn ngẫu nhiên (nhiễu trắng AWGN) tùy theo mức SNR mạnh hay yếu.

### 4. `estimator.py` — Bộ não ước lượng kênh
- Lớp `LSEstimator`:
  - `estimate_ls_2d`: Thực hiện phép chia tại Pilot rồi nội suy 2D (tần số và thời gian).
  - `estimate_dft_2d`: Nhận kết quả LS, IFFT về miền trễ, dùng cửa sổ cắt bỏ nhiễu ngoài vùng trễ, rồi FFT ngược lại.

### 5. `run_thesis_experiments.py` — File chạy tự động 4 thí nghiệm
- Khi bạn chạy lệnh `py run_thesis_experiments.py`, nó sẽ lần lượt thực hiện:
  - **Thí nghiệm 1:** So sánh các thuật toán ở kênh tĩnh và kênh chậm $\rightarrow$ Vẽ ra ảnh `fig1` (MSE) và `fig2` (BER).
  - **Thí nghiệm 2:** Thử nghiệm xe chạy từ chậm đến nhanh ($f_d$ từ $0$ đến $0.10$) $\rightarrow$ Vẽ ra ảnh `fig3`.
  - **Thí nghiệm 3:** Thử nghiệm cắm Pilot dày hay thưa $(D_f, D_t)$ để xem cái nào lợi hơn giữa tốc độ truyền và tỷ lệ lỗi $\rightarrow$ Vẽ ra ảnh `fig4` và `fig5`.
  - **Thí nghiệm 4:** Bấm giờ xem thuật toán nào chạy nhanh hơn và xuất bảng khuyến nghị pilot ra file CSV.

---

## PHẦN 5: BỘ CÂU HỎI VẤN ĐÁP "CỨU MẠNG" KHI THẦY HỎI

Hãy đọc đi đọc lại 4 câu hỏi này trước khi bước vào phòng bảo vệ:

### ❓ Câu 1: "Tại sao em lại vẽ đồ thị theo SNR từ 0 đến 24 dB?"
> **Trả lời:** *"Dạ thưa thầy, SNR thể hiện chất lượng đường truyền (tỷ lệ tín hiệu trên nhiễu). Ở 0 dB, nhiễu rất mạnh tương đương với tín hiệu, nên BER rất cao (~20-30%). Khi tăng SNR lên 24 dB, tín hiệu vượt trội hoàn toàn so với nhiễu, giúp các thuật toán ước lượng chính xác hơn và BER giảm xuống tiệm cận 0."*

### ❓ Câu 2: "Tại sao đường cong DFT với L = 3 lại tốt hơn L = 16?"
> **Trả lời:** *"Dạ thưa thầy:
> - $L = 16$ là trường hợp thực tế khi máy thu chỉ biết kênh nằm trong độ dài CP ($N_{CP} = 16$).
> - $L = 3$ là trường hợp lý tưởng khi máy thu biết chính xác kênh chỉ có 3 đường phản xạ ($PDP = [0.7, 0.2, 0.1]$). Khi đó, cửa sổ $L = 3$ gạt bỏ sạch sẽ nhiễu ở 61 mẫu còn lại, trong khi cửa sổ $L = 16$ vẫn bị lọt lại nhiễu của 13 mẫu dư thừa. Vì vậy đường $L = 3$ luôn đạt hiệu năng tốt nhất làm mốc chuẩn lý thuyết ạ."*

### ❓ Câu 3: "Khi xe chạy tốc độ rất cao ($f_{d,norm} = 0.10$), tại sao tăng SNR lên cao mà BER vẫn không giảm (bị sàn lỗi Error Floor)?"
> **Trả lời:** *"Dạ thưa thầy, khi di chuyển quá nhanh, kênh truyền biến đổi ngay trong lòng một ký tự OFDM làm mất tính trực giao giữa các sóng mang con (gây hiện tượng nhiễu ICI). Lúc này, lỗi sinh ra không phải do nhiễu trắng AWGN nữa mà do kênh biến đổi quá nhanh khiến bộ ước lượng không theo kịp. Do đó dù có tăng công suất phát (tăng SNR) thì lỗi do Doppler gây ra vẫn tồn tại, tạo thành sàn lỗi bão hòa ạ."*

### ❓ Câu 4: "Từ kết quả đồ án, em rút ra kết luận gì về việc chọn khoảng cách Pilot?"
> **Trả lời:** *"Dạ thưa thầy, có một sự đánh đổi (Trade-off) rõ rệt:
> - Nếu cắm Pilot quá dày ($D_f=2, D_t=1$): Ước lượng rất chính xác, chịu được tốc độ cao, nhưng Pilot chiếm tới 25% băng thông, làm giảm tốc độ truyền dữ liệu có ích.
> - Nếu cắm Pilot quá thưa ($D_f=8, D_t=4$): Tiết kiệm băng thông nhưng khi xe chạy nhanh sẽ bị lỗi mất dữ liệu.
> - **Khuyến nghị của em:** Với người đi bộ chậm chỉ cần dùng $D_f=4, D_t=4$ (Overhead ~8.9%); còn với môi trường di chuyển đô thị cần dùng $D_f=4, D_t=2$ (Overhead ~14.3%) để vừa tiết kiệm băng thông vừa đảm bảo chất lượng kết nối ạ."*

---
*Tài liệu này được biên soạn riêng để bạn hiểu tận gốc bản chất và tự tin 100% trong buổi báo cáo!*
