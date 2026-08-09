import numpy as np
from config import N_FFT, N_CP, MODULATION, PILOT_SPACING, PILOT_VALUE, CHANNEL_PDP
from channel import MultipathRayleighChannel
from ofdm_system import OFDMSystem
from estimator import LSEstimator

def test_noiseless_loopback():
    print("--- Running Test 1: Noiseless Ideal Loopback ---")
    ofdm = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, modulation=MODULATION)
    estimator = LSEstimator()
    channel = MultipathRayleighChannel([1.0])
    
    # Ideal flat channel h = [1.0]
    h_taps = np.array([1.0 + 0.0j])
    H_true = channel.get_frequency_response(h_taps, N_FFT)
    
    # Generate random bits
    num_data_subcarriers = ofdm.num_data
    bits_per_symbol = ofdm.bits_per_symbol
    tx_bits = np.random.randint(0, 2, num_data_subcarriers * bits_per_symbol)
    
    # 1. Map bits
    tx_data_symbols = ofdm.map_bits(tx_bits)
    
    # 2. Insert pilots
    X_freq = ofdm.insert_pilots(tx_data_symbols)
    
    # 3. Modulate (IFFT + CP)
    x_cp = ofdm.modulate(X_freq)
    
    # 4. Pass through noiseless channel
    y_cp = channel.apply_channel(x_cp, h_taps, snr_db=100.0) # High SNR = noiseless
    Y_freq = ofdm.demodulate(y_cp)
    
    # 5. Estimate channel
    H_est = estimator.estimate(Y_freq, ofdm.pilot_indices, PILOT_VALUE, N_FFT)
    
    # 6. Equalize data
    Y_data = Y_freq[ofdm.data_indices]
    H_est_data = H_est[ofdm.data_indices]
    rx_data_symbols = estimator.equalize(Y_data, H_est_data)
    
    # 7. Demap bits
    rx_bits = ofdm.demap_symbols(rx_data_symbols)
    
    bit_errors = np.sum(tx_bits != rx_bits)
    mse = estimator.compute_mse(H_true, H_est)
    
    print(f"Noiseless Loopback Bit Errors: {bit_errors}")
    print(f"Noiseless Loopback Channel MSE: {mse:.6e}")
    assert bit_errors == 0, "Loopback should have 0 bit errors!"
    assert mse < 1e-10, "Loopback MSE should be practically 0!"
    print(">>> Test 1 PASSED!\n")

def test_multipath_rayleigh():
    print("--- Running Test 2: Rayleigh Channel Simulation ---")
    ofdm = OFDMSystem(n_fft=N_FFT, n_cp=N_CP, pilot_spacing=PILOT_SPACING, modulation=MODULATION)
    channel = MultipathRayleighChannel(CHANNEL_PDP)
    estimator = LSEstimator()
    
    h_taps = channel.generate_taps()
    H_true = channel.get_frequency_response(h_taps, N_FFT)
    
    tx_bits = np.random.randint(0, 2, ofdm.num_data * ofdm.bits_per_symbol)
    tx_data = ofdm.map_bits(tx_bits)
    X_freq = ofdm.insert_pilots(tx_data)
    x_cp = ofdm.modulate(X_freq)
    
    snr_db = 20.0
    y_cp = channel.apply_channel(x_cp, h_taps, snr_db)
    Y_freq = ofdm.demodulate(y_cp)
    
    H_est = estimator.estimate(Y_freq, ofdm.pilot_indices, PILOT_VALUE, N_FFT)
    
    # Perfect CSI equalization
    rx_data_perfect = estimator.equalize(Y_freq[ofdm.data_indices], H_true[ofdm.data_indices])
    bits_perfect = ofdm.demap_symbols(rx_data_perfect)
    err_perfect = np.sum(tx_bits != bits_perfect)
    
    # LS CSI equalization
    rx_data_ls = estimator.equalize(Y_freq[ofdm.data_indices], H_est[ofdm.data_indices])
    bits_ls = ofdm.demap_symbols(rx_data_ls)
    err_ls = np.sum(tx_bits != bits_ls)
    
    mse = estimator.compute_mse(H_true, H_est)
    
    print(f"SNR: {snr_db} dB")
    print(f"Perfect CSI Bit Errors: {err_perfect} / {len(tx_bits)}")
    print(f"LS Estimation Bit Errors: {err_ls} / {len(tx_bits)}")
    print(f"LS Channel MSE: {mse:.6e}")
    print(">>> Test 2 PASSED!\n")

if __name__ == '__main__':
    test_noiseless_loopback()
    test_multipath_rayleigh()
