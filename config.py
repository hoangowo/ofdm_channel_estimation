import numpy as np

# OFDM System Parameters
N_FFT = 64              # Total number of subcarriers (FFT size)
N_CP = 16               # Cyclic Prefix length in samples
MODULATION = 'QPSK'     # Modulation scheme ('QPSK' or '16QAM')

# Pilot Parameters
PILOT_SPACING = 4       # Subcarrier distance between adjacent pilots (D_f)
PILOT_VALUE = 1.0 + 0.0j# Pilot symbol value (BPSK)

# Simulation Parameters
SNR_DB_RANGE = np.arange(0, 26, 2)  # SNR range in dB (0, 2, 4, ..., 24)
NUM_OFDM_SYMBOLS = 2000              # Number of OFDM symbols per SNR level

# Channel Parameters
CHANNEL_PDP = np.array([0.7, 0.2, 0.1])  # 3-tap normalized Power Delay Profile
