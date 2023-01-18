import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

X_00 = pd.read_excel(r'C:/Users/abhij/Desktop/WaveModel/StiffPlate/UT2DLA0S16-1T0003_M00025.xlsx').to_numpy()
X_05 = pd.read_excel(r'C:/Users/abhij/Desktop/WaveModel/Dam5mm/UT2DLA0S16-1T0003_M00025.xlsx').to_numpy()
X_10 = pd.read_excel(r'C:/Users/abhij/Desktop/WaveModel/Dam10mm/UT2DLA0S16-1T0003_M00025.xlsx').to_numpy()
X_15 = pd.read_excel(r'C:/Users/abhij/Desktop/WaveModel/Dam15mm/UT2DLA0S1-16T0003_M00025.xlsx').to_numpy()
X_20 = pd.read_excel(r'C:/Users/abhij/Desktop/WaveModel/Dam20mm/UT2DLA0S1-16T0003_M00025.xlsx').to_numpy()

sn = 10
x_sen=np.array([X_00[:,(19-sn)], X_05[:,sn], X_10[:,(19-sn)], X_15[:,sn], X_20[:,sn]])
X_sen = np.transpose(x_sen)
# normalized function
def normalize_fn(a):
    norm = []
    for j in range(len(a)):
        i = a[j]
        nor = 2*(i-np.min(a))/(np.max(a)-np.min(a))-1
        norm.append(nor)
    return norm
X_norm = []
for i in range(5):
    x_norm = normalize_fn(X_sen[:,i])
    X_norm.append(x_norm)
X_norm=np.array(X_norm)
X_norm = np.transpose(X_norm)
# FFT of normalized signal
Fs = 1e8                     #sampling Frequency
samplingInterval = 1.0 / Fs  # sample Interval
start_time = 0.0
end_time = 0.0003
time = np.arange(start_time, end_time, samplingInterval);
# Frequency domain representation
fourierTransform = np.fft.fft(X_norm[:,1])/len(X_norm[:,1])           # Normalize amplitude
fourierTransform = fourierTransform[range(int(len(X_norm[:,1])/2))]   # Exclude sampling frequency
tpCount     = len(X_norm[:,1])
values      = np.arange(int(tpCount/2))
timePeriod  = tpCount/Fs
frequencies = values/timePeriod
plt.plot(frequencies, abs(fourierTransform))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.xlim(0,1e6)
plt.ylim(0,0.01)
plt.show()


# FIR filter design
# Create a FIR filter and apply it to signal.
#------------------------------------------------
from numpy import sin, arange, pi
from scipy.signal import lfilter, firwin
# The Nyquist rate of the signal.
sample_rate = 1e8
nyq_rate = sample_rate / 2.
# The cutoff frequency of the filter: 6KHz
cutoff_hz = 295000
# Length of the filter (number of coefficients, i.e. the filter order + 1)
numtaps = 10001
# Use firwin to create a lowpass FIR filter
fir_coeff = firwin(numtaps, cutoff_hz/nyq_rate)
# Use lfilter to filter the signal with the FIR filter
signal =X_norm[:,1]
filtered_signal = lfilter(fir_coeff, 1.0,signal )
#------------------------------------------------
# Plot the original and filtered signals.
#------------------------------------------------
# The first N-1 samples are "corrupted" by the initial conditions
warmup = numtaps - 1
# The phase delay of the filtered signal
delay = (warmup / 2) / sample_rate
#------------------------------------------------
# Print values
#------------------------------------------------
def print_values(label, values):
    var = "float32_t %s[%d]" % (label, len(values))
    print("%-30s = {%s}" % (var, ', '.join(["%+.10f" % x for x in values])))
 
#print_values('signal', signal)
#print_values('fir_coeff', fir_coeff)
#print_values('filtered_signal', filtered_signal)
plt.plot(filtered_signal)
plt.show()
plt.plot(signal)
plt.show()
# FFT of filtered signal 
Fs = 1e8                     #sampling Frequency
samplingInterval = 1.0 / Fs  # sample Interval
start_time = 0.0
end_time = 0.0003
time = np.arange(start_time, end_time, samplingInterval);
# Frequency domain representation
fourierTransform = np.fft.fft(filtered_signal)/len(filtered_signal)      # Normalize amplitude
fourierTransform = fourierTransform[range(int(len(filtered_signal)/2))]   # Exclude sampling frequency
tpCount     = len(X_norm[:,1])
values      = np.arange(int(tpCount/2))
timePeriod  = tpCount/Fs
frequencies = values/timePeriod
plt.plot(frequencies, abs(fourierTransform))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.xlim(0,1e6)
plt.ylim(0,0.01)
plt.show()


import pywt
x = np.arange(512)
y = np.sin(2*np.pi*x/32)
coef, freqs=pywt.cwt(X_norm[:,1],np.arange(1,129),'gaus1')
plt.matshow(coef) # doctest: +SKIP
plt.show() # doctest: +SKIP
t = np.linspace(-1, 1, 200, endpoint=False)
sig  = np.cos(2 * np.pi * 7 * t) + np.real(np.exp(-7*(t-0.4)**2)*np.exp(1j*2*np.pi*2*(t-0.4)))
widths = np.arange(1, 31)
cwtmatr, freqs = pywt.cwt(X_norm[:,1], widths, 'mexh')
plt.imshow(cwtmatr, extent=[-1, 1, 1, 31], cmap='PRGn', aspect='auto',
           vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())  # doctest: +SKIP
plt.show() # doctest: +SKIP
