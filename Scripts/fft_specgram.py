# %%

from speechbrain.lobes.features import Fbank
import torch
import matplotlib.pyplot as plt
from speechbrain.dataio.dataio import read_audio

signal, _ = read_audio('/home/aviadb/Desktop/thesis/code/speech_-0.82918_0.55279_-0.082918.flac')
signal, _ = read_audio('/home/aviadb/Desktop/study/525.804_Thesis_II/research/Scripts/spk1_snt1.wav')

signal.shape

#%%
signal = signal.unsqueeze(0)

signal.shape

# %%
# signal2, _ = read_audio('/home/aviadb/Desktop/thesis/code/speech_-0.82918_0.55279_-0.082918.flac')

# signal2.shape

# #%%
# signal2 = signal2.unsqueeze(0)
# signal2.shape

# # %%
# sig = signal2.squeeze(0)[:,0]
# sig.shape

#%%

from speechbrain.processing.features import spectral_magnitude
from speechbrain.processing.features import Filterbank
from speechbrain.processing.features import STFT

compute_fbanks = Filterbank(n_mels=40, n_fft=1024)
compute_STFT = STFT(sample_rate=16000, win_length=25, hop_length=10, n_fft=1024)
STFT = compute_STFT(signal)
mag = spectral_magnitude(STFT)
fbanks = compute_fbanks(mag)

freq_mat = compute_fbanks.all_freqs_mat
fb_mat = compute_fbanks.fbank_matrix

freq_mat

print(freq_mat.shape)
plt.plot(freq_mat[0,:])
plt.plot(freq_mat[1,:])
plt.show()

plt.plot(freq_mat[1,:])
plt.show()
# fb_mat

plt.figure(figsize=(8, 4), dpi=100)
plt.plot(fb_mat)
plt.grid()
plt.show()

# print(freq_mat[0,:])

#%%
print(compute_fbanks.f_central)
print(compute_fbanks.band)

plt.figure(figsize=(8.5, 5), dpi=100)
plt.plot(fb_mat)
plt.grid()
plt.show()
# print(STFT.shape)
# print(mag.shape)
# print(fbanks.shape)

# plt.imshow(fbanks.squeeze(0).t(), cmap='hot', interpolation='nearest', origin='lower')
# plt.xlabel('Time')
# plt.ylabel('Frequency')
# plt.show()

#%%


sig = signal2.unsqueeze(0)
FB = Fbank()
fffb = FB(sig)

fffb.shape
# %%
