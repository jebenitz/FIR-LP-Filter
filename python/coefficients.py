import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

fs = 100000 #sample freq
fc = 12500 #cutoff freq

M = 16  #number of even taps
omegac = (fc/fs) *2*np.pi


def hd (n,omegac,M):
    return np.sin(omegac*(n-(M-1)/2))/(np.pi*(n-(M-1)/2))

n = np.arange(M)

h = hd(n, omegac, M ) * np.hamming(M) #Hamming window

h = h / np.sum(h)  #normalized value forces H(0)=1
h0 = np.sum(h) 

print ("coefficents: ",h)  #coefficients values
print ("H(0)=",h0 ) # H(0) = 1.0 (check)


#### Q1.15 formating and hex conversion ####
q = 2**15

hq32 = np.round(q*h).astype(np.int32) 
assert np.all(hq32 >= -32768) and np.all(hq32 <= 32767), "overflow/wrap detected in Q1.15 coefficients"
hq15 = hq32.astype(np.int16)  

hhex = [format(int(v) & 0xFFFF,'04X') for v in hq15]

print("coefficients Q15:", hq15) 
print("HEX: ", hhex)

hround = hq15 / (2**15) #rounded values for  plots comparison 
print("Rounded: ", hround)


#### saving the coefficients list ####

np.savetxt("python/coff.csv",h, fmt="%f", delimiter= ",")


with open("python/coeffhex.vh", "w") as f:
    for i, val in enumerate(hhex):
        f.write(f"localparam COEFF_{i} = 16'h{val};\n")


### Freq response plot ###
freqs, H = signal.freqz(h, worN=512) #worN freq points between 0 and pi values
freqs, HR = signal.freqz(hround, worN=512)


plt.plot(freqs, 20*np.log10(np.abs(H)),label ="|H| calculated")
plt.plot(freqs, 20*np.log10(np.abs(HR)),label ="|H| truncated")

plt.annotate("$\omega_c$", xy =[omegac+0.05,-88] )
plt.axvline(x=omegac, linestyle='--', linewidth=1, color="gray")


plt.xlabel("$\omega$ (rad/sample)")
plt.ylabel("|H| (dB)")
plt.title("FIR frequency response")
plt.grid(True)
plt.legend()
plt.savefig("python/response.png")


