import pandas as pd 
import matplotlib.pyplot as plt

df = pd.read_csv("quantized.csv")

plt.figure(figsize=(10, 5))
plt.plot(df["t"], df["original"], label="Original")
plt.plot(df["t"], df["q_restored"], label="Quantized")
plt.title("Original vs 16-bit quantized sine wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("sine_wave_q.png")
plt.show()


plt.figure(figsize=(10, 5))
plt.plot(df["t"], df["error"])
plt.xlabel("Time [s]")
plt.ylabel("Error")
plt.title("Quantization error")
plt.grid(True)
plt.tight_layout()
plt.savefig("quantization_error.png")
plt.show()

