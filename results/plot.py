import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("sin.csv")

plt.figure(figsize=(10, 5))
plt.plot(df["t"], df["x"])

plt.title("Sine wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.grid(True)

plt.savefig("sine_wave.png")

plt.show()

