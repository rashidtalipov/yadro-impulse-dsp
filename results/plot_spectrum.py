import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FS = 100.0
FS2 = 200.0
EPS = 1e-12


def compute_spectrum_raw(x, fs):
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    xw = x * np.hanning(len(x))
    spectrum = np.fft.rfft(xw)
    freqs = np.fft.rfftfreq(len(xw), d=1.0 / fs)
    return freqs, np.abs(spectrum)


def compute_spectrum_db(x, fs):
    freqs, magnitude = compute_spectrum_raw(x, fs)
    return freqs, 20.0 * np.log10(magnitude / np.max(magnitude) + EPS)


def plot_single_spectrum(csv_file, value_column, fs, title, output_png, xmax=None):
    df = pd.read_csv(csv_file)
    freqs, mag_db = compute_spectrum_db(df[value_column], fs)

    plt.figure(figsize=(10, 5))
    plt.plot(freqs, mag_db)
    plt.title(title)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [dB, normalized]")
    plt.grid(True)
    plt.xlim(0, xmax if xmax is not None else fs / 2)
    plt.ylim(-120, 5)
    plt.tight_layout()
    plt.savefig(output_png)
    plt.show()


def plot_output_comparison():
    fir = pd.read_csv("interpolated_fir.csv")
    q15 = pd.read_csv("interpolated_q15.csv")

    fir_freqs, fir_mag = compute_spectrum_raw(fir["x"], FS2)
    _, q15_mag = compute_spectrum_raw(q15["q_restored"], FS2)

    ref = np.max(fir_mag)
    fir_db = 20 * np.log10(fir_mag / ref + EPS)
    q15_db = 20 * np.log10(q15_mag / ref + EPS)

    plt.figure(figsize=(10, 5))
    plt.plot(fir_freqs, fir_db, label="Floating-point FIR")
    plt.plot(fir_freqs, q15_db, label="Q15 fixed-point FIR", alpha=0.7)
    plt.title("Output spectrum: floating-point FIR vs Q15 fixed-point FIR")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [dB, normalized to float peak]")
    plt.grid(True)
    plt.legend()
    plt.ylim(-120, 5)
    plt.xlim(0, FS2 / 2)
    plt.tight_layout()
    plt.savefig("spectrum_output_compare.png")
    plt.show()


plot_single_spectrum(
    csv_file="sin.csv",
    value_column="x",
    fs=FS,
    title="Input sine spectrum",
    output_png="spectrum_input.png"
)

plot_single_spectrum(
    csv_file="interpolated_fir.csv",
    value_column="x",
    fs=FS2,
    title="Output spectrum after floating-point FIR interpolation",
    output_png="spectrum_interpolated_fir.png"
)

plot_single_spectrum(
    csv_file="interpolated_q15.csv",
    value_column="q_restored",
    fs=FS2,
    title="Output spectrum after Q15 fixed-point FIR interpolation",
    output_png="spectrum_interpolated_q15.png"
)

plot_output_comparison()