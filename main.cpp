#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include <filesystem>

const double PI = 3.14159265358979323846;

const double C_freq = 0.5;
const double FS = 100.0;
const double DURATION =  10.0;

const int FIR_NUM_COEFFS = 61;
const double FS2 = 2 * FS;
const int Q15_SCALE = 1 << 15;
const int COEFF_FRAC_BITS = 15;
const int COEFF_SCALE = 1 << COEFF_FRAC_BITS;

const std::string RESULTS_DIR = "results";

std::string result_path(const std::string& filename) {
    return RESULTS_DIR + "/" + filename;
}

// Пункт 1. Генерация сигнала
std::vector<double> make_sin(double freq_hz, double sample_rate_hz) {
    int count = static_cast<int>(DURATION*sample_rate_hz);
    std::vector<double> signal(count);

    for (int n = 0; n < count; n++) {
        double t = n / sample_rate_hz;
        signal[n] = sin(2 * PI * freq_hz * t);
    }

    return signal;
}

void save_csv(const std::vector<double>& signal, 
              double sample_rate_hz,
              const std::string& filename){
                
    std::ofstream file(filename);
    file << "t,x\n";
    
    for (size_t n = 0; n < signal.size(); n++){
        double t = n / sample_rate_hz;
        file << t << "," << signal[n] << "\n";
    }

}

// Пункт 2. Квантование     
std::vector<int16_t> quantize_q15(const std::vector<double>& signal){
    std::vector<int16_t> quantized_signal(signal.size());

    for (size_t i = 0; i < signal.size(); i++){
        double x = std::clamp(signal[i], -1.0, 1.0);
        int32_t q = static_cast<int32_t>(std::round(x * Q15_SCALE));
        q = std::clamp(q, static_cast<int32_t>(INT16_MIN), static_cast<int32_t>(INT16_MAX));
        quantized_signal[i] = static_cast<int16_t>(q);
    }

    return quantized_signal;
}

void save_quantized_csv(const std::vector<double>& original,
                        const std::vector<int16_t>& quantized,
                        double sample_rate_hz,
                        const std::string& filename) {
    std::ofstream file(filename);
    file << "t,original,q_int,q_restored,error\n";

    for (size_t n = 0; n < original.size(); n++) {
        double t = n / sample_rate_hz;
        double restored = quantized[n] / static_cast<double>(Q15_SCALE);
        double error = original[n] - restored;
        
        file << t << ","
             << original[n] << ","
             << quantized[n] << ","
             << restored << ","
             << error << "\n";
    }
}

void analyze_quantization_error(const std::vector<double>& original,
                                const std::vector<int16_t>& quantized){
    
    double max_abs_error = 0.0;
    double mean_error = 0.0;
    double signal_power = 0.0;
    double noise_power = 0.0;

    for (size_t i = 0; i < original.size(); i++){
        double restored = quantized[i] / static_cast<double>(Q15_SCALE);
        double error = original[i] - restored;

        mean_error += error;
        signal_power += original[i] * original[i];
        noise_power += error * error;

        if (std::abs(error) > max_abs_error){
            max_abs_error = std::abs(error);
        }
    }

    mean_error /= original.size();
    signal_power /= original.size();
    noise_power /= original.size();

    double rms_error = std::sqrt(noise_power);
    double snr_db = 10 * std::log10(signal_power / noise_power);

    std::ofstream report(result_path("quantization_metrics.txt"));

    report << "Quantization error analysis\n";
    report << "Max absolute error: " << max_abs_error << "\n";
    report << "Mean error: " << mean_error << "\n";
    report << "RMS error: " << rms_error << "\n";
    report << "SNR (dB): " << snr_db << "\n";

    report.close();

    std::cout << "Quantization metrics saved to quantization_metrics.txt\n";
}

// Пункт 3. Интреполяция плавающей точки 
double sinc(double x) {
    if (std::abs(x) < 1e-12) {
        return 1.0;
    }

    return std::sin(PI * x) / (PI * x);
}

std::vector<double> design_fir_interpolator_x2(size_t num_coeffs = FIR_NUM_COEFFS){
    if (num_coeffs % 2 == 0){
        num_coeffs += 1;
    }

    std::vector<double> h(num_coeffs);
    const double L = 2;
    const double cutoff = 0.25;
    const int center = static_cast<int>(num_coeffs / 2);
    
    for (size_t n = 0; n < num_coeffs; n++) {
        int k = static_cast<int>(n) - center;

        double ideal = L * 2.0 * cutoff * sinc(2.0 * cutoff * k);

        double window = 0.54 - 0.46 * std::cos(
            2.0 * PI * static_cast<double>(n) / static_cast<double>(num_coeffs - 1)
        );

        h[n] = ideal * window;
    }

    return h;
}

std::vector<double> interpolate_fir_x2(const std::vector<double>& signal) {
    if (signal.empty()) {
        return {};
    }

    const size_t L = 2;

    std::vector<double> upsampled(signal.size() * L, 0.0);

    for (size_t i = 0; i < signal.size(); i++) {
        upsampled[i * L] = signal[i];
    }

    std::vector<double> h = design_fir_interpolator_x2(FIR_NUM_COEFFS);
    int delay = static_cast<int>(h.size() / 2);

    std::vector<double> result(upsampled.size(), 0.0);

    for (size_t n = 0; n < upsampled.size(); n++) {
        double acc = 0.0;

        for (size_t k = 0; k < h.size(); k++) {
            int idx = static_cast<int>(n) + static_cast<int>(k) - delay;

            if (idx >= 0 && idx < static_cast<int>(upsampled.size())) {
                acc += h[k] * upsampled[idx];
            }
        }

        result[n] = acc;
    }

    return result;
}

// Пункт 5. Интерполяция квантованного сигнала  фикс. точки
std::vector<int32_t> quantize_fir_coefficients_q15(const std::vector<double>& h) {
    std::vector<int32_t> h_q15(h.size());

    for (size_t i = 0; i < h.size(); i++) {
        h_q15[i] = static_cast<int32_t>(std::round(h[i] * COEFF_SCALE));
    }

    return h_q15;
}

int16_t saturate_to_int16(int64_t value) {
    if (value > INT16_MAX) {
        return INT16_MAX;
    }

    if (value < INT16_MIN) {
        return INT16_MIN;
    }

    return static_cast<int16_t>(value);
}

std::vector<int16_t> interpolate_fir_x2_q15(const std::vector<int16_t>& signal) {
    if (signal.empty()) {
        return {};
    }

    const size_t L = 2;

    std::vector<int16_t> upsampled(signal.size() * L, 0);

    for (size_t i = 0; i < signal.size(); i++) {
        upsampled[i * L] = signal[i];
    }

    std::vector<double> h_float = design_fir_interpolator_x2(FIR_NUM_COEFFS);
    std::vector<int32_t> h_q15 = quantize_fir_coefficients_q15(h_float);

    int delay = static_cast<int>(h_q15.size() / 2);

    std::vector<int16_t> result(upsampled.size(), 0);

    for (size_t n = 0; n < upsampled.size(); n++) {
        int64_t acc = 0;

        for (size_t k = 0; k < h_q15.size(); k++) {
            int idx = static_cast<int>(n) + static_cast<int>(k) - delay;

            if (idx >= 0 && idx < static_cast<int>(upsampled.size())) {
                acc += static_cast<int64_t>(h_q15[k]) *
                       static_cast<int64_t>(upsampled[idx]);
            }
        }

        int64_t rounded;

        if (acc >= 0) {
            rounded = (acc + COEFF_SCALE / 2) / COEFF_SCALE;
        } else {
            rounded = (acc - COEFF_SCALE / 2) / COEFF_SCALE;
        }

        result[n] = saturate_to_int16(rounded);
    }

    return result;
}

void save_q15_csv(const std::vector<int16_t>& signal,
                  double sample_rate_hz,
                  const std::string& filename) {
    std::ofstream file(filename);
    file << "t,q_int,q_restored\n";

    for (size_t n = 0; n < signal.size(); n++) {
        double t = n / sample_rate_hz;
        double restored = signal[n] / static_cast<double>(Q15_SCALE);

        file << t << ","
             << signal[n] << ","
             << restored << "\n";
    }
}

// Пункт 4. Анализ качества интреполяции в зависимости от частоты
struct InterpolationMetrics {
    double max_abs_error;
    double mean_error;
    double rmse;
    double snr_db;
};

InterpolationMetrics compare_signals(const std::vector<double>& reference,
                                      const std::vector<double>& test,
                                      size_t skip_edges = 0) {
    size_t n = std::min(reference.size(), test.size());

    if (2 * skip_edges >= n) {
        skip_edges = 0; // не учитываем крайние точки, так как FIR-фильтр дает искажения
    }

    double max_abs_error = 0.0;
    double mean_error = 0.0;
    double signal_power = 0.0;
    double noise_power = 0.0;

    size_t count = 0;

    for (size_t i = skip_edges; i < n - skip_edges; i++) {
        double error = reference[i] - test[i];

        mean_error += error;
        signal_power += reference[i] * reference[i];
        noise_power += error * error;

        max_abs_error = std::max(max_abs_error, std::abs(error));
        count++;
    }

    mean_error /= count;
    signal_power /= count;
    noise_power /= count;

    InterpolationMetrics metrics;
    metrics.max_abs_error = max_abs_error;
    metrics.mean_error = mean_error;
    metrics.rmse = std::sqrt(noise_power);

    if (noise_power > 0.0) {
        metrics.snr_db = 10.0 * std::log10(signal_power / noise_power);
    } else {
        metrics.snr_db = 0.0;
    }

    return metrics;
}

InterpolationMetrics compare_q15_with_reference(const std::vector<double>& reference,
                                                 const std::vector<int16_t>& test_q15,
                                                 size_t skip_edges = 0) {
    size_t n = std::min(reference.size(), test_q15.size());

    if (2 * skip_edges >= n) {
        skip_edges = 0;
    }

    double max_abs_error = 0.0;
    double mean_error = 0.0;
    double signal_power = 0.0;
    double noise_power = 0.0;

    size_t count = 0;

    for (size_t i = skip_edges; i < n - skip_edges; i++) {
        double restored = test_q15[i] / static_cast<double>(Q15_SCALE);
        double error = reference[i] - restored;

        mean_error += error;
        signal_power += reference[i] * reference[i];
        noise_power += error * error;

        max_abs_error = std::max(max_abs_error, std::abs(error));
        count++;
    }

    mean_error /= count;
    signal_power /= count;
    noise_power /= count;

    InterpolationMetrics metrics;
    metrics.max_abs_error = max_abs_error;
    metrics.mean_error = mean_error;
    metrics.rmse = std::sqrt(noise_power);

    if (noise_power > 0.0) {
        metrics.snr_db = 10.0 * std::log10(signal_power / noise_power);
    } else {
        metrics.snr_db = 0.0;
    }

    return metrics;
}

void analyze_interpolation_quality_by_frequency() {
    std::vector<double> frequencies = {
        0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0
    };

    std::ofstream file(result_path("interpolation_quality.csv"));

    file << "freq_hz,"
         << "float_rmse,"
         << "float_max_abs_error,"
         << "float_mean_error,"
         << "float_snr_db,"
         << "q15_rmse,"
         << "q15_max_abs_error,"
         << "q15_mean_error,"
         << "q15_snr_db\n";

    size_t skip_edges = (FIR_NUM_COEFFS - 1) / 2;

    for (double freq : frequencies) {
        auto original = make_sin(freq, FS);
        auto ideal = make_sin(freq, FS2);

        auto interpolated_float = interpolate_fir_x2(original);
        auto float_metrics = compare_signals(ideal, interpolated_float, skip_edges);

        auto quantized = quantize_q15(original);
        auto interpolated_q15 = interpolate_fir_x2_q15(quantized);
        auto q15_metrics = compare_q15_with_reference(ideal, interpolated_q15, skip_edges);

        file << freq << ","
             << float_metrics.rmse << ","
             << float_metrics.max_abs_error << ","
             << float_metrics.mean_error << ","
             << float_metrics.snr_db << ","
             << q15_metrics.rmse << ","
             << q15_metrics.max_abs_error << ","
             << q15_metrics.mean_error << ","
             << q15_metrics.snr_db << "\n";
    }

    file.close();

    std::cout << "Interpolation quality saved to interpolation_quality.csv\n";
}

int main()
{
    std::filesystem::create_directories(RESULTS_DIR);

    // Пункт 1. Генерация сигнала
    auto signal = make_sin(C_freq, FS);
    save_csv(signal, FS, result_path("sin.csv"));
    std::cout << "saved sin.csv\n";

    // Пункт 2. Квантование синуса в Q15 и сохранение отчета об ошибках
    auto q_signal = quantize_q15(signal);
    save_quantized_csv(signal, q_signal, FS, result_path("quantized.csv"));
    std::cout << "saved quantized.csv\n";

    analyze_quantization_error(signal, q_signal);

    // Пункт 3. Интерполяция при помощи FIR-фильтра
    auto interpolated_fir = interpolate_fir_x2(signal);
    save_csv(interpolated_fir, FS2, result_path("interpolated_fir.csv"));
    std::cout << "saved interpolated_fir.csv\n";  

    // Пункт 4. Анализ качества интреполяции при разных частотах
    analyze_interpolation_quality_by_frequency();

    // Пункт 5. Интерполяция при помощи FIR-фильтра для фикс. точки
    auto interpolated_q15 = interpolate_fir_x2_q15(q_signal);
    save_q15_csv(interpolated_q15, FS2, result_path("interpolated_q15.csv"));
    std::cout << "saved interpolated_q15.csv\n";
    
    return 0;
}



