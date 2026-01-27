// idetector.hpp - Detector interface and common definitions

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>
#include <iomanip>

#include "../global_function.hpp"

namespace fs = std::filesystem;

namespace Detectors {

enum class FieldComponent {
    Ex, Ey, Ez,
    Hx, Hy, Hz,
    E_magnitude,
    H_magnitude,
    Sx, Sy, Sz,
    S_magnitude,
    Energy_E,
    Energy_H,
    Energy_total
};

inline const char* field_component_name(FieldComponent fc) {
    switch (fc) {
    case FieldComponent::Ex: return "Ex";
    case FieldComponent::Ey: return "Ey";
    case FieldComponent::Ez: return "Ez";
    case FieldComponent::Hx: return "Hx";
    case FieldComponent::Hy: return "Hy";
    case FieldComponent::Hz: return "Hz";
    case FieldComponent::E_magnitude: return "E_mag";
    case FieldComponent::H_magnitude: return "H_mag";
    case FieldComponent::Sx: return "Sx";
    case FieldComponent::Sy: return "Sy";
    case FieldComponent::Sz: return "Sz";
    case FieldComponent::S_magnitude: return "S_mag";
    case FieldComponent::Energy_E: return "U_E";
    case FieldComponent::Energy_H: return "U_H";
    case FieldComponent::Energy_total: return "U_total";
    }
    return "Unknown";
}

struct IDetector {
    virtual ~IDetector() = default;

    virtual void record_after_E(std::size_t n, real dt,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz) = 0;

    virtual void record_after_H(std::size_t n, real dt,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz) {
        (void)n; (void)dt;
        (void)Ex; (void)Ey; (void)Ez;
        (void)Hx; (void)Hy; (void)Hz;
    }

    virtual std::string name() const { return "IDetector"; }
    virtual void finalize() {}
};

struct DetectorConfig {
    std::string name = "detector";
    std::size_t save_every = 1;
    std::size_t n_steps = 0;
    bool write_float64 = true;
    std::string frame_pattern = "frame_%04d.raw";
};

inline real get_field_value(FieldComponent component, std::size_t idx,
    const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
    const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz)
{
    switch (component) {
    case FieldComponent::Ex: return Ex[idx];
    case FieldComponent::Ey: return Ey[idx];
    case FieldComponent::Ez: return Ez[idx];
    case FieldComponent::Hx: return Hx[idx];
    case FieldComponent::Hy: return Hy[idx];
    case FieldComponent::Hz: return Hz[idx];
    case FieldComponent::E_magnitude:
        return std::sqrt(Ex[idx]*Ex[idx] + Ey[idx]*Ey[idx] + Ez[idx]*Ez[idx]);
    case FieldComponent::H_magnitude:
        return std::sqrt(Hx[idx]*Hx[idx] + Hy[idx]*Hy[idx] + Hz[idx]*Hz[idx]);
    case FieldComponent::Sx:
        return Ey[idx]*Hz[idx] - Ez[idx]*Hy[idx];
    case FieldComponent::Sy:
        return Ez[idx]*Hx[idx] - Ex[idx]*Hz[idx];
    case FieldComponent::Sz:
        return Ex[idx]*Hy[idx] - Ey[idx]*Hx[idx];
    case FieldComponent::S_magnitude: {
        real Sx = Ey[idx]*Hz[idx] - Ez[idx]*Hy[idx];
        real Sy = Ez[idx]*Hx[idx] - Ex[idx]*Hz[idx];
        real Sz = Ex[idx]*Hy[idx] - Ey[idx]*Hx[idx];
        return std::sqrt(Sx*Sx + Sy*Sy + Sz*Sz);
    }
    default:
        return 0.0;
    }
}

inline void write_binary_value(std::ofstream& ofs, real value, bool write_float64) {
    if (write_float64) {
        double v = static_cast<double>(value);
        ofs.write(reinterpret_cast<const char*>(&v), sizeof(double));
    } else {
        float v = static_cast<float>(value);
        ofs.write(reinterpret_cast<const char*>(&v), sizeof(float));
    }
}

inline bool create_detector_directory(const fs::path& dir) {
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        std::cerr << "[ERR] Failed to create directory " << dir << ": " << ec.message() << "\n";
        return false;
    }
    return true;
}

} // namespace Detectors
