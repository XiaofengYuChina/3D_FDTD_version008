// idetector.hpp - Detector interface and common definitions
//
// This file defines the base interface for all electromagnetic field detectors
// in the FDTD simulation.

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

// ==================== Field component selection ====================
enum class FieldComponent {
    Ex, Ey, Ez,    // Electric field components
    Hx, Hy, Hz,    // Magnetic field components
    E_magnitude,   // |E| = sqrt(Ex^2 + Ey^2 + Ez^2)
    H_magnitude,   // |H| = sqrt(Hx^2 + Hy^2 + Hz^2)
    Sx, Sy, Sz,    // Poynting vector components
    S_magnitude,   // |S| magnitude
    Energy_E,      // Electric energy density
    Energy_H,      // Magnetic energy density
    Energy_total   // Total energy density
};

// Convert FieldComponent to string
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

// ==================== Detector interface ====================
struct IDetector {
    virtual ~IDetector() = default;

    // Called after E-field update
    virtual void record_after_E(std::size_t n, real dt,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz) = 0;

    // Called after H-field update (optional, default does nothing)
    virtual void record_after_H(std::size_t n, real dt,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz) {
        (void)n; (void)dt;
        (void)Ex; (void)Ey; (void)Ez;
        (void)Hx; (void)Hy; (void)Hz;
    }

    // Get detector name for logging
    virtual std::string name() const { return "IDetector"; }

    // Finalize (called at end of simulation)
    virtual void finalize() {}
};

// ==================== Common detector configuration ====================
struct DetectorConfig {
    std::string name = "detector";           // Detector name (used for output directory)
    std::size_t save_every = 1;              // Save every N steps
    std::size_t n_steps = 0;                 // Total simulation steps (for metadata)
    bool write_float64 = true;               // true = double, false = float
    std::string frame_pattern = "frame_%04d.raw";  // Output file pattern
};

// ==================== Utility functions ====================

// Extract field value at a grid point
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

// Write binary value (float32 or float64)
inline void write_binary_value(std::ofstream& ofs, real value, bool write_float64) {
    if (write_float64) {
        double v = static_cast<double>(value);
        ofs.write(reinterpret_cast<const char*>(&v), sizeof(double));
    } else {
        float v = static_cast<float>(value);
        ofs.write(reinterpret_cast<const char*>(&v), sizeof(float));
    }
}

// Create output directory safely
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
