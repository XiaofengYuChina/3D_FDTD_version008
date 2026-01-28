// tls_materials.hpp - Two-Level System Material Library
//
// Usage: Define materials in LIBRARY, reference by name in STRUCTURES:
//   {"cylinder", 3.3, {...}, "gain_1500nm"}

#pragma once

#include <string>
#include <vector>

namespace TLSMaterials {

struct TLSMaterial {
    std::string name;
    double lambda0;             // Transition wavelength (m)
    double gamma;               // Polarization damping rate (1/s)
    double tau;                 // Upper level lifetime (s)
    double N0;                  // Dipole density (atoms/m^3)
    double inversion_fraction;  // Initial inversion: 0.5=equilibrium, >0.5=gain, 1.0=full
};

// Format: {name, lambda0, gamma, tau, N0, inversion_fraction}
inline const std::vector<TLSMaterial> LIBRARY = {
    {"gain_1500nm", 1500e-9, 7e12, 1e-10, 1e23, 0.5},
};

inline const TLSMaterial* find(const std::string& name) {
    if (name.empty()) return nullptr;
    for (const auto& mat : LIBRARY)
        if (mat.name == name) return &mat;
    return nullptr;
}

inline bool exists(const std::string& name) {
    return find(name) != nullptr;
}

inline std::vector<std::string> get_all_names() {
    std::vector<std::string> names;
    for (const auto& mat : LIBRARY)
        names.push_back(mat.name);
    return names;
}

} // namespace TLSMaterials
