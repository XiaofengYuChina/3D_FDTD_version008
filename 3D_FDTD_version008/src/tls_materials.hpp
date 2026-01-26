// tls_materials.hpp — Two-Level System Material Library
//
// Define all your TLS (gain medium) materials here.
// Each material has a unique name that can be referenced in STRUCTURES.
//
// Usage:
//   1. Define materials in TLS_MATERIAL_LIBRARY below
//   2. Reference by name in user_config.hpp STRUCTURES:
//      {"cylinder", 3.3, {...}, "gain_1500nm"}

#pragma once

#include <string>
#include <vector>

namespace TLSMaterials {

// ============================================================================
//                         TLS Material Definition
// ============================================================================

struct TLSMaterial {
    std::string name;           // Unique name to reference this material
    double lambda0;             // Transition wavelength (m)
    double gamma;               // Polarization damping rate (1/s)
    double tau;                 // Upper level lifetime (s)
    double N0;                  // Dipole density (atoms/m³)
    double inversion_fraction;  // Initial population inversion (0 to 1)
                                //   0.5 = thermal equilibrium
                                //   > 0.5 = population inversion (gain)
                                //   1.0 = full inversion
};

// ============================================================================
//                         TLS Material Library
// ============================================================================
//
// Define your TLS materials here. Each material needs:
//   - name: unique identifier string
//   - lambda0: transition wavelength in meters
//   - gamma: polarization damping rate in 1/s
//   - tau: upper level lifetime in seconds
//   - N0: dipole density in atoms/m³
//   - inversion_fraction: initial population inversion (0 to 1)
//
// Format: {name, lambda0, gamma, tau, N0, inversion_fraction}

inline const std::vector<TLSMaterial> LIBRARY = {

    // ===== 1500 nm Gain Materials =====
    {"gain_1500nm",
        1500e-9,    // lambda0: 1500 nm
        7e12,       // gamma: 7 THz damping
        1e-12,      // tau: 1 ps lifetime
        1e25,       // N0: 10^25 atoms/m³
        1.0         // full inversion
    },

    // ===== Example: 1300 nm Gain Material =====
    // {"gain_1300nm",
    //     1300e-9,    // lambda0: 1300 nm
    //     5e12,       // gamma: 5 THz damping
    //     2e-12,      // tau: 2 ps lifetime
    //     5e24,       // N0: 5×10^24 atoms/m³
    //     0.8         // 80% inversion
    // },

    // ===== Example: Weak Gain Material =====
    // {"weak_gain",
    //     1550e-9,    // lambda0: 1550 nm
    //     1e13,       // gamma: 10 THz damping
    //     5e-12,      // tau: 5 ps lifetime
    //     1e24,       // N0: 10^24 atoms/m³
    //     0.6         // 60% inversion
    // },

    // ===== Example: Strong Gain Material =====
    // {"strong_gain",
    //     1500e-9,    // lambda0: 1500 nm
    //     5e12,       // gamma: 5 THz damping
    //     0.5e-12,    // tau: 0.5 ps lifetime (short = strong coupling)
    //     2e25,       // N0: 2×10^25 atoms/m³
    //     1.0         // full inversion
    // },

};

// ============================================================================
//                         Helper Functions
// ============================================================================

// Find TLS material by name (returns nullptr if not found)
inline const TLSMaterial* find(const std::string& name) {
    if (name.empty()) return nullptr;
    for (const auto& mat : LIBRARY) {
        if (mat.name == name) return &mat;
    }
    return nullptr;
}

// Check if a material exists
inline bool exists(const std::string& name) {
    return find(name) != nullptr;
}

// Get all material names (for listing/debugging)
inline std::vector<std::string> get_all_names() {
    std::vector<std::string> names;
    for (const auto& mat : LIBRARY) {
        names.push_back(mat.name);
    }
    return names;
}

} // namespace TLSMaterials
