// point_field_detector.hpp - Point field time series detector

#pragma once

#include "idetector.hpp"
#include <sstream>

namespace Detectors {

struct PointFieldDetectorConfig {
    std::string name = "point_probe";
    real x = 0.0, y = 0.0, z = 0.0;  // Position in physical coordinates (meters)
    std::vector<FieldComponent> components = {FieldComponent::Ez};
    std::size_t save_every = 1;
    std::size_t n_steps = 0;
    bool write_float64 = true;
};

struct PointFieldDetector final : public IDetector {
    fs::path det_dir;
    std::size_t NxT{}, NyT{}, NzT{};
    std::size_t i0{}, j0{}, k0{};
    std::vector<FieldComponent> components;
    std::vector<std::ofstream> output_files;
    std::size_t save_every{1};
    std::size_t n_steps{0};
    real dt_sim{0};
    bool write_float64{true};
    real probe_x_physical{}, probe_y_physical{}, probe_z_physical{};
    real pml_offset_x{}, pml_offset_y{}, pml_offset_z{};
    std::string detector_name{"PointFieldDetector"};

    PointFieldDetector() = default;

    ~PointFieldDetector() {
        for (auto& f : output_files) {
            if (f.is_open()) f.close();
        }
    }

    void initialize(const fs::path& out_root, const PointFieldDetectorConfig& config,
                    std::size_t NxT_, std::size_t NyT_, std::size_t NzT_,
                    std::size_t npml,
                    const GridSpacing& grid,
                    real dt,
                    real domain_min_x, real domain_min_y, real domain_min_z)
    {
        det_dir = out_root / config.name;
        NxT = NxT_; NyT = NyT_; NzT = NzT_;
        components = config.components;
        save_every = config.save_every;
        n_steps = config.n_steps;
        write_float64 = config.write_float64;
        dt_sim = dt;

        pml_offset_x = grid.pml_offset_x();
        pml_offset_y = grid.pml_offset_y();
        pml_offset_z = grid.pml_offset_z();

        real phys_x = config.x - domain_min_x;
        real phys_y = config.y - domain_min_y;
        real phys_z = config.z - domain_min_z;

        i0 = grid.physical_to_index_x(phys_x);
        j0 = grid.physical_to_index_y(phys_y);
        k0 = grid.physical_to_index_z(phys_z);

        i0 = std::max(std::size_t(1), std::min(i0, NxT - 2));
        j0 = std::max(std::size_t(1), std::min(j0, NyT - 2));
        k0 = std::max(std::size_t(1), std::min(k0, NzT - 2));

        probe_x_physical = grid.cell_center_physical_x(i0);
        probe_y_physical = grid.cell_center_physical_y(j0);
        probe_z_physical = grid.cell_center_physical_z(k0);

        create_detector_directory(det_dir);

        output_files.resize(components.size());
        for (std::size_t c = 0; c < components.size(); ++c) {
            std::string filename = std::string(field_component_name(components[c])) + "_ts.bin";
            fs::path path = det_dir / filename;
            output_files[c].open(path, std::ios::binary);
            if (!output_files[c]) {
                std::cerr << "[ERR] Failed to open " << path << "\n";
            }
        }

        write_metadata();

        std::ostringstream oss;
        oss << "PointFieldDetector@(" << i0 << "," << j0 << "," << k0 << ")[";
        for (std::size_t c = 0; c < components.size(); ++c) {
            if (c > 0) oss << ",";
            oss << field_component_name(components[c]);
        }
        oss << "]";
        detector_name = oss.str();

        std::cout << "[PointFieldDetector] " << det_dir << "\n";
        std::cout << "  Position: (" << probe_x_physical * 1e9 << ", "
                  << probe_y_physical * 1e9 << ", " << probe_z_physical * 1e9 << ") nm\n";
        std::cout << "  Grid indices: (" << i0 << ", " << j0 << ", " << k0 << ")\n";
        std::cout << "  Components: ";
        for (std::size_t c = 0; c < components.size(); ++c) {
            if (c > 0) std::cout << ", ";
            std::cout << field_component_name(components[c]);
        }
        std::cout << "\n";
    }

    void record_after_E(std::size_t n, real dt,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz) override
    {
        (void)dt;
        if (n % save_every != 0) return;

        std::size_t idx = idx3(i0, j0, k0, NyT, NzT);
        for (std::size_t c = 0; c < components.size(); ++c) {
            if (!output_files[c]) continue;
            real v = get_field_value(components[c], idx, Ex, Ey, Ez, Hx, Hy, Hz);
            write_binary_value(output_files[c], v, write_float64);
        }
    }

    std::string name() const override { return detector_name; }

    void finalize() override {
        for (auto& f : output_files) {
            if (f.is_open()) f.close();
        }
    }

private:
    void write_metadata() {
        fs::path json_path = det_dir / "metadata.json";
        std::ofstream ofs(json_path);
        if (!ofs) return;

        ofs << std::setprecision(15);
        ofs << "{\n";
        ofs << "  \"detector_type\": \"PointFieldDetector\",\n";
        ofs << "  \"NxT\": " << NxT << ",\n";
        ofs << "  \"NyT\": " << NyT << ",\n";
        ofs << "  \"NzT\": " << NzT << ",\n";
        ofs << "  \"i0\": " << i0 << ",\n";
        ofs << "  \"j0\": " << j0 << ",\n";
        ofs << "  \"k0\": " << k0 << ",\n";
        ofs << "  \"probe_x_physical_m\": " << probe_x_physical << ",\n";
        ofs << "  \"probe_y_physical_m\": " << probe_y_physical << ",\n";
        ofs << "  \"probe_z_physical_m\": " << probe_z_physical << ",\n";
        ofs << "  \"probe_x_physical_nm\": " << probe_x_physical * 1e9 << ",\n";
        ofs << "  \"probe_y_physical_nm\": " << probe_y_physical * 1e9 << ",\n";
        ofs << "  \"probe_z_physical_nm\": " << probe_z_physical * 1e9 << ",\n";
        ofs << "  \"pml_offset_x_m\": " << pml_offset_x << ",\n";
        ofs << "  \"pml_offset_y_m\": " << pml_offset_y << ",\n";
        ofs << "  \"pml_offset_z_m\": " << pml_offset_z << ",\n";
        ofs << "  \"save_every\": " << save_every << ",\n";
        ofs << "  \"n_steps\": " << n_steps << ",\n";
        ofs << "  \"dt\": " << std::setprecision(18) << dt_sim << ",\n";
        ofs << "  \"dtype\": \"" << (write_float64 ? "float64" : "float32") << "\",\n";

        ofs << "  \"components\": [";
        for (std::size_t c = 0; c < components.size(); ++c) {
            if (c > 0) ofs << ", ";
            ofs << "\"" << field_component_name(components[c]) << "\"";
        }
        ofs << "],\n";

        ofs << "  \"files\": [";
        for (std::size_t c = 0; c < components.size(); ++c) {
            if (c > 0) ofs << ", ";
            ofs << "\"" << field_component_name(components[c]) << "_ts.bin\"";
        }
        ofs << "],\n";

        ofs << "  \"note\": \"time series of field components at single point\"\n";
        ofs << "}\n";
    }
};

inline std::unique_ptr<PointFieldDetector> make_point_field_detector(
    const fs::path& out_root,
    const PointFieldDetectorConfig& config,
    std::size_t NxT, std::size_t NyT, std::size_t NzT,
    std::size_t npml,
    const GridSpacing& grid,
    real dt,
    real domain_min_x = 0.0, real domain_min_y = 0.0, real domain_min_z = 0.0)
{
    auto det = std::make_unique<PointFieldDetector>();
    det->initialize(out_root, config, NxT, NyT, NzT, npml, grid, dt,
                    domain_min_x, domain_min_y, domain_min_z);
    return det;
}

inline std::unique_ptr<PointFieldDetector> make_simple_probe(
    const fs::path& out_root,
    const std::string& name,
    real x, real y, real z,
    FieldComponent component,
    std::size_t NxT, std::size_t NyT, std::size_t NzT,
    std::size_t npml,
    const GridSpacing& grid,
    real dt,
    std::size_t save_every = 1,
    std::size_t n_steps = 0,
    real domain_min_x = 0.0, real domain_min_y = 0.0, real domain_min_z = 0.0)
{
    PointFieldDetectorConfig config;
    config.name = name;
    config.x = x;
    config.y = y;
    config.z = z;
    config.components = {component};
    config.save_every = save_every;
    config.n_steps = n_steps;

    return make_point_field_detector(out_root, config, NxT, NyT, NzT, npml, grid, dt,
                                     domain_min_x, domain_min_y, domain_min_z);
}

} // namespace Detectors
