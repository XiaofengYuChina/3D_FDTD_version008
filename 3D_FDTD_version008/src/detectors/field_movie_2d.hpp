// field_movie_2d.hpp - 2D field movie detector

#pragma once

#include "idetector.hpp"
#include <sstream>
#include <cstdio>

namespace Detectors {

struct FieldMovie2DConfig {
    std::string name = "field_movie";
    FieldComponent component = FieldComponent::Ez;
    int slice_plane = 0;             // 0=XY, 1=XZ, 2=YZ
    real slice_position = 0.0;       // Position in physical coords (meters)
    std::size_t save_every = 1;
    std::size_t n_steps = 0;
    bool write_float64 = true;
    std::string frame_pattern = "frame_%04d.raw";
};

struct FieldMovie2D final : public IDetector {
    fs::path det_dir;
    std::size_t NxT{}, NyT{}, NzT{};
    std::size_t slice_index{};
    int slice_plane{};
    FieldComponent component{FieldComponent::Ez};
    std::size_t save_every{1};
    std::size_t n_steps{0};
    std::size_t frame_id{0};
    bool write_float64{true};
    std::string frame_pattern{"frame_%04d.raw"};
    real pml_offset_x{}, pml_offset_y{}, pml_offset_z{};
    real slice_physical_coord{};
    std::size_t Nx_phys{}, Ny_phys{}, Nz_phys{}, npml{};
    real Lx_phys{}, Ly_phys{}, Lz_phys{};
    std::string detector_name{"FieldMovie2D"};

    FieldMovie2D() = default;

    const char* slice_plane_name() const {
        switch (slice_plane) {
        case 0: return "XY";
        case 1: return "XZ";
        case 2: return "YZ";
        }
        return "Unknown";
    }

    void initialize(const fs::path& out_root, const FieldMovie2DConfig& config,
                    std::size_t NxT_, std::size_t NyT_, std::size_t NzT_,
                    std::size_t Nx_phys_, std::size_t Ny_phys_, std::size_t Nz_phys_,
                    std::size_t npml_,
                    const GridSpacing& grid,
                    real Lx_phys_, real Ly_phys_, real Lz_phys_,
                    real domain_min_x, real domain_min_y, real domain_min_z)
    {
        det_dir = out_root / config.name;
        NxT = NxT_; NyT = NyT_; NzT = NzT_;
        Nx_phys = Nx_phys_; Ny_phys = Ny_phys_; Nz_phys = Nz_phys_;
        npml = npml_;
        slice_plane = config.slice_plane;
        component = config.component;
        save_every = config.save_every;
        n_steps = config.n_steps;
        write_float64 = config.write_float64;
        frame_pattern = config.frame_pattern;
        Lx_phys = Lx_phys_; Ly_phys = Ly_phys_; Lz_phys = Lz_phys_;

        pml_offset_x = grid.pml_offset_x();
        pml_offset_y = grid.pml_offset_y();
        pml_offset_z = grid.pml_offset_z();

        real slice_phys = config.slice_position;
        switch (slice_plane) {
        case 0:
            slice_phys -= domain_min_z;
            slice_index = grid.physical_to_index_z(slice_phys);
            slice_index = std::max(std::size_t(1), std::min(slice_index, NzT - 2));
            slice_physical_coord = grid.cell_center_physical_z(slice_index);
            break;
        case 1:
            slice_phys -= domain_min_y;
            slice_index = grid.physical_to_index_y(slice_phys);
            slice_index = std::max(std::size_t(1), std::min(slice_index, NyT - 2));
            slice_physical_coord = grid.cell_center_physical_y(slice_index);
            break;
        case 2:
            slice_phys -= domain_min_x;
            slice_index = grid.physical_to_index_x(slice_phys);
            slice_index = std::max(std::size_t(1), std::min(slice_index, NxT - 2));
            slice_physical_coord = grid.cell_center_physical_x(slice_index);
            break;
        }

        create_detector_directory(det_dir);
        write_metadata();

        std::ostringstream oss;
        oss << "FieldMovie2D[" << field_component_name(component) << ","
            << slice_plane_name() << "@" << slice_index << "]";
        detector_name = oss.str();

        std::cout << "[FieldMovie2D] " << det_dir << " (" << field_component_name(component)
                  << " on " << slice_plane_name() << " plane at index " << slice_index << ")\n";
    }

    void record_after_E(std::size_t n, real dt,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz) override
    {
        (void)dt;
        if (n % save_every != 0) return;

        char fname[128];
        std::snprintf(fname, sizeof(fname), frame_pattern.c_str(), static_cast<int>(frame_id++));
        fs::path out_path = det_dir / fname;

        std::ofstream ofs(out_path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[ERR] Failed to open " << out_path << "\n";
            return;
        }

        switch (slice_plane) {
        case 0:
            for (std::size_t i = 0; i < NxT; ++i) {
                for (std::size_t j = 0; j < NyT; ++j) {
                    std::size_t idx = idx3(i, j, slice_index, NyT, NzT);
                    real v = get_field_value(component, idx, Ex, Ey, Ez, Hx, Hy, Hz);
                    write_binary_value(ofs, v, write_float64);
                }
            }
            break;
        case 1:
            for (std::size_t i = 0; i < NxT; ++i) {
                for (std::size_t k = 0; k < NzT; ++k) {
                    std::size_t idx = idx3(i, slice_index, k, NyT, NzT);
                    real v = get_field_value(component, idx, Ex, Ey, Ez, Hx, Hy, Hz);
                    write_binary_value(ofs, v, write_float64);
                }
            }
            break;
        case 2:
            for (std::size_t j = 0; j < NyT; ++j) {
                for (std::size_t k = 0; k < NzT; ++k) {
                    std::size_t idx = idx3(slice_index, j, k, NyT, NzT);
                    real v = get_field_value(component, idx, Ex, Ey, Ez, Hx, Hy, Hz);
                    write_binary_value(ofs, v, write_float64);
                }
            }
            break;
        }
    }

    std::string name() const override { return detector_name; }

private:
    void write_metadata() {
        fs::path json_path = det_dir / "metadata.json";
        std::ofstream ofs(json_path);
        if (!ofs) return;

        ofs << std::setprecision(15);
        ofs << "{\n";
        ofs << "  \"detector_type\": \"FieldMovie2D\",\n";
        ofs << "  \"field_component\": \"" << field_component_name(component) << "\",\n";
        ofs << "  \"slice_plane\": \"" << slice_plane_name() << "\",\n";
        ofs << "  \"slice_index\": " << slice_index << ",\n";
        ofs << "  \"slice_physical_m\": " << slice_physical_coord << ",\n";
        ofs << "  \"slice_physical_nm\": " << slice_physical_coord * 1e9 << ",\n";
        ofs << "  \"NxT\": " << NxT << ",\n";
        ofs << "  \"NyT\": " << NyT << ",\n";
        ofs << "  \"NzT\": " << NzT << ",\n";
        ofs << "  \"Nx_phys\": " << Nx_phys << ",\n";
        ofs << "  \"Ny_phys\": " << Ny_phys << ",\n";
        ofs << "  \"Nz_phys\": " << Nz_phys << ",\n";
        ofs << "  \"npml\": " << npml << ",\n";
        ofs << "  \"save_every\": " << save_every << ",\n";
        ofs << "  \"n_steps\": " << n_steps << ",\n";
        ofs << "  \"frame_pattern\": \"" << frame_pattern << "\",\n";
        ofs << "  \"pml_offset_x_m\": " << pml_offset_x << ",\n";
        ofs << "  \"pml_offset_y_m\": " << pml_offset_y << ",\n";
        ofs << "  \"pml_offset_z_m\": " << pml_offset_z << ",\n";
        ofs << "  \"Lx_phys_m\": " << Lx_phys << ",\n";
        ofs << "  \"Ly_phys_m\": " << Ly_phys << ",\n";
        ofs << "  \"Lz_phys_m\": " << Lz_phys << ",\n";
        ofs << "  \"dtype\": \"" << (write_float64 ? "float64" : "float32") << "\",\n";

        switch (slice_plane) {
        case 0:
            ofs << "  \"slice_dim1\": " << NxT << ",\n";
            ofs << "  \"slice_dim2\": " << NyT << ",\n";
            ofs << "  \"dim1_label\": \"x\",\n";
            ofs << "  \"dim2_label\": \"y\",\n";
            break;
        case 1:
            ofs << "  \"slice_dim1\": " << NxT << ",\n";
            ofs << "  \"slice_dim2\": " << NzT << ",\n";
            ofs << "  \"dim1_label\": \"x\",\n";
            ofs << "  \"dim2_label\": \"z\",\n";
            break;
        case 2:
            ofs << "  \"slice_dim1\": " << NyT << ",\n";
            ofs << "  \"slice_dim2\": " << NzT << ",\n";
            ofs << "  \"dim1_label\": \"y\",\n";
            ofs << "  \"dim2_label\": \"z\",\n";
            break;
        }

        ofs << "  \"note\": \"row-major storage, dim1 varies fastest\"\n";
        ofs << "}\n";
    }
};

inline std::unique_ptr<FieldMovie2D> make_field_movie_2d(
    const fs::path& out_root,
    const FieldMovie2DConfig& config,
    std::size_t NxT, std::size_t NyT, std::size_t NzT,
    std::size_t Nx_phys, std::size_t Ny_phys, std::size_t Nz_phys,
    std::size_t npml,
    const GridSpacing& grid,
    real Lx_phys, real Ly_phys, real Lz_phys,
    real domain_min_x = 0.0, real domain_min_y = 0.0, real domain_min_z = 0.0)
{
    auto det = std::make_unique<FieldMovie2D>();
    det->initialize(out_root, config, NxT, NyT, NzT,
                    Nx_phys, Ny_phys, Nz_phys, npml,
                    grid, Lx_phys, Ly_phys, Lz_phys,
                    domain_min_x, domain_min_y, domain_min_z);
    return det;
}

} // namespace Detectors
