#include "cov/mo_diagram.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace {

cov::MolecularOrbital mo(double energy, float occ, const char* sym) {
    cov::MolecularOrbital orbital;
    orbital.energy_hartree = energy;
    orbital.occupation = occ;
    orbital.symmetry = sym;
    return orbital;
}

} // namespace

int main() {
    cov::Wavefunction wf;
    wf.atoms.resize(10);

    // Cp- Cartesian reference values from the uploaded Gaussian output,
    // orbitals 1..27. The first five are deep core and should be hidden by the
    // default valence diagram; 6..18 are occupied valence; 19..27 are frontier
    // virtuals. This fixture contains both tolerance-level degeneracies and
    // extremely close but symmetry-incompatible levels.
    wf.orbitals = {
        mo(-10.056478,2,"A1'"), mo(-10.056121,2,"E1'"), mo(-10.056120,2,"E1'"),
        mo(-10.055483,2,"E2'"), mo(-10.055482,2,"E2'"),
        mo(-0.741616,2,"A1'"), mo(-0.574601,2,"E1'"), mo(-0.574600,2,"E1'"),
        mo(-0.412893,2,"E2'"), mo(-0.412892,2,"E2'"), mo(-0.399340,2,"A1'"),
        mo(-0.254871,2,"E1'"), mo(-0.254868,2,"E1'"),
        mo(-0.237168,2,"E2'"), mo(-0.237168,2,"E2'"), mo(-0.236093,2,"A2\""),
        mo(-0.068048,2,"E1\""), mo(-0.068048,2,"E1\""),
        mo(0.133376,0,"E1'"), mo(0.133376,0,"E1'"), mo(0.137838,0,"A1'"),
        mo(0.138589,0,"E2'"), mo(0.138590,0,"E2'"), mo(0.224213,0,"E2'"),
        mo(0.224222,0,"E2'"), mo(0.225642,0,"A2\""), mo(0.236579,0,"E1'")
    };

    cov::MODiagramOptions options;
    options.selected_index = 16; // raw MO17
    options.degeneracy.tolerance_hartree = 1.0e-5;
    options.filter.core_energy_cutoff_hartree = -1.5;
    options.filter.virtual_window_hartree = 1.5;
    options.max_levels = 23;
    options.max_virtual_levels = 10;
    options.energy_axis_mode = cov::EnergyAxisMode::NonlinearFocus;
    options.nonlinear_minimum_gap_weight = 0.070;
    options.width = 1200;
    options.height = 900;

    const cov::MODiagramData data = cov::build_mo_diagram_data(wf, options);
    // The fixture has exactly 5 deep-core MOs and 22 valence/frontier MOs.
    if (data.selection.included_indices.size() != 22 || data.selection.hidden_count != 5 ||
        data.selection.valence_occupied_count != 13 || data.selection.frontier_virtual_count != 9) {
        std::cerr << "unexpected Cp valence selection: shown="
                  << data.selection.included_indices.size() << " hidden="
                  << data.selection.hidden_count << " occupied-valence="
                  << data.selection.valence_occupied_count << " frontier-virtual="
                  << data.selection.frontier_virtual_count << '\n';
        return 1;
    }

    if (data.metadata[16].display_label != "17-a" ||
        data.metadata[17].display_label != "17-b") {
        std::cerr << "Cp MO17/MO18 degeneracy labels failed\n";
        return 2;
    }
    if (data.metadata[16].symmetry != "E1\"" || data.metadata[17].symmetry != "E1\"") {
        std::cerr << "Cp E1 double-prime source symmetry was not preserved\n";
        return 3;
    }
    if (data.annotations[16].family != "unavailable" ||
        data.annotations[16].bonding_class != cov::BondingClass::Unclassified) {
        std::cerr << "Cp chemistry was over-inferred from E1 double-prime / occupancy\n";
        return 4;
    }

    // A1' MO21 and E2' MO22 are very close, but incompatible symmetry must
    // keep them non-degenerate; MO22/MO23 remain the E2' degenerate pair.
    if (data.metadata[20].degeneracy_size != 1 ||
        data.metadata[21].display_label != "22-a" ||
        data.metadata[22].display_label != "22-b") {
        std::cerr << "near-but-nondegenerate symmetry gate failed\n";
        return 5;
    }

    const auto notation = cov::parse_symmetry_notation("A1g");
    if (notation.base != "A" || notation.subscript != "1g" || !notation.superscript.empty()) {
        std::cerr << "A1g symmetry formatting parse failed\n";
        return 6;
    }
    const auto prime = cov::parse_symmetry_notation("E2''");
    if (prime.base != "E" || prime.subscript != "2" || prime.superscript != "′′") {
        std::cerr << "E2'' symmetry formatting parse failed\n";
        return 7;
    }

    struct LinearNotation { const char* raw; const char* base; const char* sub; const char* super; };
    for (const auto& expected : std::vector<LinearNotation>{
             {"SG","Σ","",""},{"SGG","Σ","g",""},{"SGU","Σ","u",""},
             {"PI","Π","",""},{"PIG","Π","g",""},{"PIU","Π","u",""},
             {"DLTA","Δ","",""},{"DLTG","Δ","g",""},{"DLTU","Δ","u",""},
             {"PHI","Φ","",""},{"PHIU","Φ","u",""},
             {"Sigma_g^+","Σ","g","+"},{"Σu−","Σ","u","−"},
             {"sigma.g","Σ","g",""},{"Γg","Γ","g",""},
             {"SGunknown","SGunknown","",""}}) {
        const auto parsed=cov::parse_symmetry_notation(expected.raw);
        if (parsed.raw!=expected.raw || parsed.base!=expected.base ||
            parsed.subscript!=expected.sub || parsed.superscript!=expected.super) {
            std::cerr << "linear irrep notation/provenance failed for " << expected.raw << '\n';
            return 17;
        }
    }

    // For the default 1200x900 export, every distinct non-degenerate level in
    // this 6..27 Cp window must receive at least ~50 px vertical separation.
    // Degenerate members remain exactly co-linear by sharing display_coordinate.
    std::vector<double> distinct_coordinates;
    for (const auto& level : data.levels) {
        const double c = cov::energy_display_coordinate(level.layout_energy_hartree,
                                                        data.energy_transform);
        if (distinct_coordinates.empty() ||
            std::abs(c - distinct_coordinates.back()) > 1.0e-12) {
            distinct_coordinates.push_back(c);
        }
    }
    for (std::size_t i = 1; i < distinct_coordinates.size(); ++i) {
        const double pixel_gap = (distinct_coordinates[i] - distinct_coordinates[i - 1]) * 720.0;
        if (pixel_gap < 49.0) {
            std::cerr << "crowded Cp level gap below doubled readability target: "
                      << pixel_gap << " px\n";
            return 8;
        }
    }

    const auto out_dir = std::filesystem::current_path() / "visual_artifacts";
    std::filesystem::create_directories(out_dir);
    const auto base = out_dir / "cp_reference_adaptive_v3";
    const auto result = cov::export_mo_diagram_bundle(wf, options, base);
    if (!result.svg || !result.png || !result.json || !result.csv) {
        std::cerr << "Cp visual export failed: " << result.error << '\n';
        return 9;
    }

    std::ifstream svg_file(base.string() + ".mo.svg", std::ios::binary);
    std::ostringstream svg_buf;
    svg_buf << svg_file.rdbuf();
    const std::string svg = svg_buf.str();
    if (svg.find("Nonlinear energy axis") == std::string::npos ||
        svg.find("data-symmetry=\"E1&quot;\"") == std::string::npos ||
        svg.find("font-size=\"6.3\"") != std::string::npos ||
        svg.find("baseline-shift=") != std::string::npos) {
        std::cerr << "clean SVG / symmetry metadata contract failed\n";
        return 10;
    }
    if (svg.find(">17-a<") != std::string::npos || svg.find(">17-b<") != std::string::npos) {
        std::cerr << "MO numbering leaked into exported figure\n";
        return 11;
    }
    if (svg.find("font-size=\"10\">MO [list] ") == std::string::npos) {
        std::cerr << "missing name must retain an explicit canonical list identity\n";
        return 12;
    }
    if (svg.find(">↑</text>") != std::string::npos ||
        svg.find(">↓</text>") != std::string::npos ||
        svg.find("stroke-linecap=\"round\"") == std::string::npos) {
        std::cerr << "SVG electron arrows must be vector primitives, not font glyphs\n";
        return 13;
    }

    std::ifstream json_file(base.string() + ".mo.json", std::ios::binary);
    std::ostringstream json_buf;
    json_buf << json_file.rdbuf();
    const std::string json = json_buf.str();
    if (json.find("\"transform_name\": \"adaptive-log-gap-v3\"") == std::string::npos ||
        json.find("\"label\": \"17-a\"") == std::string::npos ||
        json.find("\"orbital_family\": \"unavailable\"") == std::string::npos ||
        json.find("\"bonding_class\": \"unclassified\"") == std::string::npos) {
        std::cerr << "Cp machine metadata markers missing\n";
        return 14;
    }

    std::cout << "cp_reference_layout_smoke ok; shown=22 hidden=5; artifacts="
              << out_dir.string() << '\n';
    return 0;
}
