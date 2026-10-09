#include "cov/orbital_ui.hpp"
#include "cov/orbital_inspection_ui.hpp"
#include "cov/nbo_aomo_text.hpp"
#include "cov/nbo_ui.hpp"
#include "cov/validation.hpp"
#include "cov/ui_forensic_helpers.hpp"

#include "cov/mo_diagram.hpp"
#include "cov/composition_display_policy.hpp"
#include "cov/mo_group_display_json.hpp"
#include "cov/mo_diagram_layout.hpp"
#include "cov/orbital_ui_text.hpp"
#include "cov/local_orbital_symmetry.hpp"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace cov::ui {
namespace {

const NboAomoName* canonical_name_for(const Wavefunction& w,std::size_t index,
                                     const OrbitalUIState& state){
    if(state.nbo_ui && state.nbo_ui->aomo.names &&
       index<state.nbo_ui->aomo.names->canonical.size())
        return &state.nbo_ui->aomo.names->canonical[index];
    const auto names=canonical_mo_names(w);
    return index<names->canonical.size()?&names->canonical[index]:nullptr;
}

std::string displayed_canonical_name(const Wavefunction& w,std::size_t index,
                                     const OrbitalUIState& state){
    return canonical_mo_display_label(w,index,canonical_name_for(w,index,state));
}

constexpr ImU32 kSigmaColour = IM_COL32(57, 210, 232, 255);
constexpr ImU32 kPiColour = IM_COL32(226, 93, 220, 255);
constexpr ImU32 kDeltaColour = IM_COL32(244, 155, 65, 255);
constexpr ImU32 kPhiColour = IM_COL32(55, 199, 170, 255);
constexpr ImU32 kBondingColour = IM_COL32(77, 218, 145, 255);
constexpr ImU32 kAntibondingColour = IM_COL32(244, 93, 105, 255);
constexpr ImU32 kNonbondingColour = IM_COL32(235, 181, 65, 255);
constexpr ImU32 kSymmetryColour = IM_COL32(177, 125, 245, 255);
constexpr ImU32 kMulticentreColour = IM_COL32(238, 194, 89, 255);
constexpr ImU32 kNumericColour = IM_COL32(93, 174, 255, 255);
constexpr ImU32 kUnavailableColour = IM_COL32(128, 149, 177, 255);

ImVec4 text_colour(const ImU32 colour) {
    return ImGui::ColorConvertU32ToFloat4(colour);
}

void labelled_value(const char* label, const std::string& value, const ImU32 colour) {
    cov::validation::field(label,value);
    ImGui::Text("%s:", label);
    ImGui::SameLine();
    ImGui::TextColored(text_colour(colour), "%s", value.c_str());
}

bool usable_symmetry_text(const std::string& label) {
    return !label.empty() && label!="?" && label!="N/A" && label!="UND";
}
void labelled_number(const char* label, const std::string& value) {
    cov::validation::field(label,value);
    labelled_value(label, value, kNumericColour);
}

std::string fixed_number(const double value, const int precision) {
    std::ostringstream result;
    result.setf(std::ios::fixed, std::ios::floatfield);
    result.precision(precision);
    result << value;
    return result.str();
}

ImU32 family_colour(const std::string_view family) {
    if (family == "sigma") return kSigmaColour;
    if (family == "pi") return kPiColour;
    if (family == "delta") return kDeltaColour;
    if (family == "phi") return kPhiColour;
    return kUnavailableColour;
}

ImU32 bonding_colour(const BondingClass value) {
    switch (value) {
        case BondingClass::Bonding: return kBondingColour;
        case BondingClass::Antibonding: return kAntibondingColour;
        case BondingClass::Nonbonding: return kNonbondingColour;
        default: return kUnavailableColour;
    }
}

std::string superscript_number_ui(const std::size_t value) {
    static constexpr const char* digits[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
    std::string result;
    for (const char digit : std::to_string(value)) result += digits[digit - '0'];
    return result;
}

std::string subscript_number_ui(const int value) {
    static constexpr const char* digits[] = {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};
    std::string result;
    for (const char digit : std::to_string(std::max(0, value))) result += digits[digit - '0'];
    return result;
}

std::string pi_descriptor_ui(const DelocalisedPiDescriptor& descriptor) {
    std::string result=std::string("Π")+
        superscript_number_ui(descriptor.participating_atoms)+
        subscript_number_ui(static_cast<int>(
            std::lround(descriptor.participating_electrons)));
    const std::string topology_suffix=compact_pi_topology_suffix(descriptor);
    if (!topology_suffix.empty()) result+=' '+topology_suffix;
    return result;
}

const char* pi_topology_value(const DelocalisedPiDescriptor& descriptor,
                              const Language language) {
    if (!descriptor.topology_available) return "N/A";
    switch (descriptor.topology) {
        case DelocalisedPiTopology::Path:
            return orbital_tr(OrbitalText::TopologyPath, language);
        case DelocalisedPiTopology::Cycle:
            return orbital_tr(OrbitalText::TopologyCycle, language);
        case DelocalisedPiTopology::BranchedResonance:
            return orbital_tr(OrbitalText::TopologyBranchedResonance, language);
        case DelocalisedPiTopology::Spiro:
            return orbital_tr(OrbitalText::TopologySpiro, language);
        case DelocalisedPiTopology::HapticMetal:
            return orbital_tr(OrbitalText::TopologyHapticMetal, language);
        case DelocalisedPiTopology::SymmetryDirectSum:
            return orbital_tr(OrbitalText::TopologySymmetryDirectSum, language);
        case DelocalisedPiTopology::MultiChannel:
            return orbital_tr(OrbitalText::TopologyMultiChannel, language);
        default:
            return "N/A";
    }
}

std::string pi_channel_detail(const DelocalisedPiDescriptor& descriptor,
                              const Language language) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(3);
    for (std::size_t index=0;index<descriptor.orientation_channel_details.size();++index) {
        if (index) out << " | ";
        const auto& channel=descriptor.orientation_channel_details[index];
        out << (index+1u) << ": "
            << orbital_tr(OrbitalText::ChannelAtoms, language) << ' ';
        for (std::size_t atom=0;atom<channel.atoms.size();++atom) {
            if (atom) out << '-';
            out << (channel.atoms[atom]+1u);
        }
        out << "; n=(" << channel.direction[0] << ',' << channel.direction[1]
            << ',' << channel.direction[2] << "); "
            << orbital_tr(OrbitalText::Coherence, language) << '='
            << channel.coherence;
    }
    return out.str();
}

const char* intermediate_toggle_label(const Language language) {
    return orbital_tr(OrbitalText::HideIntermediateFrameworkMOs, language);
}

[[maybe_unused]] void draw_pi_ring_evidence(const DelocalisedPiDescriptor& descriptor,
                          const Language language) {
    const auto& graph=descriptor.topology_graph;
    // A molecule-wide connectivity graph is not evidence for a ring in this
    // selected electronic family. Do not emit an empty template for acyclic
    // families (including an ordinary two-centre pi bond).
    if(graph.source==PiTopologyGraphSource::Unavailable ||
       (!descriptor.cyclic_topology && graph.channel_ring_witnesses.empty() &&
        !std::any_of(descriptor.orientation_channel_details.begin(),
                     descriptor.orientation_channel_details.end(),
                     [](const auto& channel){return channel.cyclic;})))return;
    const auto source=graph.source==PiTopologyGraphSource::MayerDistanceModel
        ?OrbitalText::MayerDistanceModel:OrbitalText::CovalentDistanceModel;
    ImGui::Separator();
    ImGui::TextUnformatted(orbital_tr(OrbitalText::RingTopologyEvidence,language));
    const std::string source_text=std::string(orbital_tr(OrbitalText::ConnectivityModel,language))+": "+orbital_tr(source,language);
    cov::validation::field(orbital_tr(OrbitalText::ConnectivityModel,language),orbital_tr(source,language));
    ImGui::TextWrapped("%s",source_text.c_str());
    cov::validation::item("pi.rings."+descriptor.family_id+".source");
    ImGui::TextWrapped("%s",orbital_tr(OrbitalText::RingGraphModelExplanation,language));
    if(graph.channel_ring_witnesses.empty() && descriptor.orientation_channels>1u) {
        ImGui::TextWrapped("%s",orbital_tr(OrbitalText::NoChannelRingWitness,language));
        cov::validation::item("pi.rings."+descriptor.family_id+".scope");
    }
    for(std::size_t i=0;i<graph.channel_ring_witnesses.size();++i) {
        const auto& witness=graph.channel_ring_witnesses[i];
        const auto base="pi.rings."+descriptor.family_id+"."+std::to_string(i);
        labelled_number(orbital_tr(OrbitalText::RingUnionCentre,language),
                        std::to_string(witness.ring_union.hub+1u));
        labelled_number(orbital_tr(OrbitalText::OrientationChannels,language),
                        std::to_string(witness.channel_indices[0]+1u)+" / "+std::to_string(witness.channel_indices[1]+1u));
        for(std::size_t ring=0;ring<2;++ring) {
            std::ostringstream path;
            path<<orbital_tr(OrbitalText::RingPath,language)<<' '<<(ring+1u)<<": ";
            const auto& atoms=witness.ring_union.cycles[ring];
            for(std::size_t j=0;j<atoms.size();++j) {if(j)path<<" - ";path<<(atoms[j]+1u);}
            const auto text=path.str();
            cov::validation::field(orbital_tr(OrbitalText::RingPath,language),text);
            ImGui::TextWrapped("%s",text.c_str());
            cov::validation::item(base+".path"+std::to_string(ring));
        }
        labelled_value(orbital_tr(OrbitalText::AdditionalRingConnection,language),
                       orbital_tr(witness.ring_union.additional_connection_without_hub?OrbitalText::Yes:OrbitalText::No,language),kPiColour);
        cov::validation::item(base+".additional-connection");
    }
}

const char* intermediate_toggle_tooltip(const Language language) {
    return orbital_tr(OrbitalText::HideIntermediateFrameworkMOsTooltip, language);
}

const char* filter_name(const OrbitalFilterMode mode, const Language language) {
    switch (mode) {
        case OrbitalFilterMode::All: return tr(Text::FilterAll, language);
        case OrbitalFilterMode::Occupied: return tr(Text::FilterOccupied, language);
        case OrbitalFilterMode::Virtual: return tr(Text::FilterVirtual, language);
        case OrbitalFilterMode::Core: return tr(Text::FilterCore, language);
        case OrbitalFilterMode::Valence: return tr(Text::FilterValence, language);
        default: return tr(Text::FilterAuto, language);
    }
}

const char* spin_name_ui(const Spin spin, const Language language) {
    return tr(spin == Spin::Beta ? Text::Beta : Text::Alpha, language);
}

const char* family_symbol_ui(const std::string& family) {
    if (family == "sigma") return "σ";
    if (family == "pi") return "π";
    if (family == "delta") return "δ";
    if (family == "phi") return "φ";
    return "N/A";
}

std::string angular_character_text(const OrbitalChannelDistribution& channel) {
    if(channel.status!=ChemistryStatus::Determined && channel.status!=ChemistryStatus::Percentages)return "N/A";
    std::string result;
    const char* labels[]={"σ","π","δ","φ","?"};
    const double values[]={channel.sigma,channel.pi,channel.delta,channel.phi,channel.undetermined};
    for(std::size_t i=0;i<5;++i) {
        if(!std::isfinite(values[i]) || values[i]<=0.0)continue;
        if(!result.empty())result+=" / ";
        result+=std::string(labels[i])+" "+fixed_number(100.0*values[i],1)+"%";
    }
    return result.empty()?"N/A":result;
}

const char* bonding_ui(const BondingClass value, const Language language) {
    return localised_bonding_class(value, language);
}

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool search_matches(const OrbitalMetadata& item, const char* query) {
    if (!query || !*query) return true;
    std::ostringstream haystack;
    haystack << item.raw_mo_number << ' ' << item.display_label << ' '
             << item.symmetry << ' ' << item.symmetry_view.label << ' '
             << (item.spin == Spin::Beta ? "beta" : "alpha") << ' ';
    switch (item.region) {
        case OrbitalRegion::Core: haystack << "core"; break;
        case OrbitalRegion::Valence: haystack << "valence"; break;
        default: haystack << "virtual"; break;
    }
    return lower_ascii(haystack.str()).find(lower_ascii(query)) != std::string::npos;
}

std::optional<std::size_t> previous_occupied(const std::vector<MolecularOrbital>& orbitals,
                                             std::size_t from,
                                             const Spin spin,
                                             const double threshold) {
    if (orbitals.empty()) return std::nullopt;
    if (from > orbitals.size()) from = orbitals.size();
    while (from > 0) {
        --from;
        if (orbitals[from].spin == spin && orbitals[from].occupation > threshold) return from;
    }
    return std::nullopt;
}

std::optional<std::size_t> next_virtual(const std::vector<MolecularOrbital>& orbitals,
                                        std::size_t from,
                                        const Spin spin,
                                        const double threshold) {
    for (std::size_t i = from + 1; i < orbitals.size(); ++i) {
        if (orbitals[i].spin == spin && orbitals[i].occupation <= threshold) return i;
    }
    return std::nullopt;
}

std::optional<std::size_t> frontier_for_selected_spin(const FrontierOrbitals& frontier,
                                                      const Spin spin,
                                                      const bool lumo) {
    if (!frontier.separate_spin_sets) return lumo ? frontier.lumo : frontier.homo;
    if (spin == Spin::Beta) return lumo ? frontier.beta_lumo : frontier.beta_homo;
    return lumo ? frontier.alpha_lumo : frontier.alpha_homo;
}

void quick_nav_button(const char* label,
                      const std::optional<std::size_t> target,
                      OrbitalUIActions& actions,
                      const float width) {
    if (!target) ImGui::BeginDisabled();
    if (ImGui::Button(label, ImVec2(width, 0.0f)) && target) actions.select_orbital = *target;
    if (!target) ImGui::EndDisabled();
}

void energy_unit_combo(OrbitalUIState& state) {
    constexpr EnergyUnit units[] = {
        EnergyUnit::Hartree,
        EnergyUnit::ElectronVolt,
        EnergyUnit::JoulePerMol,
        EnergyUnit::KilojoulePerMol,
        EnergyUnit::CaloriePerMol,
        EnergyUnit::KilocaloriePerMol,
    };
    const bool unit_open=ImGui::BeginCombo("##energy_unit", energy_unit_symbol(state.energy_unit));
    cov::validation::item("browser.unit");
    if (unit_open) {
        for (const EnergyUnit unit : units) {
            const bool selected = state.energy_unit == unit;
            if (ImGui::Selectable(energy_unit_symbol(unit), selected)) state.energy_unit = unit;
            cov::validation::item("browser.unit."+std::to_string(static_cast<int>(unit)));
            if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
}

void filter_combo(OrbitalUIState& state, const Language language) {
    const bool filter_open=ImGui::BeginCombo("##orbital_filter", filter_name(state.filter.mode, language));
    cov::validation::item("browser.filter");
    if (filter_open) {
        constexpr OrbitalFilterMode modes[] = {
            OrbitalFilterMode::AutoReasonable, OrbitalFilterMode::All,
            OrbitalFilterMode::Occupied, OrbitalFilterMode::Virtual,
            OrbitalFilterMode::Core, OrbitalFilterMode::Valence,
        };
        for (const auto mode : modes) {
            const bool selected = state.filter.mode == mode;
            if (ImGui::Selectable(filter_name(mode, language), selected)) state.filter.mode = mode;
            cov::validation::item("browser.filter."+std::to_string(static_cast<int>(mode)));
            if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
}

float rich_symmetry_width(const SymmetryNotation& notation) {
    if (notation.base.empty()) return ImGui::CalcTextSize("N/A").x;
    const float base = ImGui::CalcTextSize(notation.base.c_str()).x;
    const float sub = notation.subscript.empty() ? 0.0f : ImGui::CalcTextSize(notation.subscript.c_str()).x * 0.72f;
    const float sup = notation.superscript.empty() ? 0.0f : ImGui::CalcTextSize(notation.superscript.c_str()).x * 0.72f;
    return base + std::max(sub, sup);
}

void draw_rich_symmetry(const std::string& raw,
                        const ImU32 colour = IM_COL32(220, 228, 240, 255)) {
    const SymmetryNotation notation = parse_symmetry_notation(raw);
    if (notation.base.empty()) {
        ImGui::TextColored(text_colour(kUnavailableColour), "N/A");
        return;
    }
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    ImFont* font = ImGui::GetFont();
    const float size = ImGui::GetFontSize();
    draw->AddText(font, size, pos, colour, notation.base.c_str());
    const float x = pos.x + ImGui::CalcTextSize(notation.base.c_str()).x;
    if (!notation.subscript.empty()) {
        draw->AddText(font, size * 0.72f, ImVec2(x, pos.y + size * 0.36f),
                      colour, notation.subscript.c_str());
    }
    if (!notation.superscript.empty()) {
        draw->AddText(font, size * 0.72f, ImVec2(x, pos.y - size * 0.08f),
                      colour, notation.superscript.c_str());
    }
    ImGui::Dummy(ImVec2(rich_symmetry_width(notation), size * 1.08f));
}

const char* symmetry_scope_tag(const OrbitalSymmetryExplanation& explanation,Language language) {
    if (orbital_symmetry_is_candidate(explanation)) return orbital_tr(OrbitalText::SymmetryCandidateTag,language);
    if (orbital_symmetry_is_local(explanation)) return orbital_tr(OrbitalText::SymmetryLocalTag,language);
    if (explanation.origin==OrbitalSymmetryOrigin::Producer) return orbital_tr(OrbitalText::SymmetrySourceTag,language);
    if (explanation.origin==OrbitalSymmetryOrigin::MolecularOperations) return orbital_tr(OrbitalText::SymmetryMolecularTag,language);
    if (explanation.origin==OrbitalSymmetryOrigin::MixedMembers) return orbital_tr(OrbitalText::SymmetryMixedTag,language);
    return "N/A";
}

const char* symmetry_scope_origin(const OrbitalSymmetryExplanation& e,Language language) {
    switch(e.origin) {
        case OrbitalSymmetryOrigin::LocalMetricProjection: return orbital_tr(OrbitalText::SymmetryProjectionOrigin,language);
        case OrbitalSymmetryOrigin::LocalMetricDecomposition: return orbital_tr(OrbitalText::SymmetryProjectionOrigin,language);
        case OrbitalSymmetryOrigin::LocalDimensionCandidate: return orbital_tr(OrbitalText::SymmetryDimensionOrigin,language);
        case OrbitalSymmetryOrigin::PiPartnerCandidate: return orbital_tr(OrbitalText::SymmetryPiOrigin,language);
        case OrbitalSymmetryOrigin::SpinCounterpartCandidate: return orbital_tr(OrbitalText::SymmetrySpinOrigin,language);
        case OrbitalSymmetryOrigin::MolecularOperations: return orbital_tr(OrbitalText::SymmetryMolecularOrigin,language);
        default: return symmetry_scope_tag(e,language);
    }
}

std::string symmetry_members(const std::vector<std::size_t>& members,const Wavefunction& wf) {
    std::ostringstream out;
    for (const auto index:members) {
        if (out.tellp()>0) out<<", ";
        out<<canonical_mo_display_label(wf,index);
    }
    return out.str();
}

void draw_symmetry_scope_details(const OrbitalSymmetryExplanation& e,const Wavefunction& wf,
                                Language language,const std::string& prefix,
                                const bool show_reference_shell_traces) {
    // Scores, coordinate frames and residuals remain in the evidence export.
    if(!orbital_symmetry_is_local(e) || orbital_symmetry_is_candidate(e))return;
    labelled_value(aomo_text(language,"Local coordination symmetry"),
                   e.label.empty()?"?":e.label,kSymmetryColour);
    cov::validation::item(prefix+".label");
    if(!e.point_group.empty())
        labelled_value(orbital_tr(OrbitalText::PointGroup,language),
                       point_group_display(e.point_group),kSymmetryColour);
    if(!e.orbital_indices.empty())
        labelled_value(orbital_tr(OrbitalText::SymmetryTargetMembers,language),
                       symmetry_members(e.orbital_indices,wf),kNumericColour);
    if(show_reference_shell_traces && e.local_decomposition &&
       e.local_decomposition->status==MetricSubspaceStatus::Available) {
        const auto& projection=*e.local_decomposition;
        if(projection.represented_spin_orbital_rank>0) {
            constexpr const char* shells[]={"s","p","d","f","g"};
            for(std::size_t l=0;l<5;++l) {
                double trace=0;
                for(std::size_t component=0;component<2*l+1;++component)
                    trace+=projection.component_projection_traces[l*l+component];
                if(std::isfinite(trace))
                    labelled_number(shells[l],fixed_number(100.0*trace/
                        static_cast<double>(projection.represented_spin_orbital_rank),3)+"%");
            }
        }
    }
}

int group_member_index(const std::string& label) {
    const auto dash = label.find('-');
    if (dash == std::string::npos || dash + 1 >= label.size()) return 0;
    const char suffix = label[dash + 1];
    return suffix >= 'a' && suffix <= 'z' ? static_cast<int>(suffix - 'a') : 0;
}

std::size_t group_base_raw(const OrbitalMetadata& item) {
    if (item.degeneracy_size <= 1) return item.raw_mo_number;
    const auto member = static_cast<std::size_t>(std::max(0, group_member_index(item.display_label)));
    return item.raw_mo_number >= member ? item.raw_mo_number - member : item.raw_mo_number;
}


std::string delocalised_member_list(const Wavefunction& wavefunction,
                                    const OrbitalUIState& state,
                                    const DelocalisedPiDescriptor& descriptor) {
    std::ostringstream result;
    for (const auto orbital_index : descriptor.orbital_indices) {
        if (result.tellp() > 0) result << ", ";
        result << displayed_canonical_name(wavefunction,orbital_index,state);
    }
    return result.str();
}

std::string delocalised_atom_list(const DelocalisedPiDescriptor& descriptor) {
    std::ostringstream result;
    for (const auto atom_index : descriptor.atom_indices) {
        if (result.tellp() > 0) result << ", ";
        result << atom_index + 1u;
    }
    return result.str();
}

ImU32 pi_interaction_colour(PiInteractionKind kind);

bool pi_relation_visible(const PiInteractionDescriptor& relation,bool brief) {
    if(relation.orbital_evidence && !relation.orbital_evidence->channel.mode_assessment.two_endpoint_relation)
        return false;
    if(!brief)return true;
    return relation.orbital_evidence && !relation.orbital_evidence->weak &&
        relation.orbital_evidence->channel.mode_assessment.ordinary_display_eligible;
}

std::string pi_relation_scope(const PiInteractionDescriptor& relation,Language language) {
    if(!relation.orbital_evidence)return {};
    const auto& mode=relation.orbital_evidence->channel.mode_assessment;
    return aomo_text(language,mode.channel_family=="valence-d-pi"?"Valence d-pi channel":
        mode.channel_family=="auxiliary-p-pi"?"Additional p-pi channel":"Additional radial channel");
}

std::string pi_relation_identity(const PiInteractionDescriptor& relation) {
    if(!relation.orbital_evidence)return {};
    const auto& c=relation.orbital_evidence->channel;
    std::string identity=c.channel_id;
    for(const auto& mode:c.mode_assessment.shared_mode_ids)identity+="|"+mode;
    return identity;
}

std::string pi_relation_mode_json(const PiInteractionDescriptor& relation) {
    if(!relation.orbital_evidence)return "null";
    const auto& m=relation.orbital_evidence->channel.mode_assessment;
    const auto strings=[](const std::vector<std::string>& items) {
        std::string result="[";
        for(std::size_t i=0;i<items.size();++i){if(i)result+=',';result+=validation::quote(items[i]);}
        return result+"]";
    };
    return "{\"verified\":"+std::string(forensic::boolean(m.verified))+
        ",\"direction_verified\":"+forensic::boolean(m.direction_verified)+
        ",\"direction\":"+validation::quote(m.direction)+
        ",\"family\":"+validation::quote(m.channel_family)+
        ",\"reason\":"+validation::quote(m.reason)+
        ",\"multi_group_relation\":"+forensic::boolean(m.multi_group_relation)+
        ",\"two_endpoint_relation\":"+forensic::boolean(m.two_endpoint_relation)+
        ",\"shared_mode_ids\":"+strings(m.shared_mode_ids)+
        ",\"matched_edge_ids\":"+strings(m.matched_edge_ids)+
        ",\"lower_coverage\":"+forensic::number(m.lower_coverage)+
        ",\"upper_coverage\":"+forensic::number(m.upper_coverage)+
        ",\"shared_fragment_contraction_norm_hartree\":"+forensic::number(m.shared_fragment_contraction_norm_hartree)+"}";
}

const char* pi_network_direction_key(const PiModeNetworkAssessment& network) {
    if(network.direction=="ligand_to_centre")return "Pi donation";
    if(network.direction=="centre_to_ligand")return "Pi back-donation";
    if(network.direction=="bidirectional")return "Pi back-donation and donation";
    if(network.direction=="occupied_space_mixing")return "Occupied-space mixing";
    return "Pi mixing";
}

std::string pi_network_scope(const PiModeNetworkAssessment& network,Language language) {
    const char* key=network.ligand_space_kind=="transverse-p-basis"?"Transverse p reference space":
        network.channel_family=="valence-d-pi"?"Valence d-pi channel":
        network.channel_family=="auxiliary-p-pi"?"Additional p-pi channel":"Additional radial channel";
    return aomo_text(language,key);
}

bool pi_network_visible(const PiModeNetworkAssessment& network,bool brief,
                        const std::optional<std::size_t> selected={}) {
    if(!network.verified)return false;
    if(brief && !network.ordinary_display_eligible)return false;
    if(!selected)return true;
    return std::any_of(network.nodes.begin(),network.nodes.end(),[&](const auto& node) {
        return (!brief||node.primary) &&
            std::find(node.members.begin(),node.members.end(),*selected)!=node.members.end();
    });
}

bool pi_network_has_visible_pair(const MODiagramData& data,const PiModeNetworkAssessment& network,
                                bool brief,const std::optional<std::size_t> selected={}) {
    return std::any_of(data.pi_interactions.begin(),data.pi_interactions.end(),[&](const auto& relation) {
        if(!pi_relation_visible(relation,brief)||!relation.orbital_evidence)return false;
        const auto& channel=relation.orbital_evidence->channel;
        if(channel.channel_id!=network.channel_id ||
           std::find(channel.mode_assessment.shared_mode_ids.begin(),channel.mode_assessment.shared_mode_ids.end(),
               network.mode_id)==channel.mode_assessment.shared_mode_ids.end())return false;
        return !selected || std::find(relation.lower_orbitals.begin(),relation.lower_orbitals.end(),*selected)!=relation.lower_orbitals.end() ||
            std::find(relation.upper_orbitals.begin(),relation.upper_orbitals.end(),*selected)!=relation.upper_orbitals.end();
    });
}

std::string coordination_skeleton_name(const MODiagramData& data,
                                      const Wavefunction& wavefunction) {
    if(data.ligand_field_metal_atom>=wavefunction.atoms.size())return {};
    std::string name=wavefunction.atoms[data.ligand_field_metal_atom].symbol;
    std::map<std::string,std::size_t> neighbours;
    for(const auto atom:data.ligand_field_ligand_atoms)
        if(atom<wavefunction.atoms.size())++neighbours[wavefunction.atoms[atom].symbol];
    for(const auto& [symbol,count]:neighbours) {
        name+=symbol;
        if(count>1)name+=std::to_string(count);
    }
    return name;
}

void draw_level_details(const MODiagramData& data,
                        const MODiagramLevel& level,
                        const Wavefunction& wavefunction,
                        const OrbitalUIState& state,
                        const Language language,
                        const std::size_t orbital_index) {
    const OrbitalMetadata& metadata=orbital_index<data.metadata.size()
        ?data.metadata[orbital_index]:level.metadata;
    const auto* verified_name=canonical_name_for(wavefunction,orbital_index,state);
    const std::string displayed_symmetry=canonical_mo_current_irrep(wavefunction,orbital_index,verified_name);
    const auto* actual=orbital_index<wavefunction.orbitals.size()?&wavefunction.orbitals[orbital_index]:nullptr;
    const auto ml=metal_ligand_detail_availability(wavefunction,data,level);
    const auto composition_policy=composition_display_policy(level.composition,ml);
    const auto& selected_annotation=orbital_index<data.annotations.size()
        ?data.annotations[orbital_index]:OrbitalAnnotation{};
    ImGui::TextUnformatted(displayed_canonical_name(wavefunction,orbital_index,state).c_str());
    ImGui::Separator();
    labelled_number(tr(Text::RawMO, language), canonical_mo_source_label(wavefunction,orbital_index));

    labelled_number(data.using_ro_common_energy?aomo_text(language,"Common expectation energy"):
                    tr(Text::ExactEnergy, language),
                    format_energy(metadata.energy_hartree, state.energy_unit, 6));
    if(data.using_ro_common_energy && actual)
        labelled_number(aomo_text(language,"Source orbital energy"),
                        format_energy(actual->energy_hartree,state.energy_unit,6));
    labelled_number(tr(Text::Occupation, language),
        actual && actual->occupation_provenance!=DataProvenance::Unavailable && std::isfinite(actual->occupation)
            ?fixed_number(actual->occupation,6):"N/A");
    if(data.ro_common_energy.restricted_open_shell.verified)
        ImGui::TextUnformatted(aomo_text(language,"Shared orbital space"));
    else ImGui::Text("%s: %s", tr(Text::Spin, language),
        actual && actual->spin_provenance!=DataProvenance::Unavailable?spin_name_ui(actual->spin,language):"N/A");

    if(usable_symmetry_text(displayed_symmetry)) {
        labelled_value(aomo_text(language,"Full MO symmetry"),displayed_symmetry,kSymmetryColour);
        cov::validation::item("details.symmetry.current-label");
    }
    if(usable_symmetry_text(displayed_symmetry) && verified_name &&
       verified_name->verified && !verified_name->point_group.empty())
        labelled_value(verified_name->point_group=="SO(3)"?aomo_text(language,"Rotation group"):orbital_tr(OrbitalText::PointGroup,language),
                       point_group_display(verified_name->point_group),kSymmetryColour);
    // The row may additionally have a local coordination interpretation.
    draw_symmetry_composition(wavefunction,verified_name,state,language,metadata.orbital_index);
    // It is retained as a separate scope, never substituted for the full MO.
    draw_symmetry_scope_details(level.metadata.symmetry_view,wavefunction,language,
        "details.local-symmetry",composition_policy.reference_shell_traces);
    ImGui::Separator();
    ImGui::TextWrapped(orbital_tr(OrbitalText::GroupRepresentativeData,language),
        displayed_canonical_name(wavefunction,level.metadata.orbital_index,state).c_str());
    if (!data.ligand_field_point_group.empty()) {
        labelled_value(aomo_text(language,"Local coordination skeleton"),
                       coordination_skeleton_name(data,wavefunction)+" · "+
                           point_group_display(data.ligand_field_point_group),
                       kSymmetryColour);
        cov::validation::item("details.coordination-skeleton");
        if(ImGui::IsItemHovered() && data.ligand_field_metal_atom<wavefunction.atoms.size()) {
            std::string atoms;
            for(const auto atom:data.ligand_field_ligand_atoms) {
                if(atom>=wavefunction.atoms.size())continue;
                if(!atoms.empty())atoms+=", ";
                atoms+=wavefunction.atoms[atom].symbol+std::to_string(atom+1);
            }
            ImGui::SetTooltip("%s: %s%zu; %s",aomo_text(language,"Centre and coordinating atoms"),
                wavefunction.atoms[data.ligand_field_metal_atom].symbol.c_str(),
                data.ligand_field_metal_atom+1,atoms.c_str());
        }
        if(validation::forensic_mode())
            validation::record("details.coordination-skeleton","{\"centre_atom\":"+
                std::to_string(data.ligand_field_metal_atom)+",\"coordinating_atoms\":"+
                forensic::indices(data.ligand_field_ligand_atoms)+",\"point_group\":"+
                validation::quote(data.ligand_field_point_group)+",\"scope\":\"first-coordination-shell\"}");
        if (!data.ligand_field_geometry_name.empty()) {
            labelled_value(orbital_tr(OrbitalText::CoordinationGeometry, language),
                           localised_geometry_name(
                               data.ligand_field_geometry_id,
                               data.ligand_field_geometry_name,language)+
                               " ("+data.ligand_field_geometry_id+")",
                           kSymmetryColour);
        }
        labelled_number(orbital_tr(OrbitalText::CoordinationNumber, language),
                        std::to_string(data.ligand_field_coordination_number));
    }
    if (!data.local_geometries.empty()) {
        std::set<std::size_t> relevant_centres;
        const auto* routed=state.nbo_ui?state.nbo_ui->routed:nullptr;
        const bool use_nao=routed && orbital_index<routed->mo_composition.size() &&
            routed->mo_composition[orbital_index].available() &&
            routed->mo_composition[orbital_index].provider==RoutedProvider::Nbo;
        if(use_nao) {
            for(const auto& atom:routed->mo_composition[orbital_index].value->atoms)
                if(atom.weight>=0.08)relevant_centres.insert(atom.atom);
        } else if(actual) {
            std::map<std::size_t,double> atom_weights;
            for(const auto& contribution:actual->chemistry.ao_contributions)
                atom_weights[contribution.atom_index]+=contribution.weight;
            for(const auto& [atom,weight]:atom_weights)
                if(weight>=0.08)relevant_centres.insert(atom);
        }
        std::vector<const LocalGeometryDiagramDescriptor*> relevant_geometries;
        for(const auto& geometry:data.local_geometries) {
            if(relevant_centres.contains(geometry.centre_atom) &&
               !geometry.geometry_name.empty() && !geometry.geometry_id.empty() &&
               usable_symmetry_text(geometry.point_group)) {
                relevant_geometries.push_back(&geometry);
            }
        }
        const auto relevant_geometry_count=relevant_geometries.size();
        if(relevant_geometry_count)
            labelled_number(orbital_tr(OrbitalText::LocalMolecularGeometries,language),
                            std::to_string(relevant_geometry_count));
        std::size_t shown=0u;
        for(const auto* geometry_ptr:relevant_geometries) {
            const auto& geometry=*geometry_ptr;
            std::ostringstream value;
            value << orbital_tr(OrbitalText::Atom, language) << ' '
                  << (geometry.centre_atom+1u) << ": "
                  << localised_geometry_name(geometry.geometry_id,
                                             geometry.geometry_name,language)
                  << " (" << geometry.geometry_id << ", "
                  << point_group_display(geometry.point_group) << ", CN"
                  << geometry.neighbour_atoms.size() << ")";
            labelled_value(orbital_tr(OrbitalText::LocalGeometry, language),
                           value.str(),kSymmetryColour);
            if (++shown==6u) break;
        }
        if (shown<relevant_geometry_count) {
            ImGui::TextDisabled("…");
        }
    }

    labelled_number(tr(Text::DegenerateSet, language), std::to_string(level.metadata.degeneracy_size));
    if (level.metadata.degeneracy_size > 1) {
        std::ostringstream members;
        const std::size_t member_count=level.member_indices.empty()
            ?level.metadata.degeneracy_size:level.member_indices.size();
        const std::size_t fallback_base=group_base_raw(level.metadata);
        for (std::size_t i = 0; i < member_count; ++i) {
            if (i) members << ", ";
            const std::size_t member_index=level.member_indices.empty()
                ?fallback_base+i-1u:level.member_indices[i];
            members << displayed_canonical_name(wavefunction,member_index,state);
        }
        labelled_value(tr(Text::DegenerateMembers, language), members.str(), kNumericColour);
    }

    labelled_number(orbital_tr(OrbitalText::GroupOccupation, language),
                    fixed_number(level.total_occupation,3));
    if (composition_policy.metal_ligand) {
        ImGui::TextDisabled("%s",orbital_tr(OrbitalText::MetalLigandGroupAnalysis,language));
        cov::validation::item("details.metal-ligand.scope");
        if(composition_policy.complete_metal_composition && ImGui::IsItemHovered()) {
            const auto atoms=[&](const std::vector<std::size_t>& indices) {
                std::string result;
                for(const auto index:indices) {
                    if(index>=wavefunction.atoms.size())continue;
                    if(!result.empty())result+=", ";
                    result+=wavefunction.atoms[index].symbol+std::to_string(index+1);
                }
                return result;
            };
            ImGui::SetTooltip("%s\n%s: %s\n%s: %s",
                aomo_text(language,"Percentages use the complete MO; group values are averages per member."),
                aomo_text(language,"Metal atoms"),
                atoms(data.composition_scope.centre_atoms).c_str(),
                aomo_text(language,"Ligand atoms"),
                atoms(data.composition_scope.ligand_atoms).c_str());
        }
        if(composition_policy.complete_metal_composition) {
            const auto& c=level.composition;
            const double current_weights[]={c.centre_current_s,c.centre_current_p,
                c.centre_current_d,c.centre_current_f};
            for(int l=0;l<4;++l) {
                if(l==3 && current_weights[l]<0.00005)continue;
                const auto caption=current_composition_shell_caption(wavefunction,
                    data.current_radial_shells,l);
                labelled_number(caption.c_str(),fixed_number(100*current_weights[l],2)+"%");
            }
            cov::validation::item("details.composition.current-metal");
            for(const auto& [label,weight]:std::initializer_list<std::pair<const char*,double>>{
                {"Other metal shells",c.centre_other},
                {"Ligand valence",c.ligand_valence},{"Other ligand space",c.ligand_other},
                {"Core composition",c.core},{"Uncovered composition",c.unresolved}}) {
                if(weight>=0.00005)labelled_number(aomo_text(language,label),fixed_number(100*weight,2)+"%");
            }
            if(c.other_atoms>1e-8)
                labelled_number(aomo_text(language,"Other atom space"),fixed_number(100*c.other_atoms,2)+"%");
        } else if(composition_policy.reference_populations){
        labelled_number(orbital_tr(OrbitalText::MetalSPD, language),fixed_number(100.0*level.metal_s_weight,1)+"% / "+
             fixed_number(100.0*level.metal_p_weight,1)+"% / "+
             fixed_number(100.0*level.metal_d_weight,1)+"%");
        cov::validation::item("details.metal-ligand.populations");
        labelled_number(orbital_tr(OrbitalText::LigandP, language),fixed_number(100.0*level.ligand_p_weight,1)+"%");}
        if(ml.channels){labelled_number(orbital_tr(OrbitalText::MetalLigandSigmaPiChannel, language),
             fixed_number(100.0*level.sigma_fraction,1)+"% / "+fixed_number(100.0*level.pi_fraction,1)+"%");
        cov::validation::item("details.metal-ligand.channels");
        }
        if(ml.overlap){labelled_number(orbital_tr(OrbitalText::MetalLigandOverlap, language),fixed_number(level.metal_ligand_overlap,6));
        cov::validation::item("details.metal-ligand.overlap");}
    }
    if(validation::active()) {
        validation::record("details.composition",mo_group_composition_json(level.composition));
        validation::record("details.display-decision",mo_group_display_decision_json(level.display_decision));
    }
    cov::validation::item("details.metal-ligand.end");
    const auto level_index=static_cast<std::size_t>(&level-data.levels.data());
    const bool brief=state.nbo_ui && state.nbo_ui->aomo.preset==NboAomoPreset::Teaching;
    const auto pi_label=[&](const PiInteractionDescriptor& relation)->std::string {
        if(relation.orbital_evidence && relation.orbital_evidence->channel.direction=="occupied_space_mixing")
            return aomo_text(language,"Occupied-space mixing");
        if(relation.kind==PiInteractionKind::Coupled)return aomo_text(language,"Pi mixing");
        return localised_pi_interaction_kind(relation.kind,language);
    };
    bool has_pi=false;
    std::map<std::string,const PiInteractionDescriptor*> primary_channels;
    for(const auto& interaction:data.pi_interactions) {
        if(!pi_relation_visible(interaction,brief))continue;
        if(interaction.lower_level!=level_index && interaction.upper_level!=level_index &&
           interaction.retained_level!=level_index)continue;
        has_pi=true;
        if(interaction.kind==PiInteractionKind::Donor||interaction.kind==PiInteractionKind::Acceptor)
            primary_channels.emplace(pi_relation_identity(interaction),&interaction);
    }
    std::vector<const PiModeNetworkAssessment*> selected_networks;
    for(const auto& network:data.pi_mode_networks)
        if(pi_network_visible(network,brief,orbital_index) && !pi_network_has_visible_pair(data,network,brief,orbital_index)) {
            selected_networks.push_back(&network);has_pi=true;
        }
    // A physical spin channel is evidence, never an independent whole-ligand
    // verdict. The shared spatial, total-occupation response owns the summary.
    std::vector<PiFieldGroupAssessment> joined_assessments;
    if(ml.composition)for(const auto& response:data.pi_field_response.responses) {
        if(std::none_of(response.centre_atoms.begin(),response.centre_atoms.end(),[&](auto atom){
            return std::find(data.composition_scope.centre_atoms.begin(),
                data.composition_scope.centre_atoms.end(),atom)!=data.composition_scope.centre_atoms.end();}))continue;
        const auto assessment=assess_pi_field_group(response,level.member_indices);
        if(!assessment.available || !assessment.applicable)continue;
        const char* key=nullptr;
        switch(assessment.role) {
            case PiFieldRole::AcceptorDominant:key=assessment.occupied_metal_backbond_supported?
                "Pi back-donation":"Pi acceptor field";break;
            case PiFieldRole::DonorDominant:key="Pi donor field";break;
            case PiFieldRole::Mixed:key="Mixed pi response";break;
            case PiFieldRole::Negligible:key="Pi field below display resolution";break;
            default:break;
        }
        if(key) {
            has_pi=true;
            joined_assessments.push_back(assessment);
            labelled_value(aomo_text(language,"Joined d-pi response"),aomo_text(language,key),kPiColour);
            validation::item("details.pi.joined-summary");
            if(validation::forensic_mode())validation::record("details.pi.joined-summary",
                "{\"label\":"+validation::quote(aomo_text(language,key))+
                ",\"assessment\":"+pi_field_group_assessment_json(assessment)+"}");
        }
    }
    bool show_pi_details=!brief;
    if(brief && has_pi) {
        show_pi_details=ImGui::CollapsingHeader(aomo_text(language,"Pi interaction details"));
        validation::item("details.pi.toggle");
        validation::item(show_pi_details?"details.pi.toggle.open":"details.pi.toggle.closed");
    }
    bool show_local_sources=false;
    if(show_pi_details && has_pi) {
        for(const auto& assessment:joined_assessments) {
            labelled_number(aomo_text(language,"Frozen d-pi field response"),
                format_energy(assessment.mean_shift_hartree,state.energy_unit,6));
            validation::item("details.pi.joined-response");
        }
        // Source records can contain opposite physical spin channels. They
        // remain inspectable, but never compete with the shared-space result
        // as two unqualified whole-ligand conclusions.
        show_local_sources=ImGui::CollapsingHeader(aomo_text(language,"Local channel sources"));
        validation::item("details.pi.sources.toggle");
        validation::item(show_local_sources?"details.pi.sources.open":"details.pi.sources.closed");
    }
    std::set<std::tuple<std::size_t,std::size_t,PiInteractionKind,OrbitalEnergyGapKind,std::string>> shown_relations;
    std::size_t gap_index=0;
    for (const auto* gap:orbital_energy_gaps(data)) {
        const auto& interaction=*gap;
        const bool cf=interaction.gap_kind==OrbitalEnergyGapKind::CrystalField;
        if(!cf && !pi_relation_visible(interaction,brief))continue;
        if(!cf && !show_local_sources)continue;
        const auto prefix="details.energy-gap."+std::to_string(gap_index++);
        if (interaction.lower_level!=level_index &&
            interaction.upper_level!=level_index &&
            interaction.retained_level!=level_index) continue;
        if(!shown_relations.insert({interaction.lower_level,interaction.upper_level,interaction.kind,interaction.gap_kind,
            pi_relation_identity(interaction)}).second)continue;
        ImGui::PushID(static_cast<int>(gap_index));
        if(!cf && state.nbo_ui)for(const auto* members:{&interaction.lower_orbitals,&interaction.upper_orbitals}) {
            if(members->empty())continue;
            const auto target=members->front();
            const auto label=displayed_canonical_name(wavefunction,target,state)+"##counterpart"+std::to_string(target);
            if(ImGui::SmallButton(label.c_str()))state.nbo_ui->focus.pending_canonical_selection=target;
            cov::validation::item(prefix+".endpoint."+std::to_string(target));
        }
        if(!cf) {
            const auto scope=pi_relation_scope(interaction,language);
            const bool multiple=interaction.orbital_evidence &&
                interaction.orbital_evidence->channel.mode_assessment.multi_group_relation;
            ImGui::TextDisabled("%s%s%s",scope.c_str(),multiple?" · ":"",
                multiple?aomo_text(language,"Shared orbital space"):"");
        }
        labelled_value(cf?orbital_tr(OrbitalText::CrystalFieldGap,language):aomo_text(language,"Local channel character"),
                       cf?interaction.symmetry:((interaction.orbital_evidence&&interaction.orbital_evidence->weak)?
                           std::string(aomo_text(language,"Weak interaction"))+" · ":std::string{})+pi_label(interaction),
                       pi_interaction_colour(interaction.kind));
        cov::validation::item(prefix+".kind");
        labelled_number(data.using_ro_common_energy?aomo_text(language,"Common expectation energy gap"):
                            orbital_tr(cf?OrbitalText::CrystalFieldGap:OrbitalText::PiSplitting, language),
                        format_energy(interaction.splitting_hartree,
                                      state.energy_unit,6));
        cov::validation::item(prefix+".splitting");
        // Catalogue priors stay in the structured result, not the details UI.
        if(validation::forensic_mode())
            validation::record("details.energy-gap","{\"kind\":"+std::to_string(static_cast<int>(interaction.kind))+
                ",\"gap_kind\":"+validation::quote(orbital_energy_gap_kind_name(interaction.gap_kind))+
                ",\"label\":"+validation::quote(cf?interaction.symmetry:pi_label(interaction))+
                ",\"lower_orbitals\":"+forensic::indices(interaction.lower_orbitals)+
                ",\"upper_orbitals\":"+forensic::indices(interaction.upper_orbitals)+
                ",\"channel_identity\":"+validation::quote(pi_relation_identity(interaction))+
                ",\"channel_scope\":"+validation::quote(pi_relation_scope(interaction,language))+
                ",\"ordinary_display_eligible\":"+forensic::boolean(pi_relation_visible(interaction,true))+
                ",\"mode\":"+pi_relation_mode_json(interaction)+
                ",\"splitting_hartree\":"+forensic::number(interaction.splitting_hartree)+"}");
        else if(validation::active())
            cov::validation::record("details.energy-gap",orbital_energy_gap_json(interaction,state.energy_unit));
        ImGui::PopID();
    }

    if(show_local_sources)for(const auto* network:selected_networks) {
        ImGui::PushID(network->mode_id.c_str());
        const std::string heading=pi_network_scope(*network,language)+" · "+
            aomo_text(language,"Multi-group mixing")+"##mode-network";
        const bool open=ImGui::TreeNode(heading.c_str());
        validation::item("details.pi-network."+network->mode_id);
        validation::item("details.pi-network."+network->mode_id+(open?".open":".closed"));
        if(validation::forensic_mode())validation::record("details.pi-network",
            "{\"expanded\":"+std::string(forensic::boolean(open))+",\"network\":"+
            pi_mode_network_assessment_json(*network)+"}");
        if(open) {
            // A mode with several participating groups is one network. It has
            // no unique pairwise MO gap and must not be expanded into all pairs.
            const bool physical=network->ligand_space_kind!="transverse-p-basis";
            if(physical)labelled_value(aomo_text(language,"Local channel character"),
                aomo_text(language,pi_network_direction_key(*network)),kPiColour);
            for(const auto& node:network->nodes) {
                if(node.members.empty() || (brief&&!node.primary))continue;
                const auto target=node.members.front();
                if(target>=wavefunction.orbitals.size())continue;
                ImGui::PushID(static_cast<int>(node.group_index));
                const auto label=displayed_canonical_name(wavefunction,target,state);
                if(ImGui::SmallButton(label.c_str()) && state.nbo_ui)
                    state.nbo_ui->focus.pending_canonical_selection=target;
                validation::item("details.pi-network.endpoint."+network->mode_id+"."+std::to_string(target));
                ImGui::SameLine();
                ImGui::TextDisabled("%s: %.2f%% / %.2f%%",aomo_text(language,"Mode contribution: metal / ligand"),
                    100*node.centre_weight,100*node.ligand_weight);
                if(validation::forensic_mode())validation::record("details.pi-network.node",
                    "{\"mode_id\":"+validation::quote(network->mode_id)+",\"members\":"+
                    forensic::indices(node.members)+",\"metal\":"+forensic::number(node.centre_weight)+
                    ",\"ligand\":"+forensic::number(node.ligand_weight)+"}");
                ImGui::PopID();
            }
            ImGui::TreePop();
        }
        ImGui::PopID();
    }

    ImGui::Separator();
    if(actual && (actual->chemistry.channel.status==ChemistryStatus::Determined ||
                  actual->chemistry.channel.status==ChemistryStatus::Percentages))
        labelled_value(aomo_text(language,"All-pair angular character"),
                       angular_character_text(actual->chemistry.channel),kNumericColour);
    if(selected_annotation.bonding_class!=BondingClass::Unclassified)
        labelled_value(aomo_text(language,"Skeleton group bonding character"),bonding_ui(selected_annotation.bonding_class,language),
                       bonding_colour(selected_annotation.bonding_class));

    if (selected_annotation.multicentre.available) {
        labelled_value(tr(Text::MulticentreBond,language),
                       selected_annotation.multicentre.label,kMulticentreColour);
    }
    if (selected_annotation.delocalised_pi.available) {
        const auto& descriptor = selected_annotation.delocalised_pi;
        labelled_value(tr(Text::DelocalisedPiSystem, language), pi_descriptor_ui(descriptor), kPiColour);
        labelled_value(orbital_tr(OrbitalText::MemberMOs, language),
                       delocalised_member_list(wavefunction,state,descriptor), kPiColour);
        std::string atom_value = std::to_string(descriptor.participating_atoms);
        if (!descriptor.atom_indices.empty()) atom_value += " · " + delocalised_atom_list(descriptor);
        labelled_number(orbital_tr(OrbitalText::ParticipatingAtoms, language),atom_value);
        labelled_number(orbital_tr(OrbitalText::ParticipatingElectrons, language),
                        std::to_string(static_cast<int>(std::lround(descriptor.participating_electrons))));
        if(descriptor.topology_available)labelled_number(orbital_tr(OrbitalText::OrientationChannels, language),
                            std::to_string(descriptor.orientation_channels));
        if (descriptor.topology_available &&
            !descriptor.orientation_channel_details.empty()) {
            ImGui::TextWrapped("%s",pi_channel_detail(descriptor,language).c_str());
        }
        if(descriptor.topology_available)labelled_value(orbital_tr(OrbitalText::Topology, language),
                       pi_topology_value(descriptor,language),
                       descriptor.topology_available?kPiColour:kUnavailableColour);
        // Connectivity provenance stays in analysis data, not in a repeated
        // skeleton-ring template in every selected orbital's details.
    }
}

void draw_level_tooltip(const MODiagramData& data,
                        const MODiagramLevel& level,
                        const Wavefunction& wavefunction,
                        const OrbitalUIState& state,
                        const Language language,
                        const std::size_t orbital_index) {
    if(orbital_index>=data.metadata.size() || orbital_index>=wavefunction.orbitals.size())return;
    const auto& metadata=data.metadata[orbital_index];
    const auto& actual=wavefunction.orbitals[orbital_index];
    const ImVec2 work_size=ImGui::GetMainViewport()->WorkSize;
    ImGui::SetNextWindowSizeConstraints(ImVec2(0,0),ImVec2(std::max(120.0f,work_size.x-24.0f),std::max(120.0f,work_size.y-24.0f)));
    ImGui::BeginTooltip();
    ImGui::PushTextWrapPos(ImGui::GetFontSize()*24.0f);
    ImGui::TextUnformatted(displayed_canonical_name(wavefunction,orbital_index,state).c_str());

    labelled_number(tr(Text::ExactEnergy,language),format_energy(metadata.energy_hartree,state.energy_unit,8));
    labelled_number(tr(Text::Occupation,language),
        actual.occupation_provenance!=DataProvenance::Unavailable && std::isfinite(actual.occupation)?fixed_number(actual.occupation,6):"N/A");
    ImGui::Text("%s: %s",tr(Text::Spin,language),
        actual.spin_provenance!=DataProvenance::Unavailable?spin_name_ui(actual.spin,language):"N/A");
    const auto symmetry=canonical_mo_current_irrep(wavefunction,orbital_index,
        canonical_name_for(wavefunction,orbital_index,state));
    if(usable_symmetry_text(symmetry))
        labelled_value(aomo_text(language,"Full MO symmetry"),symmetry,kSymmetryColour);
    draw_symmetry_composition(wavefunction,canonical_name_for(wavefunction,orbital_index,state),
        state,language,orbital_index,{},false);
    if(orbital_index<data.annotations.size()) {
        const auto& annotation=data.annotations[orbital_index];
        labelled_value(aomo_text(language,"All-pair angular character"),
                       angular_character_text(actual.chemistry.channel),kNumericColour);
        if(annotation.bonding_class!=BondingClass::Unclassified)
            labelled_value(aomo_text(language,"Skeleton group bonding character"),bonding_ui(annotation.bonding_class,language),bonding_colour(annotation.bonding_class));
    }
    if(level.member_indices.size()>1)
        ImGui::Text(orbital_tr(OrbitalText::LevelGroupContainsMOs,language),level.member_indices.size());
    ImGui::Separator();
    ImGui::TextWrapped("%s",orbital_tr(OrbitalText::OrbitalDetailsHint,language));
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
}

float map_energy_y(const double energy,
                   const EnergyTransform& transform,
                   const float top,
                   const float bottom) {
    if (transform.knots.empty()) return 0.5f * (top + bottom);
    const double c0 = transform.knots.front().coordinate;
    const double c1 = transform.knots.back().coordinate;
    const double c = energy_display_coordinate(energy, transform);
    const double t = (c - c0) / std::max(1.0e-12, c1 - c0);
    return bottom - static_cast<float>(t) * (bottom - top);
}

double energy_at_fraction(const double fraction, const EnergyTransform& transform) {
    if (transform.knots.empty()) return 0.0;
    const double c0 = transform.knots.front().coordinate;
    const double c1 = transform.knots.back().coordinate;
    return energy_from_display_coordinate(c0 + fraction * (c1 - c0), transform);
}

struct DiagramPoint {
    const MODiagramLevel* level = nullptr;
    ImVec2 left{};
    ImVec2 right{};
    float y = 0.0f;
};

struct DiagramMemberPoint {
    const MODiagramLevel* level = nullptr;
    std::size_t orbital_index = 0;
    ImVec2 left{};
    ImVec2 right{};
    float y = 0.0f;
};

struct PiDiagramGroup {
    const DelocalisedPiDescriptor* descriptor = nullptr;
    std::vector<const DiagramPoint*> points;
};

std::string pi_group_key(const DelocalisedPiDescriptor& descriptor) {
    if (!descriptor.family_id.empty()) return descriptor.family_id;
    std::ostringstream key;
    key << descriptor.label << ':' << descriptor.participating_atoms << ':'
        << std::lround(descriptor.participating_electrons);
    for (const auto atom_index : descriptor.atom_indices) key << ":a" << atom_index;
    for (const auto orbital_index : descriptor.orbital_indices) key << ':' << orbital_index;
    return key.str();
}

void draw_pi_groups(ImDrawList* draw,
                    const std::vector<DiagramPoint>& points,
                    const ImVec2 canvas_max,
                    const float ui_scale) {
    std::map<std::string, PiDiagramGroup> groups;
    for (const auto& point : points) {
        if (!point.level || !point.level->annotation.delocalised_pi.available) continue;
        const auto& descriptor = point.level->annotation.delocalised_pi;
        if (descriptor.participating_atoms < 3u) continue;
        auto& group = groups[pi_group_key(descriptor)];
        if (!group.descriptor) group.descriptor = &descriptor;
        group.points.push_back(&point);
    }

    float bracket_x = canvas_max.x - 13.0f * ui_scale;
    for (const auto& [key, group] : groups) {
        (void)key;
        if (!group.descriptor || group.points.size() < 2u) continue;
        if (!group.descriptor->orbital_indices.empty()) {
            const bool complete=std::all_of(
                group.descriptor->orbital_indices.begin(),
                group.descriptor->orbital_indices.end(),
                [&](const std::size_t orbital_index) {
                    return std::any_of(
                        group.points.begin(),group.points.end(),
                        [&](const DiagramPoint* point) {
                            return point && point->level &&
                                mo_diagram_level_covers_orbital(
                                    *point->level,orbital_index);
                        });
                });
            if (!complete) continue;
        }
        float upper = group.points.front()->y;
        float lower = group.points.front()->y;
        for (const auto* point : group.points) {
            upper = std::min(upper, point->y);
            lower = std::max(lower, point->y);
        }
        upper -= 5.0f * ui_scale;
        lower += 5.0f * ui_scale;
        if (lower - upper < 14.0f * ui_scale) {
            const float middle = 0.5f * (upper + lower);
            upper = middle - 7.0f * ui_scale;
            lower = middle + 7.0f * ui_scale;
        }

        const float cap = 7.0f * ui_scale;
        draw->AddLine(ImVec2(bracket_x, upper), ImVec2(bracket_x, lower), kPiColour, 2.0f * ui_scale);
        draw->AddLine(ImVec2(bracket_x - cap, upper), ImVec2(bracket_x, upper), kPiColour, 2.0f * ui_scale);
        draw->AddLine(ImVec2(bracket_x - cap, lower), ImVec2(bracket_x, lower), kPiColour, 2.0f * ui_scale);

        const std::string label = pi_descriptor_ui(*group.descriptor);
        const ImVec2 label_size = ImGui::CalcTextSize(label.c_str());
        draw->AddText(ImVec2(bracket_x - cap - label_size.x - 5.0f * ui_scale,
                             0.5f * (upper + lower) - 0.5f * label_size.y),
                      kPiColour, label.c_str());
        bracket_x -= std::max(20.0f * ui_scale, label_size.x + 12.0f * ui_scale);
    }
}

ImU32 pi_interaction_colour(const PiInteractionKind kind) {
    switch (kind) {
        case PiInteractionKind::Donor: return IM_COL32(242,163,64,255);
        case PiInteractionKind::Acceptor: return IM_COL32(70,206,218,255);
        case PiInteractionKind::WeakNearNonbonding:
            return IM_COL32(174,154,190,255);
        default: return kPiColour;
    }
}

void draw_ligand_field_pi_interactions(
    ImDrawList* draw,
    const MODiagramData& data,
    const std::vector<DiagramPoint>& points,
    const ImVec2 canvas_max,
    const EnergyUnit unit,
    const Language language,
    const float ui_scale,const bool brief) {
    float bracket_x=canvas_max.x-13.0f*ui_scale;
    std::size_t gap_index=0;
    for (const auto* gap:orbital_energy_gaps(data)) {
        const auto& interaction=*gap;
        const bool cf=interaction.gap_kind==OrbitalEnergyGapKind::CrystalField;
        if(!cf && !pi_relation_visible(interaction,brief))continue;
        if (interaction.retained_level>=points.size()) continue;
        const ImU32 colour=pi_interaction_colour(interaction.kind);
        const char* symbol=cf?"ΔCF":"Δπ";
        const ImVec2 symbol_size=ImGui::CalcTextSize(symbol);
        const float cap=7.0f*ui_scale;
        float hit_top=0.0f;
        float hit_bottom=0.0f;

        if ((!cf && interaction.kind==PiInteractionKind::WeakNearNonbonding) ||
            !interaction.lower_visible || !interaction.upper_visible ||
            interaction.lower_level>=points.size() ||
            interaction.upper_level>=points.size()) {
            const float y=points[interaction.retained_level].y;
            hit_top=y-7.0f*ui_scale;
            hit_bottom=y+7.0f*ui_scale;
            draw->AddLine(ImVec2(bracket_x,hit_top),
                          ImVec2(bracket_x,hit_bottom),colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,hit_top),
                          ImVec2(bracket_x,hit_top),colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,hit_bottom),
                          ImVec2(bracket_x,hit_bottom),
                          colour,2.0f*ui_scale);
        } else {
            float upper=std::min(points[interaction.lower_level].y,
                                 points[interaction.upper_level].y)-5.0f*ui_scale;
            float lower=std::max(points[interaction.lower_level].y,
                                 points[interaction.upper_level].y)+5.0f*ui_scale;
            if (lower-upper<14.0f*ui_scale) {
                const float middle=0.5f*(upper+lower);
                upper=middle-7.0f*ui_scale;
                lower=middle+7.0f*ui_scale;
            }
            hit_top=upper;
            hit_bottom=lower;
            draw->AddLine(ImVec2(bracket_x,upper),ImVec2(bracket_x,lower),
                          colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,upper),ImVec2(bracket_x,upper),
                          colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,lower),ImVec2(bracket_x,lower),
                          colour,2.0f*ui_scale);
        }
        const float symbol_x=bracket_x-cap-symbol_size.x-5.0f*ui_scale;
        const float symbol_y=0.5f*(hit_top+hit_bottom)-0.5f*symbol_size.y;
        draw->AddText(ImVec2(symbol_x,symbol_y),colour,symbol);
        const auto gap_id="diagram.energy-gap."+std::to_string(gap_index++);
        cov::validation::hit(gap_id,ImVec2(symbol_x,hit_top),ImVec2(bracket_x,hit_bottom));
        std::ostringstream gap_draw;
        gap_draw << "{\"gap\":" << orbital_energy_gap_json(interaction,unit)
                 << ",\"symbol\":\"" << symbol << "\",\"symbol_x\":" << symbol_x
                 << ",\"symbol_y\":" << symbol_y << ",\"bracket_x\":" << bracket_x
                 << ",\"top\":" << hit_top << ",\"bottom\":" << hit_bottom << '}';
        cov::validation::record("draw.energy-gap",gap_draw.str());

        if (ImGui::IsItemHovered()) {
            const ImVec2 mouse=ImGui::GetIO().MousePos;
            const float pad=6.0f*ui_scale;
            if (mouse.x>=symbol_x-pad && mouse.x<=bracket_x+pad &&
                mouse.y>=hit_top-pad && mouse.y<=hit_bottom+pad) {
                ImGui::BeginTooltip();
                if (!interaction.symmetry.empty()) {
                    labelled_value(tr(Text::Symmetry, language),
                                   interaction.symmetry,kSymmetryColour);
                }
                 labelled_value(orbital_tr(OrbitalText::Interaction, language),
                                cf?orbital_tr(OrbitalText::CrystalFieldGap,language):
                                    localised_pi_interaction_kind(interaction.kind,language),colour);
                labelled_number(orbital_tr(OrbitalText::Splitting, language),
                    format_energy(interaction.splitting_hartree,unit,6));
                labelled_number(orbital_tr(OrbitalText::SplittingHartree, language),
                    fixed_number(interaction.splitting_hartree,10));
                ImGui::EndTooltip();
            }
        }
        bracket_x-=32.0f*ui_scale;
    }
}

void draw_arrow(ImDrawList* draw, const ImVec2 start, const bool up, const ImU32 colour) {
    const float dy = up ? -13.0f : 13.0f;
    const ImVec2 tip(start.x, start.y + dy);
    draw->AddLine(start, tip, colour, 1.6f);
    const float head = up ? 4.0f : -4.0f;
    draw->AddTriangleFilled(tip, ImVec2(tip.x - 3.5f, tip.y + head),
                            ImVec2(tip.x + 3.5f, tip.y + head), colour);
}

std::string compact_metadata(const Wavefunction& wavefunction,const OrbitalMetadata& item, const EnergyUnit unit,
                             const NboAomoName* verified) {
    std::ostringstream out;
    out << "label=" << canonical_mo_display_label(wavefunction,item.orbital_index,verified)
        << "; source_identity=" << canonical_mo_source_label(wavefunction,item.orbital_index)
        << "; grouping_label=" << item.display_label << "; raw_mo=" << item.raw_mo_number
        << "; internal_index=" << item.orbital_index
        << "; energy_hartree=" << item.energy_hartree
        << "; energy_display=" << convert_hartree(item.energy_hartree, unit)
        << ' ' << energy_unit_symbol(unit)
        << "; occupation=" << item.occupation << "; source symmetry=" << item.symmetry
        << "; view symmetry=" << orbital_symmetry_compact_text(item.symmetry_view)
        << "; point_group=" << item.symmetry_view.point_group
        << "; symmetry_origin=" << orbital_symmetry_origin_name(item.symmetry_view.origin)
        << "; symmetry_source_path=" << item.symmetry_view.source_path
        << "; symmetry_source_lines=" << item.symmetry_view.source_line_begin << '-' << item.symmetry_view.source_line_end
        << "; producer_full_group=" << item.symmetry_view.producer_detected_group
        << "; producer_abelian_group=" << item.symmetry_view.producer_abelian_group
        << "; current display symmetry="
        << canonical_mo_current_irrep(wavefunction,item.orbital_index,verified)
        << "; current display point_group="
        << (verified&&verified->verified?point_group_display(verified->point_group):std::string{});
    if(verified)out<<"; display_name_metadata="<<serialize_orbital_name_json(*verified);
    return out.str();
}

bool same_degeneracy_settings(const DegeneracySettings& left,
                              const DegeneracySettings& right) noexcept {
    return left.tolerance_hartree==right.tolerance_hartree &&
        left.require_same_spin==right.require_same_spin &&
        left.require_compatible_symmetry==right.require_compatible_symmetry &&
        left.maximum_group_size==right.maximum_group_size;
}

bool same_filter_settings(const OrbitalFilterSettings& left,
                          const OrbitalFilterSettings& right) noexcept {
    return left.mode==right.mode &&
        left.occupation_threshold==right.occupation_threshold &&
        left.virtual_window_hartree==right.virtual_window_hartree &&
        left.core_energy_cutoff_hartree==right.core_energy_cutoff_hartree;
}

bool same_diagram_options(const MODiagramOptions& left,
                          const MODiagramOptions& right) noexcept {
    return left.routed_identity==right.routed_identity && left.nbo_source==right.nbo_source &&
        left.use_ro_common_energy==right.use_ro_common_energy &&
        left.ro_common_energy==right.ro_common_energy &&
        left.mode==right.mode &&
        left.aomo_scope==right.aomo_scope &&
        left.show_core_background==right.show_core_background &&
        left.show_fragment_background==right.show_fragment_background &&
        left.display_centre_atoms==right.display_centre_atoms &&
        left.energy_unit==right.energy_unit &&
        left.energy_axis_mode==right.energy_axis_mode &&
        same_degeneracy_settings(left.degeneracy,right.degeneracy) &&
        same_filter_settings(left.filter,right.filter) &&
        left.selected_index==right.selected_index &&
        left.neighbourhood==right.neighbourhood &&
        left.max_levels==right.max_levels &&
        left.max_virtual_levels==right.max_virtual_levels &&
        left.hide_ligand_centred_intermediates==
            right.hide_ligand_centred_intermediates &&
        left.nonlinear_minimum_gap_weight==
            right.nonlinear_minimum_gap_weight &&
        left.weak_pi_split_hartree==right.weak_pi_split_hartree &&
        left.weak_crystal_field_split_hartree==right.weak_crystal_field_split_hartree &&
        left.weak_crystal_field_overlap==right.weak_crystal_field_overlap &&
        left.weak_metal_ligand_overlap==
            right.weak_metal_ligand_overlap &&
        left.width==right.width && left.height==right.height &&
        left.include_hidden_in_metadata==
            right.include_hidden_in_metadata;
}

bool diagram_cache_matches(const OrbitalUIDiagramCache& cache,
                           const Wavefunction& wavefunction,
                           const MODiagramOptions& options) noexcept {
    return cache.data.has_value() && cache.options.has_value() &&
        cache.wavefunction==&wavefunction &&
        cache.orbital_data==wavefunction.orbitals.data() &&
        cache.atom_count==wavefunction.atoms.size() &&
        cache.orbital_count==wavefunction.orbitals.size() &&
        same_diagram_options(*cache.options,options);
}

bool diagram_cache_wavefunction_matches(
    const OrbitalUIDiagramCache& cache,
    const Wavefunction& wavefunction) noexcept {
    return cache.options.has_value() &&
        cache.wavefunction==&wavefunction &&
        cache.orbital_data==wavefunction.orbitals.data() &&
        cache.atom_count==wavefunction.atoms.size() &&
        cache.orbital_count==wavefunction.orbitals.size();
}

bool browser_cache_matches(const OrbitalUIBrowserCache& cache,
                           const Wavefunction& wavefunction,
                           const DegeneracySettings& degeneracy,
                           const OrbitalFilterSettings& filter) noexcept {
    return cache.frontier.has_value() &&
        cache.wavefunction==&wavefunction &&
        cache.orbital_data==wavefunction.orbitals.data() &&
        cache.atom_count==wavefunction.atoms.size() &&
        cache.orbital_count==wavefunction.orbitals.size() &&
        cache.metadata.size()==wavefunction.orbitals.size() &&
        cache.degeneracy.has_value() && cache.filter.has_value() &&
        same_degeneracy_settings(*cache.degeneracy,degeneracy) &&
        same_filter_settings(*cache.filter,filter);
}

void draw_diagram_details_window(const MODiagramData& data,
    const Wavefunction& wavefunction,const std::size_t selected_index,
    OrbitalUIState& state,const Language language,const float ui_scale) {
    if(!state.show_diagram_details) {
        if(validation::forensic_mode())validation::record("forensic.diagram.details",
            "{\"schema\":1,\"open\":false,\"selected_canonical_index\":"+
            std::to_string(selected_index)+",\"active_view\":"+forensic::active(state.active_view)+"}");
        return;
    }
    const auto* viewport=ImGui::GetMainViewport();
    const ImVec2 work=viewport->WorkSize;
    ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x+work.x*.5f,viewport->WorkPos.y+work.y*.5f),
        ImGuiCond_Appearing,ImVec2(.5f,.5f));
    ImGui::SetNextWindowSize(ImVec2(std::min(720.0f*ui_scale,work.x*.8f),work.y*.7f),ImGuiCond_FirstUseEver);
    const ImVec2 margin(std::min(12.0f,work.x*.1f),std::min(12.0f,work.y*.1f));
    const ImVec2 maximum(std::max(1.0f,work.x-2*margin.x),std::max(1.0f,work.y-2*margin.y));
    const ImVec2 minimum(std::min(180.0f,maximum.x),std::min(120.0f,maximum.y));
    ImGui::SetNextWindowSizeConstraints(minimum,maximum);
    if(state.diagram_details_bounds) {
        const auto& bounds=*state.diagram_details_bounds;
        const ImVec2 size(std::clamp(bounds[2],minimum.x,maximum.x),
                          std::clamp(bounds[3],minimum.y,maximum.y));
        const ImVec2 lower(viewport->WorkPos.x+margin.x,viewport->WorkPos.y+margin.y);
        const ImVec2 position(std::clamp(bounds[0],lower.x,lower.x+maximum.x-size.x),
                              std::clamp(bounds[1],lower.y,lower.y+maximum.y-size.y));
        if(position.x!=bounds[0] || position.y!=bounds[1])
            ImGui::SetNextWindowPos(position,ImGuiCond_Always);
    }
    const std::string details_title=std::string(orbital_tr(OrbitalText::OrbitalDetails,language))+"###cov.orbital.details";
    if(state.focus_diagram_details) {
        ImGui::SetNextWindowCollapsed(false);
        ImGui::SetNextWindowFocus();
        state.focus_diagram_details=false;
    }
    const bool visible=ImGui::Begin(details_title.c_str(),
        &state.show_diagram_details,ImGuiWindowFlags_HorizontalScrollbar);
    const auto details_pos=ImGui::GetWindowPos(),details_size=ImGui::GetWindowSize();
    state.diagram_details_bounds=std::array<float,4>{details_pos.x,details_pos.y,details_size.x,details_size.y};
    if(validation::forensic_mode()) {
        const ImVec2 details_max(details_pos.x+details_size.x,details_pos.y+details_size.y);
        const ImVec2 title_lo(details_pos.x+24.0f*ui_scale,details_pos.y);
        const ImVec2 title_hi(std::max(title_lo.x+1.0f,details_max.x-48.0f*ui_scale),
                              details_pos.y+ImGui::GetFrameHeight());
        validation::hit("diagram.details.window",details_pos,details_max);
        validation::chrome_hit("diagram.details.title",title_lo,title_hi);
        const auto row=mo_diagram_row_for_orbital(data,selected_index);
        const std::string route=uses_inspection_details(state)?"inspection":
            row && *row<data.levels.size()?"diagram-level":"canonical-outside-diagram";
        validation::record("forensic.diagram.details","{\"schema\":1,\"open\":true,\"content_visible\":"+
            std::string(forensic::boolean(visible))+",\"title\":"+validation::quote(details_title)+
            ",\"window_hit_id\":\"diagram.details.window\",\"title_hit_id\":\"diagram.details.title\""+
            ",\"window_rect\":"+forensic::rect(details_pos,details_max)+
            ",\"title_hit_rect\":"+forensic::rect(title_lo,title_hi)+
            ",\"selection_route\":"+validation::quote(route)+
            ",\"selected_canonical_index\":"+std::to_string(selected_index)+
            ",\"level_row\":"+forensic::index(row)+
            ",\"group_member_indices\":"+(row && *row<data.levels.size()?
                forensic::indices(data.levels[*row].member_indices):"[]")+
            ",\"active_view\":"+forensic::active(state.active_view)+"}");
    }
    if(visible) {
        if(ImGui::Button(orbital_tr(OrbitalText::CloseOrbitalDetails,language)))state.show_diagram_details=false;
        cov::validation::item("diagram.details.close");
        ImGui::Separator();
        ImGui::PushTextWrapPos(0.0f);
        const auto row=mo_diagram_row_for_orbital(data,selected_index);
        if(uses_inspection_details(state)) {
            draw_inspection_details(wavefunction,state,language);
        } else if(row && *row<data.levels.size()) {
            draw_level_details(data,data.levels[*row],wavefunction,state,language,selected_index);
        } else if(selected_index<wavefunction.orbitals.size()) {
            ImGui::TextUnformatted(displayed_canonical_name(wavefunction,selected_index,state).c_str());
            const auto* name=canonical_name_for(wavefunction,selected_index,state);
            const auto& source=wavefunction.orbitals[selected_index];
            labelled_number(tr(Text::ExactEnergy,language),format_energy(source.energy_hartree,state.energy_unit,8));
            if(source.occupation_provenance!=DataProvenance::Unavailable)
                labelled_number(tr(Text::Occupation,language),fixed_number(source.occupation,6));
            if(name && name->verified)
                labelled_value(aomo_text(language,"Full MO symmetry"),orbital_irrep_display_label(*name),kSymmetryColour);
            draw_symmetry_composition(wavefunction,name,state,language,selected_index);
            ImGui::TextWrapped("%s",orbital_tr(OrbitalText::OrbitalDetailsOutsideDiagram,language));
        }
        ImGui::Separator();
        if(!uses_inspection_details(state))
            ImGui::TextWrapped("%s",orbital_tr(OrbitalText::OrbitalDetailsDataScope,language));
        cov::validation::item("diagram.details.scope");
        const std::string close_bottom=std::string(orbital_tr(OrbitalText::CloseOrbitalDetails,language))+"##bottom";
        if(ImGui::Button(close_bottom.c_str()))state.show_diagram_details=false;
        cov::validation::item("diagram.details.close.bottom");
        ImGui::PopTextWrapPos();
    }
    ImGui::End();
}

} // namespace

void draw_orbital_browser(const Wavefunction& wavefunction,
                          const std::size_t selected_index,
                          OrbitalUIState& state,
                          const Language language,
                          const float ui_scale,
    OrbitalUIActions& actions) {
    if (wavefunction.orbitals.empty()) {
        ImGui::TextDisabled("%s", tr(Text::NoOrbitals, language));
        return;
    }
    state.degeneracy.tolerance_hartree = std::clamp(state.degeneracy.tolerance_hartree, 1.0e-9, 1.0e-2);
    state.filter.virtual_window_hartree = std::clamp(state.filter.virtual_window_hartree, 0.01, 20.0);
    const auto* routed=state.nbo_ui?state.nbo_ui->routed:nullptr;
    const std::string route_identity=routed?
        routed->canonical_fingerprint+":"+routed->integration_id:"";
    if (!browser_cache_matches(state.browser_cache,wavefunction,
                               state.degeneracy,state.filter) ||
        state.browser_cache.routed_identity!=route_identity) {
        state.browser_cache.wavefunction=&wavefunction;
        state.browser_cache.orbital_data=wavefunction.orbitals.data();
        state.browser_cache.atom_count=wavefunction.atoms.size();
        state.browser_cache.orbital_count=wavefunction.orbitals.size();
        state.browser_cache.degeneracy=state.degeneracy;
        state.browser_cache.filter=state.filter;
        state.browser_cache.routed_identity=route_identity;
        state.browser_cache.frontier=find_frontier_orbitals(
            wavefunction.orbitals,state.filter.occupation_threshold);
        // Selection is a transient 3-D inspection state.  None of the table
        // metadata used below depends on its selected flag, so a sentinel
        // keeps the expensive ligand-field classification cacheable while
        // row highlighting continues to use selected_index directly.
        state.browser_cache.metadata=build_orbital_metadata(
            wavefunction,wavefunction.orbitals.size(),
            state.degeneracy,state.filter);
        if(routed && routed->canonical_fingerprint==nbo_canonical_fingerprint(wavefunction)){
            auto& metadata=state.browser_cache.metadata;
            for(std::size_t i=0;i<metadata.size();++i)
                if(const auto balance=routed_mo_shell_balance(*routed,i)){
                    const bool occupied=wavefunction.orbitals[i].occupation>
                        state.filter.occupation_threshold;
                    const bool core=balance->complete&&occupied&&balance->core>=0.70&&
                        balance->core>=balance->valence+0.10;
                    metadata[i].region=core?OrbitalRegion::Core:
                        occupied||balance->valence>=0.05?OrbitalRegion::Valence:
                        OrbitalRegion::Virtual;
                    if(state.filter.mode==OrbitalFilterMode::AutoReasonable||
                       state.filter.mode==OrbitalFilterMode::Valence)
                        metadata[i].visible=!core&&(occupied||balance->valence>=0.05||
                            metadata[i].visible);
                    else if(state.filter.mode==OrbitalFilterMode::Core)
                        metadata[i].visible=core;
                }
            const auto labels=build_orbital_labels(wavefunction.orbitals,
                point_group_limited_degeneracy(wavefunction,state.degeneracy));
            for(std::size_t i=0;i<labels.size();){
                const auto end=std::min(labels.size(),i+std::max<std::size_t>(1,labels[i].group_size));
                bool any=false;for(auto j=i;j<end;++j)any=any||metadata[j].visible;
                if(any)for(auto j=i;j<end;++j)metadata[j].visible=true;
                i=end;
            }
        }
    }
    const FrontierOrbitals& frontier=*state.browser_cache.frontier;
    const auto& metadata=state.browser_cache.metadata;
    std::shared_ptr<const NboAomoNames> standalone_names;
    const NboAomoNames* aomo_names=nullptr;
    if(state.nbo_ui && state.nbo_ui->integration &&
       prepare_nbo_aomo_state(state.nbo_ui->aomo,*state.nbo_ui->integration,wavefunction))
        // The complete MO inventory keeps complete-set names. Diagram names are
        // numbered within its visible set and must not leak into this inventory.
        aomo_names=state.nbo_ui->aomo.source_names.get();
    if(!aomo_names){
        standalone_names=canonical_mo_names(wavefunction);
        aomo_names=standalone_names.get();
    }
    const auto verified_name=[&](std::size_t index)->const NboAomoName* {
        if(!aomo_names || index>=aomo_names->canonical.size())return nullptr;
        const auto& name=aomo_names->canonical[index];
        return &name;
    };
    const Spin selected_spin = selected_index < wavefunction.orbitals.size()
        ? wavefunction.orbitals[selected_index].spin : Spin::Alpha;
    const auto homo = frontier_for_selected_spin(frontier, selected_spin, false);
    const auto lumo = frontier_for_selected_spin(frontier, selected_spin, true);
    const auto hm1 = homo ? previous_occupied(wavefunction.orbitals, *homo, selected_spin,
                                              state.filter.occupation_threshold) : std::nullopt;
    const auto lp1 = lumo ? next_virtual(wavefunction.orbitals, *lumo, selected_spin,
                                         state.filter.occupation_threshold) : std::nullopt;

    const float gap = ImGui::GetStyle().ItemSpacing.x;
    const float bw = std::max(52.0f * ui_scale, (ImGui::GetContentRegionAvail().x - 3.0f * gap) * 0.25f);
    quick_nav_button(tr(Text::HOMOMinus1, language), hm1, actions, bw); ImGui::SameLine();
    quick_nav_button(tr(Text::HOMO, language), homo, actions, bw); ImGui::SameLine();
    quick_nav_button(tr(Text::LUMO, language), lumo, actions, bw); ImGui::SameLine();
    quick_nav_button(tr(Text::LUMOPlus1, language), lp1, actions, bw);

    ImGui::TextDisabled("%s", tr(Text::Search, language));
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputTextWithHint("##orbital_search", tr(Text::Search, language),
                             state.search.data(), state.search.size());
    cov::validation::item("browser.search");

    if (ImGui::BeginTable("##filter_controls", 2, ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings)) {
        ImGui::TableNextColumn(); ImGui::TextDisabled("%s", tr(Text::Filter, language));
        ImGui::SetNextItemWidth(-1.0f); filter_combo(state, language);
        ImGui::TableNextColumn(); ImGui::TextDisabled("%s", tr(Text::EnergyUnit, language));
        ImGui::SetNextItemWidth(-1.0f); energy_unit_combo(state);
        ImGui::EndTable();
    }

    const bool advanced_filters=ImGui::CollapsingHeader(
        (std::string(aomo_text(language,"Filter and grouping settings"))+"###browser.settings").c_str());
    cov::validation::item("browser.settings");
    if(advanced_filters) {
    if (state.filter.mode == OrbitalFilterMode::AutoReasonable) {
        float window = static_cast<float>(state.filter.virtual_window_hartree);
        ImGui::TextDisabled("%s (Ha)", tr(Text::HighVirtualWindow, language));
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::SliderFloat("##virtual_window", &window, 0.05f, 5.0f, "%.2f")) {
            state.filter.virtual_window_hartree = window;
        }
    }

    double tolerance = state.degeneracy.tolerance_hartree;
    ImGui::TextDisabled("%s (Ha)", tr(Text::DegeneracyTolerance, language));
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::InputDouble("##degeneracy_tol", &tolerance, 1.0e-6, 1.0e-5, "%.2e")) {
        state.degeneracy.tolerance_hartree = std::clamp(tolerance, 1.0e-9, 1.0e-2);
    }
    ImGui::Checkbox(tr(Text::GroupedLabels, language), &state.grouped_labels);
    }

    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < metadata.size(); ++i) {
        const auto* name=verified_name(i);
        const auto query=lower_ascii(state.search.data());
        const bool name_match=lower_ascii(canonical_mo_display_label(wavefunction,i,name)).find(query)!=std::string::npos ||
            lower_ascii(canonical_mo_source_label(wavefunction,i)).find(query)!=std::string::npos;
        if (metadata[i].visible && (search_matches(metadata[i], state.search.data()) || name_match)) candidates.push_back(i);
    }
    ImGui::TextDisabled("%s: %zu / %zu", tr(Text::VisibleOrbitals, language), candidates.size(), metadata.size());

    cov::validation::anchor("browser.table");
    const bool table_scroll=candidates.size()>12;
    if (ImGui::BeginTable("##orbital_table", 4,
                          ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH |
                          (table_scroll?ImGuiTableFlags_ScrollY:0) | ImGuiTableFlags_SizingStretchProp |
                          ImGuiTableFlags_NoSavedSettings, ImVec2(0.0f, table_scroll?270.0f*ui_scale:0))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_WidthFixed, 132.0f * ui_scale);
        ImGui::TableSetupColumn(tr(Text::Energy, language), ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn(tr(Text::Occupation, language), ImGuiTableColumnFlags_WidthFixed, 64.0f * ui_scale);
        ImGui::TableSetupColumn(tr(Text::Symmetry, language), ImGuiTableColumnFlags_WidthFixed, 80.0f * ui_scale);
        ImGui::TableHeadersRow();
        ImGuiListClipper clipper;
        clipper.Begin(static_cast<int>(candidates.size()));
        while (clipper.Step()) {
            for (int row = clipper.DisplayStart; row < clipper.DisplayEnd; ++row) {
                const std::size_t index = candidates[static_cast<std::size_t>(row)];
                const auto& item = metadata[index];
                const auto* name=verified_name(index);
                const std::string label=canonical_mo_display_label(wavefunction,index,name);
                ImGui::PushID(static_cast<int>(index));
                ImGui::TableNextRow(0,1.3f*ImGui::GetTextLineHeightWithSpacing()); ImGui::TableNextColumn();
                if (ImGui::Selectable(label.c_str(), index == selected_index,
                                      ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) {
                    actions.select_orbital = index;
                }
                cov::validation::item("browser.mo."+std::to_string(index));
                if(state.grouped_labels && item.degeneracy_size>1)
                    ImGui::TextDisabled("×%zu",item.degeneracy_size);
                ImGui::TableNextColumn();
                ImGui::TextColored(text_colour(kNumericColour), "%s",
                                   format_energy(item.energy_hartree, state.energy_unit, 5).c_str());
                ImGui::TableNextColumn();
                ImGui::TextColored(text_colour(kNumericColour), "%.2f", item.occupation);
                const std::string symmetry=canonical_mo_current_irrep(wavefunction,index,name);
                ImGui::TableNextColumn(); draw_rich_symmetry(symmetry,kSymmetryColour);
                cov::validation::item("browser.symmetry."+std::to_string(index)+".label");
                cov::validation::field("browser.symmetry."+std::to_string(index),symmetry);
                ImGui::PopID();
            }
        }
        ImGui::EndTable();
    }

    if (selected_index < metadata.size() && ImGui::Button(tr(Text::CopyMetadata, language))) {
        const std::string text = uses_inspection_details(state)
            ?inspection_copy_metadata(*state.active_view)
            :compact_metadata(wavefunction,metadata[selected_index],state.energy_unit,verified_name(selected_index));
        ImGui::SetClipboardText(text.c_str());
        validation::record("browser.copy","{\"text\":"+validation::quote(text)+"}");
    }
    validation::item("browser.copy");
}

void draw_pi_counterpart_navigation(const MODiagramData& data,
    const Wavefunction& wavefunction, const OrbitalUIState& state,
    Language language, OrbitalUIActions& actions) {
    const bool brief=state.nbo_ui && state.nbo_ui->aomo.preset==NboAomoPreset::Teaching;
    if(std::none_of(data.pi_interactions.begin(),data.pi_interactions.end(),[&](const auto& p){return pi_relation_visible(p,brief);}) &&
       std::none_of(data.pi_mode_networks.begin(),data.pi_mode_networks.end(),[&](const auto& n){return pi_network_visible(n,brief);}))return;
    const bool open=ImGui::TreeNode(aomo_text(language,"Pi counterparts"));
    validation::item("pi.navigation");
    if(!open)return;
    // Display grouping uses the verified channel ID. Navigation always uses
    // the frozen canonical members; display labels never resolve identities.
    std::map<std::string,std::vector<const PiInteractionDescriptor*>> channels;
    for(const auto& pair:data.pi_interactions) {
        if(!pi_relation_visible(pair,state.nbo_ui && state.nbo_ui->aomo.preset==NboAomoPreset::Teaching))continue;
        if(!pair.orbital_evidence || pair.orbital_evidence->channel.channel_id.empty())continue;
        channels[pair.orbital_evidence->channel.channel_id].push_back(&pair);
    }
    std::size_t ordinal=0;
    for(const auto& [id,pairs]:channels) {
        const auto& channel=pairs.front()->orbital_evidence->channel;
        const bool star=id.find("internal-pi-antibonding")!=std::string::npos;
        const bool internal=id.find("internal-pi-bonding")!=std::string::npos;
        ImGui::Text("%s%s%s",aomo_text(language,star?"Ligand pi*":internal?"Ligand pi":"Metal-ligand pi"),
            channel.spin=="total"?"":" · ",channel.spin=="total"?"":channel.spin.c_str());
        ImGui::TextDisabled("%s",pi_relation_scope(*pairs.front(),language).c_str());
        for(const auto* pair:pairs) {
            ImGui::PushID(id.c_str());ImGui::PushID(static_cast<int>(ordinal));
            for(const bool upper:{false,true}) {
                const auto& members=upper?pair->upper_orbitals:pair->lower_orbitals;
                if(members.empty() || members.front()>=wavefunction.orbitals.size())continue;
                if(upper){ImGui::SameLine();ImGui::TextUnformatted("↔");ImGui::SameLine();}
                const auto label=displayed_canonical_name(wavefunction,members.front(),state)+
                    (upper?"##upper":"##lower");
                if(ImGui::SmallButton(label.c_str())) {
                    actions.select_orbital=members.front();
                    validation::record("pi.navigation.selection",orbital_energy_gap_json(*pair,state.energy_unit));
                }
                validation::item(std::string("pi.navigation.")+(upper?"upper.":"lower.")+std::to_string(ordinal));
            }
            ImGui::SameLine();ImGui::TextDisabled("%s",localised_pi_interaction_kind(pair->kind,language));
            ImGui::PopID();ImGui::PopID();++ordinal;
        }
    }
    for(const auto& network:data.pi_mode_networks) {
        if(!pi_network_visible(network,brief) || pi_network_has_visible_pair(data,network,brief))continue;
        ImGui::PushID(network.mode_id.c_str());
        const auto heading=pi_network_scope(network,language)+" · "+aomo_text(language,"Multi-group mixing");
        if(ImGui::TreeNode(heading.c_str())) {
            for(const auto& node:network.nodes) {
                if(node.members.empty() || (brief&&!node.primary) || node.members.front()>=wavefunction.orbitals.size())continue;
                const auto target=node.members.front();
                ImGui::PushID(static_cast<int>(node.group_index));
                if(ImGui::SmallButton(displayed_canonical_name(wavefunction,target,state).c_str())) {
                    actions.select_orbital=target;
                    validation::record("pi.navigation.network-selection",pi_mode_network_assessment_json(network));
                }
                validation::item("pi.navigation.network."+network.mode_id+"."+std::to_string(target));
                ImGui::PopID();
            }
            ImGui::TreePop();
        }
        ImGui::PopID();
    }
    ImGui::TreePop();
}

void draw_energy_diagram(const Wavefunction& wavefunction,
                         const std::size_t selected_index,
                         OrbitalUIState& state,
                         const Language language,
                         const float ui_scale,
                         OrbitalUIActions& actions) {
    actions.drawn_diagram.reset();
    if (wavefunction.orbitals.empty()) {
        ImGui::TextDisabled("%s", tr(Text::NoOrbitals, language));
        return;
    }
    const bool aomo_ready=state.nbo_ui && state.nbo_ui->integration &&
        nbo_capability(*state.nbo_ui->integration,"aomo") &&
        nbo_capability(*state.nbo_ui->integration,"aomo")->available();

    const char* mo_settings=language==Language::ChineseSimplified?
        "MO 图设置##cov.diagram.settings":language==Language::Japanese?
        "MO 図の設定##cov.diagram.settings":language==Language::French?
        "Réglages du diagramme OM##cov.diagram.settings":
        "MO diagram settings##cov.diagram.settings";
    bool show_mo_settings=!aomo_ready;
    if(aomo_ready) {
        show_mo_settings=ImGui::CollapsingHeader(mo_settings);
        validation::item("diagram.settings");
    }
    if(show_mo_settings) {
    const char* subset_label=language==Language::ChineseSimplified?"轨道范围":
        language==Language::Japanese?"軌道の範囲":
        language==Language::French?"Ensemble orbital":"Orbital scope";
    const char* subset_names=language==Language::ChineseSimplified?
        "价层轨道（σ 与 π）\0仅离域 π 子集\0仅多中心活性空间\0":
        language==Language::Japanese?
        "価電子軌道（σ と π）\0非局在 π 部分集合のみ\0多中心活性空間のみ\0":
        language==Language::French?
        "Orbitales de valence (σ et π)\0Sous-ensemble π délocalisé\0Espace actif multicentrique\0":
        "Valence orbitals (σ and π)\0Delocalised π subset only\0Multicentre active space only\0";
    int subset=state.diagram_mode==MODiagramMode::DelocalisedPiFamilyOnly?1:
        state.diagram_mode==MODiagramMode::MulticentreActiveSpaceOnly?2:0;
    ImGui::SetNextItemWidth(265.0f*ui_scale);
    const char* subset_items[3]={subset_names,nullptr,nullptr};
    subset_items[1]=subset_items[0]+std::strlen(subset_items[0])+1;
    subset_items[2]=subset_items[1]+std::strlen(subset_items[1])+1;
    const bool scope_open=ImGui::BeginCombo(subset_label,subset_items[subset]);
    cov::validation::item("diagram.scope");
    if(scope_open) {
        for(int i=0;i<3;++i) {
            if(ImGui::Selectable(subset_items[i],i==subset))
                state.diagram_mode=i==1?MODiagramMode::DelocalisedPiFamilyOnly:
                    i==2?MODiagramMode::MulticentreActiveSpaceOnly:MODiagramMode::ValenceCentral;
            cov::validation::item("diagram.scope."+std::to_string(i));
            if(i==subset)ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    ImGui::TextDisabled("%s", tr(Text::EnergyScale, language));
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::LinearEnergyScale, language), state.energy_axis_mode == EnergyAxisMode::Linear)) {
        state.energy_axis_mode = EnergyAxisMode::Linear;
    }
    cov::validation::item("diagram.linear");
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::NonlinearFocus, language), state.energy_axis_mode == EnergyAxisMode::NonlinearFocus)) {
        state.energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    }
    cov::validation::item("diagram.nonlinear");
    ImGui::SetNextItemWidth(155.0f * ui_scale);
    ImGui::SliderInt("##diagram_neighbourhood", &state.diagram_neighbourhood, 3, 32, "%d", ImGuiSliderFlags_AlwaysClamp);
    ImGui::SameLine(); ImGui::TextDisabled("%s", tr(Text::AroundSelected, language));
    if (ImGui::Checkbox(intermediate_toggle_label(language),
                        &state.hide_ligand_centred_intermediates)) {
        // The requested options are frozen below in this same frame.
    }
    cov::validation::item("diagram.compact");
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s",intermediate_toggle_tooltip(language));
    }
    }

    MODiagramOptions options;
    if(aomo_ready)prepare_nbo_aomo_state(state.nbo_ui->aomo,*state.nbo_ui->integration,wavefunction);
    const bool compact=aomo_ready?state.nbo_ui->aomo.preset==NboAomoPreset::Teaching:
        state.hide_ligand_centred_intermediates;
    if(aomo_ready) {
        const auto& aomo=state.nbo_ui->aomo;
        options.source_salc_model=aomo.source_salc_model.get();
        options.ro_common_energy=aomo.common_energy_model.get();
        if(options.ro_common_energy && options.ro_common_energy->available) {
            ImGui::Checkbox(aomo_text(language,"Compare RO energies with a common operator"),
                            &state.use_ro_common_energy);
            validation::item("diagram.ro-common-energy");
        }
        options.use_ro_common_energy=state.use_ro_common_energy;
        options.aomo_scope=1+static_cast<unsigned>(aomo.preset);
        options.show_core_background=aomo.show_core;
        options.show_fragment_background=aomo.show_fragment_background;
        if(aomo.salc_model)for(const auto& fragment:aomo.salc_model->fragments)
            if(fragment.side==0 && fragment.atoms.size()==1)
                options.display_centre_atoms.push_back(fragment.atoms.front());
        if(options.display_centre_atoms.size()!=1)options.display_centre_atoms.clear();
    }
    options.mode=state.diagram_mode;
    options.routed=state.nbo_ui?state.nbo_ui->routed:nullptr;
    options.nbo_source=state.nbo_ui?state.nbo_ui->integration:nullptr;
    if(options.routed)options.routed_identity=options.routed->canonical_fingerprint+
        ":"+options.routed->integration_id;
    options.energy_unit = state.energy_unit;
    options.energy_axis_mode = state.energy_axis_mode;
    options.degeneracy = state.degeneracy;
    options.filter = state.filter;
    // A compact ligand-field diagram is a canonical chemical summary.  The
    // MO selected for 3-D inspection must not add, remove or reprioritise its
    // rows.  The live selected_index is still used below for member
    // highlighting and CUDA dispatch; an out-of-range sentinel keeps the
    // expensive structural build independent of inspection state.
    options.selected_index = compact || aomo_ready
        ?wavefunction.orbitals.size():selected_index;
    options.neighbourhood = static_cast<std::size_t>(std::max(3, state.diagram_neighbourhood));
    options.hide_ligand_centred_intermediates=compact;
    options.nonlinear_minimum_gap_weight = 0.055;
    if (!diagram_cache_matches(state.diagram_cache,wavefunction,options)) {
        cov::validation::record("diagram.cache","{\"hit\":false,\"reason\":\"input-or-options-changed\"}");
        if(diagram_cache_wavefunction_matches(state.diagram_cache,wavefunction) &&
           state.diagram_cache.data && state.diagram_cache.options->routed_identity==options.routed_identity &&
           state.diagram_cache.options->nbo_source==options.nbo_source)
            options.pi_field_response=&state.diagram_cache.data->pi_field_response;
        state.diagram_cache.wavefunction=&wavefunction;
        state.diagram_cache.orbital_data=wavefunction.orbitals.data();
        state.diagram_cache.atom_count=wavefunction.atoms.size();
        state.diagram_cache.orbital_count=wavefunction.orbitals.size();
        state.diagram_cache.data=build_mo_diagram_data(wavefunction,options);
        options.pi_field_response=nullptr; // do not retain a pointer into replaced data
        state.diagram_cache.options=options;
        state.diagram_cache.snapshot.reset();
    } else {
        cov::validation::record("diagram.cache","{\"hit\":true}");
    }
    const std::optional<std::size_t> inspected=selected_index<wavefunction.orbitals.size()
        ?std::optional<std::size_t>(selected_index):std::nullopt;
    if (!state.diagram_cache.snapshot ||
        state.diagram_cache.snapshot->data.view->inspected_orbital_index!=inspected) {
        state.diagram_cache.snapshot=std::make_shared<const MODiagramViewSnapshot>(
            make_mo_diagram_view_snapshot(*state.diagram_cache.data,options,
                                         inspected,"interactive-canvas"));
    }
    actions.drawn_diagram=state.diagram_cache.snapshot;
    const auto& snapshot=*actions.drawn_diagram;
    const MODiagramData& data=snapshot.data;
    const MODiagramOptions& drawn_options=snapshot.options;
#ifdef COV_ENABLE_VALIDATION
    cov::validation::record("diagram.snapshot","{\"id\":"+
        cov::validation::quote(data.view->id)+",\"inspected_orbital_index\":"+
        (inspected?std::to_string(*inspected):"null")+
        ",\"selection_anchor\":"+std::to_string(data.view->selection_anchor)+
        ",\"row_count\":"+std::to_string(data.levels.size())+"}");
#endif

    if(!aomo_ready)ImGui::TextDisabled("%s", tr(Text::EnergyDiagram, language));
    const std::string selection_summary=
        localised_diagram_selection_summary(data,language);
    if(!aomo_ready)ImGui::TextDisabled("%s",selection_summary.c_str());
    cov::validation::field("selection_summary",selection_summary);
    bool integrated_aomo=false;
    if(state.nbo_ui && state.nbo_ui->integration) {
        integrated_aomo=draw_nbo_aomo_diagram(state.nbo_ui->aomo,
            *state.nbo_ui->integration,wavefunction,snapshot,language,ui_scale);
        if(!integrated_aomo) {
            ImGui::TextWrapped("%s",aomo_text(language,
                "AO/NAO decomposition unavailable; showing the canonical MO diagram."));
            cov::validation::field("aomo.fallback",state.nbo_ui->aomo.status);
        }
    }
    if(integrated_aomo){
        if(!state.nbo_ui->aomo.drawn_snapshot) {
            actions.drawn_diagram.reset();return;
        }
        if(validation::forensic_mode())validation::record("forensic.diagram",
            "{\"schema\":1,\"snapshot_id\":"+validation::quote(data.view?data.view->id:"")+
            ",\"canvas\":\"aomo\",\"nodes_record\":\"forensic.aomo\""+
            ",\"mode\":"+std::to_string(static_cast<int>(drawn_options.mode))+
            ",\"mode_name\":"+validation::quote(mo_diagram_mode_name(drawn_options.mode))+
            ",\"filter\":"+std::to_string(static_cast<int>(drawn_options.filter.mode))+
            ",\"filter_name\":"+validation::quote(filter_name(drawn_options.filter.mode,language))+
            ",\"virtual_window_hartree\":"+forensic::number(drawn_options.filter.virtual_window_hartree)+
            ",\"core_energy_cutoff_hartree\":"+forensic::number(drawn_options.filter.core_energy_cutoff_hartree)+
            ",\"neighbourhood\":"+std::to_string(drawn_options.neighbourhood)+
            ",\"energy_unit_name\":"+validation::quote(energy_unit_symbol(drawn_options.energy_unit))+
            ",\"compact\":"+forensic::boolean(drawn_options.hide_ligand_centred_intermediates)+
            ",\"selected_canonical_index\":"+std::to_string(selected_index)+
            ",\"active_view\":"+forensic::active(state.active_view)+"}");
        draw_pi_counterpart_navigation(data,wavefunction,state,language,actions);
        // The same canonical details remain reachable from the unified canvas;
        // its levels are already drawn there and must not be drawn a second time.
        if(ImGui::Button(orbital_tr(OrbitalText::OrbitalDetails,language),ImVec2(-1.0f,0.0f)))
            state.show_diagram_details=state.focus_diagram_details=true;
        cov::validation::item("diagram.details");
        draw_diagram_details_window(data,wavefunction,selected_index,state,language,ui_scale);
        return;
    }
    const float height = 430.0f * ui_scale;
    const float left_padding=124.0f*ui_scale;
    const float right_padding=46.0f*ui_scale;
    const float available_width=std::max(1.0f,ImGui::GetContentRegionAvail().x);
    std::vector<DiagramRowFootprint> footprints;
    for (const auto& level:data.levels) {
        const auto members=mo_diagram_member_views(data,level).size();
        const double stroke_width=members==0?0.0:24.0+31.0*static_cast<double>(members-1);
        footprints.push_back({map_energy_y(level.layout_energy_hartree,data.energy_transform,
            34.0f*ui_scale,height-24.0f*ui_scale),14.0*ui_scale,14.0*ui_scale,
            (stroke_width+24.0)*ui_scale});
    }
    const auto lanes=layout_diagram_lanes(footprints,
        std::max(1.0f,available_width-left_padding-right_padding),8.0*ui_scale,2.0*ui_scale);
    const float canvas_width=left_padding+static_cast<float>(lanes.width)+right_padding;
    // Dense physical energies need horizontal space, not a different energy
    // scale or smaller hit targets. Only the canvas scrolls horizontally.
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,ImVec2(0,0));
    ImGui::BeginChild("##energy_diagram_scroll",
        ImVec2(0,height+(canvas_width>available_width?ImGui::GetStyle().ScrollbarSize:0)),
        ImGuiChildFlags_None,ImGuiWindowFlags_HorizontalScrollbar);
    ImGui::PopStyleVar();
    ImGui::InvisibleButton("##energy_diagram_canvas", ImVec2(canvas_width, height), ImGuiButtonFlags_MouseButtonLeft);
    cov::validation::item("diagram.canvas");
    const ImVec2 p0 = ImGui::GetItemRectMin();
    const ImVec2 p1 = ImGui::GetItemRectMax();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(p0, p1, IM_COL32(12, 18, 27, 235), 7.0f * ui_scale);
    draw->AddRect(p0, p1, IM_COL32(43, 58, 77, 220), 7.0f * ui_scale);
    if (data.levels.empty()) {
        if(validation::forensic_mode())validation::record("forensic.diagram",
            "{\"schema\":1,\"snapshot_id\":"+validation::quote(data.view?data.view->id:"")+
            ",\"canvas\":\"canonical\",\"canvas_rect\":"+forensic::rect(p0,p1)+
            ",\"clip_rect\":"+forensic::rect(draw->GetClipRectMin(),draw->GetClipRectMax())+
            ",\"node_count\":0,\"nodes\":[]}");
        ImGui::EndChild();return;
    }

    const float top = p0.y + 34.0f * ui_scale;
    const float bottom = p1.y - 24.0f * ui_scale;
    const float left = p0.x + left_padding;
    const float axis_x = p0.x + 34.0f * ui_scale;
    draw->AddLine(ImVec2(axis_x, bottom), ImVec2(axis_x, top), IM_COL32(129,148,171,255), 1.4f * ui_scale);
    draw->AddTriangleFilled(ImVec2(axis_x, top - 5.0f * ui_scale),
                            ImVec2(axis_x - 4.0f * ui_scale, top + 3.0f * ui_scale),
                            ImVec2(axis_x + 4.0f * ui_scale, top + 3.0f * ui_scale), IM_COL32(129,148,171,255));
    draw->AddText(ImVec2(p0.x + 7.0f * ui_scale, p0.y + 7.0f * ui_scale),
                  IM_COL32(129,148,171,255),
                  drawn_options.energy_axis_mode == EnergyAxisMode::Linear
                    ? tr(Text::LinearEnergyScale, language)
                    : tr(Text::NonlinearEnergyScale, language));

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy = energy_at_fraction(fraction, data.energy_transform);
        const float y = map_energy_y(energy, data.energy_transform, top, bottom);
        draw->AddLine(ImVec2(axis_x - 4.0f * ui_scale, y), ImVec2(axis_x + 4.0f * ui_scale, y),
                      IM_COL32(100,118,140,255), 1.0f);
        const std::string tick = format_energy(energy, drawn_options.energy_unit, 3);
        draw->AddText(ImVec2(axis_x + 6.0f * ui_scale, y - ImGui::GetTextLineHeight() * 0.5f),
                      kNumericColour, tick.c_str());
    }

    std::vector<float> x;
    for (const double centre:lanes.centre_x) x.push_back(left+static_cast<float>(centre));

    std::vector<DiagramPoint> points;
    points.reserve(data.levels.size());
    std::vector<DiagramMemberPoint> member_points;
    member_points.reserve(data.metadata.size());
    const bool capture_forensic=validation::forensic_mode();
    std::string forensic_nodes="[";
    bool first_forensic_node=true;
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const auto& level = data.levels[i];
        const float y = map_energy_y(level.layout_energy_hartree, data.energy_transform, top, bottom);
        const auto member_views=mo_diagram_member_views(data,level);
        const std::size_t members=member_views.size();
        const float member_half=12.0f*ui_scale;
        const float member_spacing=31.0f*ui_scale;
        const float group_half=member_half+
            0.5f*member_spacing*static_cast<float>(members-1u);
        ImU32 colour = IM_COL32(196,210,227,255);
        if (level.metadata.region == OrbitalRegion::Virtual) colour = IM_COL32(126,143,165,255);
        if (level.homo) colour = IM_COL32(79,210,157,255);
        if (level.lumo) colour = IM_COL32(228,176,82,255);
        if (level.annotation.family != "unavailable") colour = family_colour(level.annotation.family);
        for (std::size_t member=0;member<members;++member) {
            const float member_x=x[i]+(static_cast<float>(member)-
                0.5f*static_cast<float>(members-1u))*member_spacing;
            const auto& member_view=member_views[member];
            const std::size_t member_orbital=member_view.orbital_index;
            const std::size_t spin_counterpart=member_view.spin_counterpart.value_or(
                wavefunction.orbitals.size());
            const bool counterpart_selected=member_view.selected &&
                member_view.inspected_orbital_index==spin_counterpart;
            const bool member_selected=member_view.selected;
            if (member_selected) {
                draw->AddLine(ImVec2(member_x-member_half-1.0f*ui_scale,y),
                              ImVec2(member_x+member_half+1.0f*ui_scale,y),
                              kNumericColour,5.0f*ui_scale);
            }
            draw->AddLine(ImVec2(member_x-member_half,y),
                          ImVec2(member_x+member_half,y),colour,
                          member_selected?2.8f:1.8f);
            const ElectronGlyphs electrons=member_view.electrons;
            if(capture_forensic) {
                const auto used=counterpart_selected?spin_counterpart:member_orbital;
                const auto clip_lo=draw->GetClipRectMin(),clip_hi=draw->GetClipRectMax();
                const ImVec2 node_lo(member_x-member_half,y-8.0f*ui_scale),
                    node_hi(member_x+member_half,y+8.0f*ui_scale);
                const ImVec2 hit_lo(std::max(clip_lo.x,node_lo.x),std::max(clip_lo.y,y-4.0f*ui_scale)),
                    hit_hi(std::min(clip_hi.x,node_hi.x),std::min(clip_hi.y,y+4.0f*ui_scale));
                if(!first_forensic_node)forensic_nodes+=',';first_forensic_node=false;
                forensic_nodes+="{\"id\":"+validation::quote("mo:"+std::to_string(used))+
                    ",\"hit_id\":"+validation::quote("diagram.mo."+std::to_string(used))+
                    ",\"canonical_index\":"+std::to_string(used)+
                    ",\"member_canonical_index\":"+std::to_string(member_orbital)+
                    ",\"spin_counterpart\":"+forensic::index(member_view.spin_counterpart)+
                    ",\"source_index\":"+(used<wavefunction.orbitals.size() &&
                        wavefunction.orbitals[used].source_orbital_index!=std::numeric_limits<std::size_t>::max()?
                        std::to_string(wavefunction.orbitals[used].source_orbital_index):"null")+
                    ",\"kind\":\"Canonical\",\"spin\":"+(used<wavefunction.orbitals.size()?
                        validation::quote(wavefunction.orbitals[used].spin==Spin::Beta?"Beta":"Alpha"):"null")+
                    ",\"label\":\"\",\"individual_label\":\"\",\"label_painted\":false"+
                    ",\"row\":"+std::to_string(i)+
                    ",\"member_indices\":"+forensic::indices(level.member_indices)+
                    ",\"member_spin_counterparts\":"+forensic::indices(level.member_spin_counterparts)+
                    ",\"layout_energy_hartree\":"+forensic::number(level.layout_energy_hartree)+
                    ",\"energy_hartree\":"+(used<wavefunction.orbitals.size()?
                        forensic::number(wavefunction.orbitals[used].energy_hartree):"null")+
                    ",\"occupation\":"+(used<wavefunction.orbitals.size() &&
                        wavefunction.orbitals[used].occupation_provenance!=DataProvenance::Unavailable?
                        forensic::number(wavefunction.orbitals[used].occupation):"null")+
                    ",\"alpha_arrows\":"+std::to_string(electrons.alpha)+
                    ",\"beta_arrows\":"+std::to_string(electrons.beta)+
                    ",\"selected_highlight\":"+forensic::boolean(member_selected)+
                    ",\"logical_rect\":"+forensic::rect(ImVec2(node_lo.x-p0.x,node_lo.y-p0.y),
                        ImVec2(node_hi.x-p0.x,node_hi.y-p0.y))+
                    ",\"screen_rect\":"+forensic::rect(node_lo,node_hi)+
                    ",\"hit_rect\":"+(hit_lo.x<hit_hi.x && hit_lo.y<hit_hi.y?forensic::rect(hit_lo,hit_hi):"null")+
                    ",\"viewport_intersects\":"+forensic::boolean(node_hi.x>=clip_lo.x && node_lo.x<=clip_hi.x &&
                        node_hi.y>=clip_lo.y && node_lo.y<=clip_hi.y)+
                    ",\"scroll_clipped\":"+forensic::boolean(node_lo.x<clip_lo.x || node_lo.y<clip_lo.y ||
                        node_hi.x>clip_hi.x || node_hi.y>clip_hi.y)+
                    ",\"clickable_in_viewport\":"+forensic::boolean(hit_lo.x<hit_hi.x && hit_lo.y<hit_hi.y)+"}";
            }
#ifdef COV_ENABLE_VALIDATION
            const auto used_mo=counterpart_selected?spin_counterpart:member_orbital;
            cov::validation::hit("diagram.mo."+std::to_string(used_mo),
                ImVec2(member_x-member_half,y-4.0f*ui_scale),ImVec2(member_x+member_half,y+4.0f*ui_scale));
            std::ostringstream traced;traced<<std::setprecision(17);
            traced<<"{\"used_internal_mo\":"<<used_mo<<",\"used_source_mo\":"<<(used_mo+1)
                <<",\"level_metadata_internal_mo\":"<<level.metadata.orbital_index
                <<",\"y\":"<<y<<",\"layout_energy_hartree\":"<<level.layout_energy_hartree
                <<",\"x\":"<<member_x<<",\"hit_half_width\":"<<(member_half+3.0f*ui_scale)
                <<",\"hit_half_height\":"<<(8.0f*ui_scale)
                <<",\"clip_y\":["<<draw->GetClipRectMin().y<<','<<draw->GetClipRectMax().y<<']'
                <<",\"symmetry\":"<<cov::validation::quote(level.metadata.symmetry_view.label)
                <<",\"symmetry_explanation\":"<<orbital_symmetry_json(level.metadata.symmetry_view)
                <<",\"alpha_arrows\":"<<electrons.alpha<<",\"beta_arrows\":"<<electrons.beta<<"}";
            cov::validation::record("draw.level",traced.str());
#endif
            if (electrons.alpha>0) {
                draw_arrow(draw,ImVec2(member_x-5.0f*ui_scale,y-2.0f),
                           true,IM_COL32(234,242,252,255));
            }
            if (electrons.beta>0) {
                draw_arrow(draw,ImVec2(member_x+5.0f*ui_scale,y-2.0f),
                           false,IM_COL32(234,242,252,255));
            }
            member_points.push_back({&level,
                counterpart_selected?spin_counterpart:member_orbital,
                ImVec2(member_x-member_half,y),
                ImVec2(member_x+member_half,y),y});
        }
        points.push_back({&level,ImVec2(x[i]-group_half,y),
                          ImVec2(x[i]+group_half,y),y});
    }

    draw_pi_groups(draw, points, p1, ui_scale);
    draw_ligand_field_pi_interactions(
        draw,data,points,p1,drawn_options.energy_unit,language,ui_scale,
        state.nbo_ui && state.nbo_ui->aomo.preset==NboAomoPreset::Teaching);
    if(capture_forensic)validation::record("forensic.diagram",
        "{\"schema\":1,\"snapshot_id\":"+validation::quote(data.view?data.view->id:"")+
        ",\"canvas\":\"canonical\",\"mode\":"+std::to_string(static_cast<int>(drawn_options.mode))+
        ",\"mode_name\":"+validation::quote(mo_diagram_mode_name(drawn_options.mode))+
        ",\"filter\":"+std::to_string(static_cast<int>(drawn_options.filter.mode))+
        ",\"filter_name\":"+validation::quote(filter_name(drawn_options.filter.mode,language))+
        ",\"virtual_window_hartree\":"+forensic::number(drawn_options.filter.virtual_window_hartree)+
        ",\"core_energy_cutoff_hartree\":"+forensic::number(drawn_options.filter.core_energy_cutoff_hartree)+
        ",\"neighbourhood\":"+std::to_string(drawn_options.neighbourhood)+
        ",\"compact\":"+forensic::boolean(drawn_options.hide_ligand_centred_intermediates)+
        ",\"energy_unit\":"+std::to_string(static_cast<int>(drawn_options.energy_unit))+
        ",\"energy_unit_name\":"+validation::quote(energy_unit_symbol(drawn_options.energy_unit))+
        ",\"energy_axis_mode\":"+std::to_string(static_cast<int>(drawn_options.energy_axis_mode))+
        ",\"selected_canonical_index\":"+std::to_string(selected_index)+
        ",\"active_view\":"+forensic::active(state.active_view)+
        ",\"canvas_rect\":"+forensic::rect(p0,p1)+
        ",\"clip_rect\":"+forensic::rect(draw->GetClipRectMin(),draw->GetClipRectMax())+
        ",\"scroll\":["+forensic::number(ImGui::GetScrollX())+","+forensic::number(ImGui::GetScrollY())+"]"+
        ",\"scroll_max\":["+forensic::number(ImGui::GetScrollMaxX())+","+forensic::number(ImGui::GetScrollMaxY())+"]"+
        ",\"node_count\":"+std::to_string(member_points.size())+
        ",\"nodes\":"+forensic_nodes+"]}");

    if (ImGui::IsItemHovered()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const DiagramMemberPoint* nearest = nullptr;
        float best = 1.0e9f;
        for (const auto& point : member_points) {
            if (mouse.x < point.left.x - 3.0f * ui_scale || mouse.x > point.right.x + 3.0f * ui_scale) continue;
            const float distance = std::abs(mouse.y - point.y);
            if (distance < best && distance <= 8.0f * ui_scale) { best = distance; nearest = &point; }
        }
        if (nearest && nearest->level) {
            draw_level_tooltip(data,*nearest->level,wavefunction,state,language,
                               nearest->orbital_index);
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                actions.select_orbital=nearest->orbital_index;
            }
        }
    }

    ImGui::EndChild();
    if (ImGui::Button(tr(Text::ExportBundle, language), ImVec2(-1.0f, 0.0f))) actions.export_diagram = true;
    cov::validation::item("diagram.export");
    const bool export_options=ImGui::CollapsingHeader(aomo_text(language,"Analysis data (advanced)"));
    cov::validation::item("diagram.export_options");
    if(export_options){
        if(ImGui::Button(aomo_text(language,"Export analysis data"))) actions.export_analysis=true;
        cov::validation::item("diagram.export_data");
    }
    if(ImGui::Button(orbital_tr(OrbitalText::OrbitalDetails,language),ImVec2(-1.0f,0.0f)))
        state.show_diagram_details=state.focus_diagram_details=true;
    cov::validation::item("diagram.details");
    if(!integrated_aomo && state.nbo_ui && state.nbo_ui->dataset)
        draw_nbo_focus_view(state.nbo_ui->focus,*state.nbo_ui->dataset,
                            wavefunction,snapshot,language,ui_scale);
    draw_diagram_details_window(data,wavefunction,selected_index,state,language,ui_scale);
}

} // namespace cov::ui
