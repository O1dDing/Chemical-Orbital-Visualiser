#include "cov/mo_diagram.hpp"
#include "cov/orbital_chemistry.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr std::array<std::size_t,5> kSyntheticExpectedMembers{
    15u,16u,17u,41u,42u // raw MO16, MO17/MO18 and MO42/MO43
};

cov::Primitive primitive(const float exponent) {
    cov::Primitive value;
    value.exponent=exponent;
    value.coefficient=1.0f;
    return value;
}

void append_shell(cov::Wavefunction& wf,
                  const std::uint32_t atom,
                  const std::uint8_t angular_momentum,
                  const std::uint8_t pure,
                  const float exponent) {
    const auto primitive_offset=
        static_cast<std::uint32_t>(wf.primitives.size());
    wf.primitives.push_back(primitive(exponent));
    wf.shells.push_back({
        atom,primitive_offset,1u,wf.basis_count,angular_momentum,pure});
    wf.basis_count+=cov::shell_basis_count(wf.shells.back());
}

cov::Wavefunction make_planar_c5_case(const bool distort_out_of_plane) {
    cov::Wavefunction wf;
    constexpr double pi=3.1415926535897932384626433832795;
    constexpr double radius=2.25; // bohr; nearest-neighbour C--C = 2.65 bohr
    for (std::size_t atom=0;atom<5u;++atom) {
        const double angle=2.0*pi*static_cast<double>(atom)/5.0;
        wf.atoms.push_back({
            "C",6,radius*std::cos(angle),radius*std::sin(angle),0.0,6.0});
    }
    if (distort_out_of_plane) {
        // Deliberately far outside any chemically reasonable Cp ring plane.
        // Keeping all Mayer edges present makes this a direct planarity gate,
        // rather than accidentally passing because the graph was disconnected.
        wf.atoms[0].z=4.0;
    }

    // One minimal C 2s/2p shell set per centre. The p shell is Cartesian in
    // x,y,z order, so the five z functions are the ring p_perp reference.
    std::array<std::size_t,5> pz{};
    for (std::uint32_t atom=0;atom<5u;++atom) {
        append_shell(wf,atom,0u,0u,1.2f);
        append_shell(wf,atom,1u,0u,0.8f);
        pz[atom]=wf.shells.back().basis_offset+2u;
    }

    // Add 23 polarization/diffuse degrees of freedom so there are 43 complete
    // canonical MOs. These functions must not become extra carbon p-valence
    // references: three spherical f shells contribute 21, plus two diffuse s.
    for (std::uint32_t atom=0;atom<3u;++atom) {
        append_shell(wf,atom,3u,1u,0.35f+0.03f*static_cast<float>(atom));
    }
    append_shell(wf,0u,0u,0u,0.04f);
    append_shell(wf,1u,0u,0u,0.04f);

    if (wf.basis_count!=43u) {
        std::cerr<<"internal fixture error: expected 43 basis functions, got "
                 <<wf.basis_count<<'\n';
        std::abort();
    }
    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) wf.ao_overlap[i*n+i]=1.0;
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;

    std::vector<std::size_t> complement;
    complement.reserve(n-5u);
    for (std::size_t basis=0;basis<n;++basis) {
        if (std::find(pz.begin(),pz.end(),basis)==pz.end()) {
            complement.push_back(basis);
        }
    }

    wf.orbitals.resize(n);
    std::size_t filler=0u;
    for (std::size_t i=0;i<n;++i) {
        auto& mo=wf.orbitals[i];
        mo.coefficients.assign(n,0.0f);
        mo.symmetry="A"+std::to_string(i+1u);
        if (i<15u) mo.energy_hartree=-1.50+0.05*static_cast<double>(i);
        else if (i==15u) mo.energy_hartree=-0.70;
        else if (i<=17u) mo.energy_hartree=-0.50;
        else if (i<41u) mo.energy_hartree=-0.30+
            0.018*static_cast<double>(i-18u);
        else mo.energy_hartree=0.30;
        mo.occupation=i<=17u?2.0f:0.0f;
        mo.occupation_provenance=cov::DataProvenance::Producer;

        if (std::find(
                kSyntheticExpectedMembers.begin(),
                kSyntheticExpectedMembers.end(),i)==
            kSyntheticExpectedMembers.end()) {
            mo.coefficients[complement[filler++]]=1.0f;
        }
    }

    // Real orthonormal SALCs of a regular five-membered p_perp ring:
    // k=0; cos/sin(k=1); cos/sin(k=2). Their occupations are 2,2,2,0,0.
    for (std::size_t member=0;
         member<kSyntheticExpectedMembers.size();++member) {
        auto& mo=wf.orbitals[kSyntheticExpectedMembers[member]];
        for (std::size_t atom=0;atom<5u;++atom) {
            const double angle=2.0*pi*static_cast<double>(atom)/5.0;
            double coefficient=0.0;
            if (member==0u) coefficient=1.0/std::sqrt(5.0);
            else if (member==1u) coefficient=std::sqrt(2.0/5.0)*std::cos(angle);
            else if (member==2u) coefficient=std::sqrt(2.0/5.0)*std::sin(angle);
            else if (member==3u) coefficient=std::sqrt(2.0/5.0)*std::cos(2.0*angle);
            else coefficient=std::sqrt(2.0/5.0)*std::sin(2.0*angle);
            mo.coefficients[pz[atom]]=static_cast<float>(coefficient);
        }
    }
    wf.orbitals[15].symmetry="A2\"";
    wf.orbitals[16].symmetry="E1\"";
    wf.orbitals[17].symmetry="E1\"";
    wf.orbitals[41].symmetry="E2\"";
    wf.orbitals[42].symmetry="E2\"";

    // Explicit ring-only connectivity. Non-neighbouring carbon pairs remain
    // below the graph threshold even though the pentagon is compact.
    for (std::uint32_t atom=0;atom<5u;++atom) {
        wf.bond_orders.push_back({
            atom,static_cast<std::uint32_t>((atom+1u)%5u),
            1.20,cov::DataProvenance::Derived});
    }
    wf.bond_order_provenance=cov::DataProvenance::Derived;
    return wf;
}

cov::Wavefunction make_rooted_cycle_case(const bool equivalent_rooted_rings,
                                         const bool asymmetric_hub_shell) {
    cov::Wavefunction wf;
    constexpr double pi=3.1415926535897932384626433832795;
    constexpr double ring_radius=2.65;
    constexpr double centre_offset=5.50;
    wf.atoms.push_back({"B",5,0.0,0.0,0.0,5.0});

    std::array<std::array<std::uint32_t,6>,2> rings{};
    std::array<std::array<std::size_t,6>,2> pz{};
    for (std::size_t ring=0u;ring<2u;++ring) {
        const double centre_x=ring==0u?centre_offset:-centre_offset;
        const double root_angle=ring==0u?pi:0.0;
        for (std::size_t position=0u;position<6u;++position) {
            const double angle=root_angle+
                2.0*pi*static_cast<double>(position)/6.0;
            const auto atom=static_cast<std::uint32_t>(wf.atoms.size());
            rings[ring][position]=atom;
            wf.atoms.push_back({
                "C",6,centre_x+ring_radius*std::cos(angle),
                ring_radius*std::sin(angle),0.0,6.0});
            append_shell(wf,atom,0u,0u,1.2f);
            append_shell(wf,atom,1u,0u,0.8f);
            pz[ring][position]=wf.shells.back().basis_offset+2u;
        }
    }

    // Both roots see one F at graph distance 1 and one at distance 2.  The
    // first ring has adjacent substituents {1,2}; the second has {1,4}, which
    // is not related by a rooted rotation/reflection.  A distance histogram
    // therefore collides while an ordered rooted-cycle signature must not.
    std::array<std::array<std::size_t,2>,2> substitutions{{
        {{1u,2u}},{{1u,4u}},
    }};
    if (equivalent_rooted_rings) substitutions[1]={{1u,2u}};
    for (std::size_t ring=0u;ring<2u;++ring) {
        for (const auto position:substitutions[ring]) {
            const auto carbon=rings[ring][position];
            const double centre_x=ring==0u?centre_offset:-centre_offset;
            const double dx=wf.atoms[carbon].x-centre_x;
            const double dy=wf.atoms[carbon].y;
            const double norm=std::hypot(dx,dy);
            const auto fluorine=static_cast<std::uint32_t>(wf.atoms.size());
            wf.atoms.push_back({
                "F",9,wf.atoms[carbon].x+1.75*dx/norm,
                wf.atoms[carbon].y+1.75*dy/norm,0.0,9.0});
            wf.bond_orders.push_back(
                {carbon,fluorine,0.30,cov::DataProvenance::Derived});
        }
    }
    for (std::size_t ring=0u;ring<2u;++ring) {
        for (std::size_t position=0u;position<6u;++position) {
            wf.bond_orders.push_back({
                rings[ring][position],rings[ring][(position+1u)%6u],
                1.20,cov::DataProvenance::Derived});
        }
        wf.bond_orders.push_back(
            {0u,rings[ring][0],0.70,cov::DataProvenance::Derived});
    }
    // Complete a genuinely saturated, three-dimensional common-hub first
    // shell.  The unique reflection exchanging the two rooted rings also
    // exchanges the first two terminal contacts.  Giving those contacts
    // different atom colours is therefore a direct negative control for a
    // C(aryl)2XY-like asymmetric hub environment; identical colours retain a
    // valid whole-molecule symmetry operation.
    const auto hub_terminal_a=static_cast<std::uint32_t>(wf.atoms.size());
    wf.atoms.push_back({"F",9,1.30,2.20,1.70,9.0});
    const auto hub_terminal_b=static_cast<std::uint32_t>(wf.atoms.size());
    wf.atoms.push_back({asymmetric_hub_shell?"Cl":"F",
                        asymmetric_hub_shell?17:9,
                        -1.30,-2.20,1.70,
                        asymmetric_hub_shell?17.0:9.0});
    const auto hub_terminal_c=static_cast<std::uint32_t>(wf.atoms.size());
    wf.atoms.push_back({"F",9,0.0,0.0,-2.80,9.0});
    wf.bond_orders.push_back(
        {0u,hub_terminal_a,0.70,cov::DataProvenance::Derived});
    wf.bond_orders.push_back(
        {0u,hub_terminal_b,0.70,cov::DataProvenance::Derived});
    wf.bond_orders.push_back(
        {0u,hub_terminal_c,0.70,cov::DataProvenance::Derived});
    wf.bond_order_provenance=cov::DataProvenance::Derived;

    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0u;i<n;++i) wf.ao_overlap[i*n+i]=1.0;
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;
    wf.orbitals.resize(n);
    for (std::size_t basis=0u;basis<n;++basis) {
        auto& mo=wf.orbitals[basis];
        mo.coefficients.assign(n,0.0f);
        mo.coefficients[basis]=1.0f;
        mo.energy_hartree=-1.50+0.01*static_cast<double>(basis);
        mo.occupation=2.0f;
        mo.occupation_provenance=cov::DataProvenance::Producer;
    }
    for (std::size_t ring=0u;ring<2u;++ring) {
        for (std::size_t mode=0u;mode<6u;++mode) {
            auto& mo=wf.orbitals[pz[ring][mode]];
            std::fill(mo.coefficients.begin(),mo.coefficients.end(),0.0f);
            for (std::size_t position=0u;position<6u;++position) {
                const double angle=2.0*pi*static_cast<double>(position)/6.0;
                double coefficient=0.0;
                if (mode==0u) coefficient=1.0/std::sqrt(6.0);
                else if (mode==1u) coefficient=std::sqrt(2.0/6.0)*std::cos(angle);
                else if (mode==2u) coefficient=std::sqrt(2.0/6.0)*std::sin(angle);
                else if (mode==3u) coefficient=std::sqrt(2.0/6.0)*std::cos(2.0*angle);
                else if (mode==4u) coefficient=std::sqrt(2.0/6.0)*std::sin(2.0*angle);
                else coefficient=std::cos(3.0*angle)/std::sqrt(6.0);
                mo.coefficients[pz[ring][position]]=static_cast<float>(
                    coefficient);
            }
            mo.energy_hartree=-0.60+0.12*static_cast<double>(mode);
            mo.occupation=mode<3u?2.0f:0.0f;
        }
    }
    return wf;
}

cov::Wavefunction make_uhf_o2_case() {
    cov::Wavefunction wf;
    wf.atoms={
        {"O",8,-1.15,0.0,0.0,8.0},
        {"O",8, 1.15,0.0,0.0,8.0}};
    for (std::uint32_t atom=0;atom<2u;++atom) {
        append_shell(wf,atom,0u,0u,1.2f);
        append_shell(wf,atom,1u,0u,0.8f);
    }
    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) wf.ao_overlap[i*n+i]=1.0;
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;

    // Complete orthonormal spatial set: s+/s-, px+/-, py+/-, pz+/-.  Both
    // spin blocks contain the same complete basis; the pi bonding and
    // antibonding pairs are exactly degenerate within each spin.
    const std::array<std::array<std::size_t,2>,8> pairs{{
        {0u,4u},{0u,4u},{1u,5u},{2u,6u},
        {3u,7u},{1u,5u},{2u,6u},{3u,7u}}};
    const std::array<double,8> energies{
        -1.0,-0.82,-0.70,-0.55,-0.55,0.35,0.10,0.10};
    wf.orbitals.resize(2u*n);
    for (std::size_t spin_block=0;spin_block<2u;++spin_block) {
        for (std::size_t spatial=0;spatial<n;++spatial) {
            auto& mo=wf.orbitals[spin_block*n+spatial];
            mo.spin=spin_block==0u?cov::Spin::Alpha:cov::Spin::Beta;
            mo.energy_hartree=energies[spatial]+0.01*spin_block;
            mo.symmetry=(spatial==3u || spatial==4u)?"Pi_u":
                        ((spatial==6u || spatial==7u)?"Pi_g":"Sigma");
            mo.coefficients.assign(n,0.0f);
            const double sign=spatial==1u || spatial>=5u?-1.0:1.0;
            mo.coefficients[pairs[spatial][0]]=
                static_cast<float>(1.0/std::sqrt(2.0));
            mo.coefficients[pairs[spatial][1]]=
                static_cast<float>(sign/std::sqrt(2.0));
            const bool alpha_occupied=spatial<=4u || spatial>=6u;
            const bool beta_occupied=spatial<=4u;
            mo.occupation=(spin_block==0u?alpha_occupied:beta_occupied)
                ?1.0f:0.0f;
            mo.occupation_provenance=cov::DataProvenance::Producer;
        }
    }
    wf.bond_orders.push_back(
        {0u,1u,1.5,cov::DataProvenance::Derived});
    wf.bond_order_provenance=cov::DataProvenance::Derived;
    return wf;
}

cov::Wavefunction make_branched_p_case(const bool unrestricted) {
    cov::Wavefunction wf;
    constexpr double radius=2.65;
    constexpr double pi=3.1415926535897932384626433832795;
    wf.atoms.push_back({"C",6,0.0,0.0,0.0,6.0});
    for (std::size_t branch=0u;branch<3u;++branch) {
        const double angle=2.0*pi*static_cast<double>(branch)/3.0;
        wf.atoms.push_back({
            "C",6,radius*std::cos(angle),radius*std::sin(angle),0.0,6.0});
    }

    std::array<std::size_t,4> pz{};
    for (std::uint32_t atom=0u;atom<4u;++atom) {
        append_shell(wf,atom,0u,0u,1.2f);
        append_shell(wf,atom,1u,0u,0.8f);
        pz[atom]=wf.shells.back().basis_offset+2u;
    }
    for (std::uint32_t terminal=1u;terminal<4u;++terminal) {
        wf.bond_orders.push_back(
            {0u,terminal,1.10,cov::DataProvenance::Derived});
    }
    wf.bond_order_provenance=cov::DataProvenance::Derived;

    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t basis=0u;basis<n;++basis) {
        wf.ao_overlap[basis*n+basis]=1.0;
    }
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;

    const std::size_t spin_blocks=unrestricted?2u:1u;
    wf.orbitals.resize(spin_blocks*n);
    for (std::size_t spin_block=0u;spin_block<spin_blocks;++spin_block) {
        for (std::size_t basis=0u;basis<n;++basis) {
            auto& mo=wf.orbitals[spin_block*n+basis];
            mo.spin=spin_block==0u?cov::Spin::Alpha:cov::Spin::Beta;
            mo.coefficients.assign(n,0.0f);
            mo.coefficients[basis]=1.0f;
            mo.energy_hartree=-1.60+0.07*static_cast<double>(basis)+
                0.001*static_cast<double>(spin_block);
            const bool occupied_pi=
                basis==pz[0] || basis==pz[1];
            mo.occupation=occupied_pi
                ?static_cast<float>(unrestricted?1.0:2.0)
                :0.0f;
            mo.occupation_provenance=cov::DataProvenance::Producer;
        }
    }
    return wf;
}

const cov::DelocalisedPiAssignment* branched_assignment(
    const cov::Wavefunction& wf) {
    const auto found=std::find_if(
        wf.delocalised_pi_assignments.begin(),
        wf.delocalised_pi_assignments.end(),
        [](const auto& assignment) {
            return assignment.topology==
                cov::DelocalisedPiTopology::BranchedResonance &&
                assignment.atoms.size()==4u;
        });
    return found==wf.delocalised_pi_assignments.end()?nullptr:&*found;
}

std::vector<std::size_t> detected_members(const cov::Wavefunction& wf) {
    std::vector<std::size_t> result;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        const auto& chemistry=wf.orbitals[i].chemistry;
        if (!chemistry.delocalised_family_id.empty()) result.push_back(i);
    }
    return result;
}

bool validate_cp_family(
    const cov::Wavefunction& wf,const char* context,
    const std::vector<std::size_t>& expected_indices) {
    std::vector<std::uint32_t> expected_orbitals;
    expected_orbitals.reserve(expected_indices.size());
    for (const auto index:expected_indices) {
        expected_orbitals.push_back(static_cast<std::uint32_t>(index));
    }
    const auto detected=detected_members(wf);
    if (detected!=expected_indices) {
        std::cerr<<context<<": expected delocalised-pi members";
        for (const auto index:expected_indices) std::cerr<<' '<<(index+1u);
        std::cerr<<"; detected";
        for (const auto index:detected) std::cerr<<' '<<(index+1u);
        std::cerr<<'\n';
        return false;
    }

    if (wf.delocalised_pi_assignments.size()!=1u) {
        std::cerr<<context<<": expected one wavefunction-level delocalised family, got "
                 <<wf.delocalised_pi_assignments.size()<<'\n';
        return false;
    }
    const auto& assignment=wf.delocalised_pi_assignments.front();
    const std::vector<std::uint32_t> expected_atoms{0u,1u,2u,3u,4u};
    if (assignment.family_id.empty() || assignment.atoms!=expected_atoms ||
        assignment.orbitals!=expected_orbitals ||
        std::abs(assignment.electron_count-6.0)>1.0e-8 ||
        assignment.orientation_channels.size()!=1u ||
        !assignment.orientation_channels.front().cyclic ||
        !assignment.cyclic_topology ||
        assignment.provenance!=cov::DataProvenance::Derived) {
        std::cerr<<context<<": wavefunction-level Π^5_6 assignment is incomplete\n";
        return false;
    }

    const auto& reference=wf.orbitals[expected_indices.front()].chemistry;
    if (reference.delocalised_family_id.empty() ||
        !reference.multicentre_label.empty() ||
        reference.family_symbol!="pi" ||
        reference.delocalised_participating_atoms!=5u ||
        std::abs(reference.delocalised_participating_electrons-6.0)>1.0e-8 ||
        reference.delocalised_participating_atom_indices!=expected_atoms ||
        reference.delocalised_family_orbitals!=expected_orbitals) {
        std::cerr<<context<<": incomplete Π^5_6 family metadata\n";
        return false;
    }

    const std::string expected_label=
        "\xCE\xA0\xE2\x81\xB5\xE2\x82\x86"; // Π⁵₆
    for (const auto index:expected_indices) {
        const auto& chemistry=wf.orbitals[index].chemistry;
        if (chemistry.delocalised_family_id!=reference.delocalised_family_id ||
            chemistry.delocalised_family_orbitals!=expected_orbitals ||
            chemistry.delocalised_participating_atom_indices!=expected_atoms ||
            !std::isfinite(chemistry.delocalised_pi_weight) ||
            chemistry.delocalised_pi_weight<=0.0 || chemistry.delocalised_pi_weight>1.0) {
            std::cerr<<context<<": family membership differs at MO"
                     <<(index+1u)<<'\n';
            return false;
        }
        // Local bond-axis character and global p-family membership can differ.
        // The annotation must preserve membership without rewriting that channel.
        auto scoped=wf.orbitals[index];
        scoped.chemistry.channel.dominant=cov::OrbitalAngularFamily::Sigma;
        const auto annotation=cov::annotate_orbital(scoped);
        if (scoped.chemistry.channel.dominant!=cov::OrbitalAngularFamily::Sigma) return false;
        if (!annotation.delocalised_pi.available ||
            annotation.delocalised_pi.label!=expected_label ||
            annotation.delocalised_pi.participating_atoms!=5u ||
            std::abs(annotation.delocalised_pi.participating_electrons-6.0)>1.0e-8 ||
            annotation.delocalised_pi.orbital_indices!=expected_indices ||
            annotation.delocalised_pi.family_id!=reference.delocalised_family_id) {
            std::cerr<<context<<": Π⁵₆ diagram descriptor missing at MO"
                     <<(index+1u)<<'\n';
            return false;
        }
    }
    return true;
}

bool validate_optional_real_input(const std::filesystem::path& path) {
    try {
        const auto wf=cov::parse_wavefunction(path);
        // The real def2-SVPD D5h calculation resolves the complete p_z
        // manifold as A2'' + E1'' + E2'': raw MO16/17/18/24/25.  The previous
        // optional oracle incorrectly expected MO42/43 (E1' sigma/Rydberg).
        return validate_cp_family(
            wf,"real Cp input",{15u,16u,17u,23u,24u});
    } catch (const std::exception& error) {
        std::cerr<<"real Cp input failed: "<<error.what()<<'\n';
        return false;
    }
}

} // namespace

int main(const int argc,char** argv) {
    auto planar=make_planar_c5_case(false);
    cov::derive_orbital_chemistry(planar);
    if (!validate_cp_family(
            planar,"synthetic planar C5",
            {kSyntheticExpectedMembers.begin(),kSyntheticExpectedMembers.end()})) {
        return EXIT_FAILURE;
    }

    // A single eta ring is already haptic; classification must not require a
    // second sandwich ring.  The Fe centre intentionally has no AO basis in
    // this topology fixture, so the result can only come from the measured
    // multi-contact graph and not from metal-specific MO coefficients.
    auto half_sandwich=make_planar_c5_case(false);
    const std::uint32_t iron=static_cast<std::uint32_t>(
        half_sandwich.atoms.size());
    half_sandwich.atoms.push_back(
        {"Fe",26,0.0,0.0,3.20,26.0});
    for (std::uint32_t carbon=0u;carbon<5u;++carbon) {
        half_sandwich.bond_orders.push_back(
            {iron,carbon,0.08,cov::DataProvenance::Derived});
    }
    cov::derive_orbital_chemistry(half_sandwich);
    const auto haptic=std::find_if(
        half_sandwich.delocalised_pi_assignments.begin(),
        half_sandwich.delocalised_pi_assignments.end(),
        [](const auto& assignment) {
            return assignment.atoms.size()==5u &&
                assignment.topology==cov::DelocalisedPiTopology::HapticMetal;
        });
    if (haptic==half_sandwich.delocalised_pi_assignments.end() ||
        cov::preferred_compact_mo_diagram_mode(half_sandwich,true)!=
            cov::MODiagramMode::DelocalisedPiFamilyOnly) {
        std::cerr<<"single eta5 ring was not classified as haptic pi\n";
        return EXIT_FAILURE;
    }

    // Larger f-block haptic contacts must use the complete element-radius
    // table.  This geometry gives d(La--C)=2.70 A: it is inside the La+C
    // first-shell envelope but outside the old generic-0.85-A metal envelope.
    auto lanthanide_half_sandwich=make_planar_c5_case(false);
    const std::uint32_t lanthanum=static_cast<std::uint32_t>(
        lanthanide_half_sandwich.atoms.size());
    lanthanide_half_sandwich.atoms.push_back(
        {"La",57,0.0,0.0,4.58,57.0});
    for (std::uint32_t carbon=0u;carbon<5u;++carbon) {
        lanthanide_half_sandwich.bond_orders.push_back(
            {lanthanum,carbon,0.08,cov::DataProvenance::Derived});
    }
    cov::derive_orbital_chemistry(lanthanide_half_sandwich);
    if (std::none_of(
            lanthanide_half_sandwich.delocalised_pi_assignments.begin(),
            lanthanide_half_sandwich.delocalised_pi_assignments.end(),
            [](const auto& assignment) {
                return assignment.atoms.size()==5u &&
                    assignment.topology==
                        cov::DelocalisedPiTopology::HapticMetal;
            })) {
        std::cerr<<"physical La eta5 contact was rejected by haptic shell\n";
        return EXIT_FAILURE;
    }

    // A non-local through-bond Mayer term is not a haptic contact.  Keep the
    // same complete ring pi manifold and electronic values, but move Fe well
    // outside the broad first-shell envelope: topology must stay an ordinary
    // cyclic pi family rather than a metal-haptic one.
    auto remote_metal=make_planar_c5_case(false);
    const std::uint32_t remote_iron=static_cast<std::uint32_t>(
        remote_metal.atoms.size());
    remote_metal.atoms.push_back(
        {"Fe",26,0.0,0.0,12.0,26.0});
    for (std::uint32_t carbon=0u;carbon<5u;++carbon) {
        remote_metal.bond_orders.push_back(
            {remote_iron,carbon,0.08,cov::DataProvenance::Derived});
    }
    cov::derive_orbital_chemistry(remote_metal);
    if (std::any_of(
            remote_metal.delocalised_pi_assignments.begin(),
            remote_metal.delocalised_pi_assignments.end(),
            [](const auto& assignment) {
                return assignment.topology==
                    cov::DelocalisedPiTopology::HapticMetal;
            })) {
        std::cerr<<"remote metal coupling became a haptic pi contact\n";
        return EXIT_FAILURE;
    }

    auto rooted_collision=make_rooted_cycle_case(false,false);
    cov::derive_orbital_chemistry(rooted_collision);
    std::size_t six_member_cycles=0u;
    for (const auto& assignment:rooted_collision.delocalised_pi_assignments) {
        if (assignment.topology==
                cov::DelocalisedPiTopology::SymmetryDirectSum ||
            assignment.atoms.size()==12u) {
            std::cerr<<"non-isomorphic rooted rings became one symmetry direct sum\n";
            return EXIT_FAILURE;
        }
        if (assignment.topology==cov::DelocalisedPiTopology::Cycle &&
            assignment.atoms.size()==6u) {
            ++six_member_cycles;
        }
    }
    if (six_member_cycles!=2u) {
        std::cerr<<"rooted-cycle negative control lost independent pi rings: "
                 <<six_member_cycles<<"; assignments=";
        for (const auto& assignment:rooted_collision.delocalised_pi_assignments) {
            std::cerr<<' '<<static_cast<int>(assignment.topology)
                     <<'/'<<assignment.atoms.size();
        }
        std::cerr<<'\n';
        return EXIT_FAILURE;
    }

    auto equivalent_hub=make_rooted_cycle_case(true,false);
    cov::derive_orbital_chemistry(equivalent_hub);
    const auto direct_sum=std::find_if(
        equivalent_hub.delocalised_pi_assignments.begin(),
        equivalent_hub.delocalised_pi_assignments.end(),
        [](const auto& assignment) {
            return assignment.topology==
                    cov::DelocalisedPiTopology::SymmetryDirectSum &&
                assignment.atoms.size()==12u &&
                assignment.orientation_channels.size()==2u;
        });
    if (direct_sum==equivalent_hub.delocalised_pi_assignments.end()) {
        std::cerr<<"whole-hub-equivalent rooted rings lost their direct sum; assignments=";
        for (const auto& assignment:equivalent_hub.delocalised_pi_assignments) {
            std::cerr<<' '<<static_cast<int>(assignment.topology)
                     <<'/'<<assignment.atoms.size()<<'/'
                     <<assignment.orientation_channels.size();
        }
        std::cerr<<'\n';
        return EXIT_FAILURE;
    }

    auto asymmetric_hub=make_rooted_cycle_case(true,true);
    cov::derive_orbital_chemistry(asymmetric_hub);
    std::size_t asymmetric_cycles=0u;
    for (const auto& assignment:asymmetric_hub.delocalised_pi_assignments) {
        if (assignment.topology==
                cov::DelocalisedPiTopology::SymmetryDirectSum ||
            assignment.atoms.size()==12u) {
            std::cerr<<"asymmetric hub shell falsely related two pi rings\n";
            return EXIT_FAILURE;
        }
        if (assignment.topology==cov::DelocalisedPiTopology::Cycle &&
            assignment.atoms.size()==6u) {
            ++asymmetric_cycles;
        }
    }
    if (asymmetric_cycles!=2u) {
        std::cerr<<"asymmetric hub control lost its independent pi rings\n";
        return EXIT_FAILURE;
    }

    auto nonplanar=make_planar_c5_case(true);
    cov::derive_orbital_chemistry(nonplanar);
    // Local aligned fragments may remain valid after puckering; the forbidden
    // claim is one globally coherent five-centre cyclic pi family.  Requiring
    // zero local pi assignments would incorrectly reject butadiene-like
    // fragments inside a non-planar ring (the same boundary used by tub-COT).
    for (const auto& assignment:nonplanar.delocalised_pi_assignments) {
        if (assignment.cyclic_topology || assignment.atoms.size()==5u) {
            std::cerr<<"non-planar C5 negative control became one cyclic/full-ring pi family\n";
            return EXIT_FAILURE;
        }
    }

    auto open_shell=make_uhf_o2_case();
    cov::derive_orbital_chemistry(open_shell);
    if (open_shell.delocalised_pi_assignments.size()!=1u) {
        std::cerr<<"UHF O2 fixture did not produce one bundled orthogonal-p family\n";
        return EXIT_FAILURE;
    }
    const auto& open_assignment=open_shell.delocalised_pi_assignments.front();
    if (open_assignment.atoms!=std::vector<std::uint32_t>{0u,1u} ||
        open_assignment.orbitals.size()!=8u ||
        std::abs(open_assignment.electron_count-6.0)>1.0e-8 ||
        open_assignment.orientation_channels.size()!=2u ||
        open_assignment.cyclic_topology ||
        open_assignment.topology!=cov::DelocalisedPiTopology::MultiChannel ||
        !open_assignment.topology_graph.channel_ring_witnesses.empty() ||
        !open_assignment.topology_graph.channel_association_search_complete ||
        open_assignment.rationale.find("2 orientation channel(s)")==
            std::string::npos) {
        std::cerr<<"UHF O2 orthogonal-p subspace was lost or mislabelled as a ring union\n";
        return EXIT_FAILURE;
    }
    double direction_dot=0.0;
    for (std::size_t axis=0;axis<3u;++axis) {
        direction_dot+=open_assignment.orientation_channels[0].direction[axis]*
            open_assignment.orientation_channels[1].direction[axis];
    }
    if (std::abs(direction_dot)>1.0e-6) {
        std::cerr<<"UHF O2 orientation channels are not orthogonal\n";
        return EXIT_FAILURE;
    }

    auto restricted_branch=make_branched_p_case(false);
    auto unrestricted_branch=make_branched_p_case(true);
    cov::derive_orbital_chemistry(restricted_branch);
    cov::derive_orbital_chemistry(unrestricted_branch);
    const auto* restricted_assignment=branched_assignment(restricted_branch);
    const auto* unrestricted_assignment=branched_assignment(unrestricted_branch);
    if (restricted_assignment==nullptr || unrestricted_assignment==nullptr) {
        std::cerr<<"equivalent RHF/UHF four-centre p stars were not both "
                    "classified as branched resonance\n";
        return EXIT_FAILURE;
    }
    constexpr double expected_branch_occupation=2.0;
    if (std::abs(restricted_assignment->branch_centre_projected_occupation-
                 expected_branch_occupation)>1.0e-10 ||
        std::abs(unrestricted_assignment->branch_centre_projected_occupation-
                 expected_branch_occupation)>1.0e-10 ||
        std::abs(restricted_assignment->branch_centre_projected_occupation-
                 unrestricted_assignment->branch_centre_projected_occupation)>
            1.0e-10 ||
        std::abs(restricted_assignment->electron_count-4.0)>1.0e-10 ||
        std::abs(unrestricted_assignment->electron_count-4.0)>1.0e-10 ||
        restricted_assignment->orbitals.size()!=4u ||
        unrestricted_assignment->orbitals.size()!=8u) {
        std::cerr<<"RHF/UHF branched-centre projected occupations differ: "
                 <<restricted_assignment->branch_centre_projected_occupation
                 <<" versus "
                 <<unrestricted_assignment->branch_centre_projected_occupation
                 <<'\n';
        return EXIT_FAILURE;
    }

    if (argc>1 && !validate_optional_real_input(argv[1])) {
        return EXIT_FAILURE;
    }

    std::cout<<"delocalised pi smoke test passed\n";
    return EXIT_SUCCESS;
}
