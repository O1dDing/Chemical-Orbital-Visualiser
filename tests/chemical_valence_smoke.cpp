#include "cov/mo_diagram.hpp"
#include "cov/orbital_chemistry.hpp"
#include "cov/orbital_view.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

cov::Primitive primitive(const float exponent) {
    cov::Primitive value;
    value.exponent=exponent;
    value.coefficient=1.0f;
    return value;
}

cov::Wavefunction make_h3_reference_case() {
    cov::Wavefunction wf;
    wf.atoms={
        {"H",1,-1.0,0.0,0.0,1.0},
        {"H",1, 0.5,0.866025403784,0.0,1.0},
        {"H",1, 0.5,-0.866025403784,0.0,1.0},
    };

    for (std::uint32_t atom=0;atom<3u;++atom) {
        const std::uint32_t primitive_offset=
            static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back(primitive(1.2f));
        wf.shells.push_back({atom,primitive_offset,1u,wf.basis_count,0u,0u});
        wf.basis_count+=1u;

        const std::uint32_t diffuse_offset=
            static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back(primitive(0.010f));
        wf.shells.push_back({atom,diffuse_offset,1u,wf.basis_count,0u,0u});
        wf.basis_count+=1u;

        const std::uint32_t p_offset=
            static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back(primitive(0.050f));
        wf.shells.push_back({atom,p_offset,1u,wf.basis_count,1u,0u});
        wf.basis_count+=3u;
    }

    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) wf.ao_overlap[i*n+i]=1.0;
    const std::array<std::size_t,3> main{0u,5u,10u};
    constexpr double s=0.20;
    for (std::size_t i=0;i<3u;++i) {
        for (std::size_t j=i+1u;j<3u;++j) {
            wf.ao_overlap[main[i]*n+main[j]]=s;
            wf.ao_overlap[main[j]*n+main[i]]=s;
        }
    }

    cov::MolecularOrbital a1;
    a1.energy_hartree=-1.022932;
    a1.occupation=2.0f;
    a1.symmetry="A1'";
    a1.coefficients.assign(n,0.0f);
    const double ca=1.0/std::sqrt(3.0+6.0*s);
    for (const auto index:main) a1.coefficients[index]=static_cast<float>(ca);

    cov::MolecularOrbital e1;
    e1.energy_hartree=-0.34;
    e1.occupation=0.0f;
    e1.symmetry="E'";
    e1.coefficients.assign(n,0.0f);
    const double ce1=1.0/std::sqrt(6.0*(1.0-s));
    e1.coefficients[main[0]]=static_cast<float>(2.0*ce1);
    e1.coefficients[main[1]]=static_cast<float>(-ce1);
    e1.coefficients[main[2]]=static_cast<float>(-ce1);

    cov::MolecularOrbital e2;
    e2.energy_hartree=-0.34;
    e2.occupation=0.0f;
    e2.symmetry="E'";
    e2.coefficients.assign(n,0.0f);
    const double ce2=1.0/std::sqrt(2.0*(1.0-s));
    e2.coefficients[main[1]]=static_cast<float>(ce2);
    e2.coefficients[main[2]]=static_cast<float>(-ce2);

    wf.orbitals={a1,e1,e2};
    for (std::size_t basis=0;basis<n;++basis) {
        if (basis==main[0] || basis==main[1] || basis==main[2]) continue;
        cov::MolecularOrbital extra;
        extra.energy_hartree=0.1+0.1*static_cast<double>(wf.orbitals.size());
        extra.occupation=0.0f;
        extra.coefficients.assign(n,0.0f);
        extra.coefficients[basis]=1.0f;
        wf.orbitals.push_back(std::move(extra));
    }

    wf.bond_orders={
        {0u,1u,0.5,cov::DataProvenance::Derived},
        {1u,2u,0.5,cov::DataProvenance::Derived},
        {0u,2u,0.5,cov::DataProvenance::Derived},
    };
    return wf;
}

cov::Wavefunction make_cation_visibility_case() {
    cov::Wavefunction wf;
    wf.atoms={{"H",1,0.0,0.0,0.0,1.0}};
    wf.charge=1;
    wf.multiplicity=2u;
    wf.charge_provenance=cov::DataProvenance::Producer;
    wf.multiplicity_provenance=cov::DataProvenance::Producer;

    const auto orbital=[](const double energy,
                          const float occupation,
                          const double deep_core,
                          const double valence,
                          const double unresolved) {
        cov::MolecularOrbital mo;
        mo.energy_hartree=energy;
        mo.occupation=occupation;
        mo.chemistry.available=true;
        mo.chemistry.deep_core_weight=deep_core;
        mo.chemistry.valence_weight=valence;
        mo.chemistry.unresolved_weight=unresolved;
        return mo;
    };
    wf.orbitals={
        orbital(-5.0,2.0f,0.92,0.02,0.06), // explicit deep core
        orbital(-0.7,2.0f,0.05,0.03,0.92), // diffuse/unresolved occupied
        orbital(-0.4,1.0f,0.00,0.02,0.98), // SOMO outside minimal rank
        orbital(-0.2,0.0f,0.00,0.03,0.97), // cation vacancy / LUMO
        orbital( 0.8,0.0f,0.00,0.01,0.99), // unrelated high virtual
    };
    return wf;
}

cov::Wavefunction make_local_linear_three_centre_case(const bool negative_edge) {
    cov::Wavefunction wf;
    // The fourth atom is a spectator: the 3c4e detector must reason over a
    // local linear unit rather than requiring the whole molecule to be
    // triatomic.
    wf.atoms={
        {"H",1,-1.5,0.0,0.0,1.0},
        {"He",2,0.0,0.0,0.0,2.0},
        {"H",1,1.5,0.0,0.0,1.0},
        {"H",1,0.0,10.0,0.0,1.0},
    };
    for (std::uint32_t atom=0;atom<4u;++atom) {
        const auto primitive_offset=static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back(primitive(1.0f));
        wf.shells.push_back({atom,primitive_offset,1u,wf.basis_count,0u,0u});
        ++wf.basis_count;
    }
    wf.ao_overlap.assign(16u,0.0);
    for (std::size_t i=0;i<4u;++i) wf.ao_overlap[i*4u+i]=1.0;

    const auto orbital=[](const double energy,const float occupation,
                          const std::array<float,4>& coefficients) {
        cov::MolecularOrbital mo;
        mo.energy_hartree=energy;
        mo.occupation=occupation;
        mo.coefficients.assign(coefficients.begin(),coefficients.end());
        return mo;
    };
    constexpr float inv_two=0.5f;
    constexpr float inv_sqrt_two=0.7071067811865475f;
    wf.orbitals={
        orbital(-0.50,2.0f,{inv_two,inv_sqrt_two,inv_two,0.0f}),
        orbital(-0.20,2.0f,{inv_sqrt_two,0.0f,-inv_sqrt_two,0.0f}),
        orbital( 0.20,0.0f,{inv_two,-inv_sqrt_two,inv_two,0.0f}),
        orbital( 1.00,0.0f,{0.0f,0.0f,0.0f,1.0f}),
    };
    wf.bond_orders={
        {0u,1u,0.45,cov::DataProvenance::Derived},
        {1u,2u,negative_edge?-0.45:0.45,cov::DataProvenance::Derived},
        {0u,2u,0.00,cov::DataProvenance::Derived},
    };
    wf.bond_order_provenance=cov::DataProvenance::Derived;
    wf.multicentre_candidates.push_back({
        0u,{1u,0u,2u},{0.50,0.25,0.25},2.0,
        cov::DataProvenance::Derived,
    });
    return wf;
}

} // namespace

int main() {
    const cov::OrbitalBondingDistribution unavailable;
    if (unavailable.status!=cov::ChemistryStatus::Unavailable ||
        unavailable.undetermined!=1.0) {
        std::cerr<<"Default unavailable bonding evidence must remain UND=1\n";
        return EXIT_FAILURE;
    }
    auto wf=make_h3_reference_case();
    cov::derive_orbital_chemistry(wf);

    for (const auto& mo:wf.orbitals) {
        const auto& distribution=mo.chemistry.bonding;
        if (!mo.chemistry.available ||
            distribution.status==cov::ChemistryStatus::NotApplicable) continue;
        const double total=distribution.bonding+distribution.antibonding+
            distribution.nonbonding+distribution.undetermined;
        if (!std::isfinite(total) || std::abs(total-1.0)>1.0e-10 ||
            distribution.undetermined<0.0 || distribution.undetermined>1.0) {
            std::cerr<<"Aggregate bonding fractions must conserve unit weight\n";
            return EXIT_FAILURE;
        }
    }

    std::vector<std::size_t> selected;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        if (wf.orbitals[i].chemistry.valence_manifold) selected.push_back(i);
    }
    if (selected!=std::vector<std::size_t>{0u,1u,2u}) {
        std::cerr<<"H3 chemical-valence manifold must be MO1 + E' pair only\n";
        return EXIT_FAILURE;
    }
    if (wf.orbitals[0].chemistry.channel.dominant!=
            cov::OrbitalAngularFamily::Sigma ||
        wf.orbitals[0].chemistry.channel.status!=
            cov::ChemistryStatus::Determined) {
        std::cerr<<"H3 occupied valence MO was not classified as sigma\n";
        return EXIT_FAILURE;
    }
    if (wf.orbitals[0].chemistry.bonding.dominant!=
            cov::OrbitalBondingRole::Bonding) {
        std::cerr<<"H3 occupied valence MO was not classified as bonding\n";
        return EXIT_FAILURE;
    }
    if (wf.orbitals[0].chemistry.multicentre_label!="3c2e" ||
        wf.orbitals[0].chemistry.multicentre_participating_atoms!=3u ||
        std::abs(wf.orbitals[0].chemistry.
                 multicentre_participating_electrons-2.0)>1.0e-8) {
        std::cerr<<"H3 3c2e family was not attached to the valence manifold\n";
        return EXIT_FAILURE;
    }

    cov::OrbitalFilterSettings filter;
    filter.mode=cov::OrbitalFilterMode::AutoReasonable;
    const auto frontier=cov::find_frontier_orbitals(wf.orbitals);
    const auto visible=cov::visible_orbital_indices(wf.orbitals,frontier,filter);
    if (visible!=selected) {
        std::cerr<<"Auto filter did not use the chemical-valence manifold\n";
        return EXIT_FAILURE;
    }

    cov::MODiagramOptions options;
    options.selected_index=0;
    const auto diagram=cov::build_mo_diagram_data(wf,options);
    if (diagram.selection.included_indices.size()!=3u ||
        diagram.levels.size()!=2u ||
        diagram.levels[1].member_indices.size()!=2u) {
        std::cerr<<"MO diagram did not preserve three H3 valence MOs in two symmetry rows\n";
        return EXIT_FAILURE;
    }
    if (!diagram.levels[0].chemistry.available ||
        diagram.levels[0].chemistry.ao_contributions.size()!=3u) {
        std::cerr<<"MO diagram did not carry chemistry/AO composition\n";
        return EXIT_FAILURE;
    }

    auto local_3c4e=make_local_linear_three_centre_case(false);
    cov::derive_orbital_chemistry(local_3c4e);
    const bool local_assignment=std::any_of(
        local_3c4e.multicentre_assignments.begin(),
        local_3c4e.multicentre_assignments.end(),[](const auto& assignment){
            return assignment.kind==cov::MulticentreKind::ThreeCentreFourElectron &&
                assignment.atoms==std::vector<std::uint32_t>{0u,1u,2u};
        });
    if (!local_assignment) {
        std::cerr<<"spectator-bearing local linear 3c4e unit was not recognised\n";
        return EXIT_FAILURE;
    }
    auto remote_terminal_density=make_local_linear_three_centre_case(false);
    remote_terminal_density.bond_orders[2].mayer_order=0.20;
    cov::derive_orbital_chemistry(remote_terminal_density);
    if (remote_terminal_density.multicentre_assignments.empty()) {
        std::cerr<<"remote terminal Mayer density suppressed a linear 3c4e unit\n";
        return EXIT_FAILURE;
    }
    auto negative_3c=make_local_linear_three_centre_case(true);
    cov::derive_orbital_chemistry(negative_3c);
    if (!negative_3c.multicentre_assignments.empty()) {
        std::cerr<<"negative Mayer edge became structural 3c support\n";
        return EXIT_FAILURE;
    }

    // A density-only three-orbital coincidence is provisional, not a 3c
    // assignment.  Once the same three atoms are bent/branched rather than a
    // supported linear or hydrogen-bridge framework, orbital chemistry must
    // discard the provisional label instead of exposing false 3c2e metadata.
    auto branched_provisional=make_local_linear_three_centre_case(false);
    branched_provisional.atoms[2].x=0.0;
    branched_provisional.atoms[2].y=1.5;
    cov::MulticentreAssignment provisional;
    provisional.kind=cov::MulticentreKind::ThreeCentreTwoElectron;
    provisional.atoms={0u,1u,2u};
    provisional.orbitals={0u,1u,2u};
    provisional.electron_count=2.0;
    provisional.rationale="density-only provisional assignment";
    provisional.provenance=cov::DataProvenance::Derived;
    branched_provisional.multicentre_assignments.push_back(provisional);
    cov::derive_orbital_chemistry(branched_provisional);
    if (!branched_provisional.multicentre_assignments.empty()) {
        std::cerr<<"branched density-only fragment retained a false 3c label\n";
        return EXIT_FAILURE;
    }
    auto producer_multicentre=make_local_linear_three_centre_case(false);
    producer_multicentre.atoms[2].x=0.0;
    producer_multicentre.atoms[2].y=1.5;
    provisional.provenance=cov::DataProvenance::Producer;
    producer_multicentre.multicentre_assignments.push_back(provisional);
    cov::derive_orbital_chemistry(producer_multicentre);
    if (producer_multicentre.multicentre_assignments.size()!=1u ||
        producer_multicentre.orbitals[0].chemistry.multicentre_label!="3c2e" ||
        !producer_multicentre.orbitals[0].chemistry.multicentre_assignment_available ||
        !cov::annotate_orbital(producer_multicentre.orbitals[0]).multicentre.available) {
        std::cerr<<"producer multicentre assignment was not preserved\n";
        return EXIT_FAILURE;
    }

    // A generic population/geometry fallback has candidate text but no
    // assignment-to-member evidence. It must not become an established bond.
    auto candidate=producer_multicentre.orbitals[0];
    candidate.chemistry.multicentre_assignment_available=false;
    if (cov::annotate_orbital(candidate).multicentre.available ||
        candidate.chemistry.multicentre_label!="3c2e") {
        std::cerr<<"candidate multicentre text was promoted or destroyed\n";
        return EXIT_FAILURE;
    }
    for(const auto& assignment:local_3c4e.multicentre_assignments)
        for(const auto index:assignment.orbitals)
            if(index<local_3c4e.orbitals.size() &&
               !cov::annotate_orbital(local_3c4e.orbitals[index]).multicentre.available) {
                std::cerr<<"qualified assignment member lost its multicentre descriptor\n";
                return EXIT_FAILURE;
            }

    auto ion=make_cation_visibility_case();
    cov::OrbitalFilterSettings ion_filter;
    ion_filter.mode=cov::OrbitalFilterMode::AutoReasonable;
    const auto ion_frontier=cov::find_frontier_orbitals(ion.orbitals);
    const auto ion_visible=cov::visible_orbital_indices(
        ion.orbitals,ion_frontier,ion_filter);
    if (ion_visible!=std::vector<std::size_t>{1u,2u}) {
        std::cerr<<"non-core occupied/SOMO visibility safeguard regression\n";
        return EXIT_FAILURE;
    }
    if (!cov::confidently_deep_core_orbital(ion.orbitals[0],ion_filter) ||
        cov::confidently_deep_core_orbital(ion.orbitals[1],ion_filter)) {
        std::cerr<<"deep-core confidence gate regression\n";
        return EXIT_FAILURE;
    }

    const auto ion_metadata=cov::build_orbital_metadata(
        ion,1u,{},ion_filter);
    if (ion_metadata.size()!=5u || ion_metadata[0].visible ||
        !ion_metadata[1].visible || !ion_metadata[2].visible ||
        !ion_metadata[3].visible || ion_metadata[4].visible) {
        std::cerr<<"positive-ion frontier-vacancy visibility regression\n";
        return EXIT_FAILURE;
    }

    cov::MODiagramOptions ion_options;
    ion_options.selected_index=1u;
    const auto ion_diagram=cov::build_mo_diagram_data(ion,ion_options);
    if (ion_diagram.selection.included_indices!=
            std::vector<std::size_t>{1u,2u,3u} ||
        ion_diagram.levels.size()!=3u) {
        std::cerr<<"central diagram silently dropped occupied/SOMO or cation vacancy\n";
        return EXIT_FAILURE;
    }

    ion_filter.mode=cov::OrbitalFilterMode::All;
    const auto all_indices=cov::visible_orbital_indices(
        ion.orbitals,ion_frontier,ion_filter);
    if (all_indices!=std::vector<std::size_t>{0u,1u,2u,3u,4u}) {
        std::cerr<<"full canonical MO accessibility regression\n";
        return EXIT_FAILURE;
    }

    std::cout<<"chemical valence/ion visibility smoke test passed\n";
    return EXIT_SUCCESS;
}
