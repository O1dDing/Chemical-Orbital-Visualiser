/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#pragma once

#include "cov/model.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/ui.hpp"

#include <cstddef>
#include <string>
#include <string_view>

namespace cov::ui {

// Localisation keys used by the FCHK-first provenance strip, ligand-field
// tooltips, delocalised-pi annotations and selected-MO chemistry panel.  Keep
// scientific identifiers (FCHK/MO/AO, point groups, energy units, Mayer,
// N/A and sigma/pi/delta/phi symbols) outside this table.
enum class OrbitalText : std::size_t {
    ProvenanceSymmetry = 0,
    ProvenanceOccupation,
    Density,
    Overlap,
    BondOrder,
    PointGroup,
    GaussianEnrichmentAttached,
    Producer,
    Derived,
    Unavailable,
    UnknownSource,

    SelectedMOChemistry,
    ChemicalValenceManifold,
    Yes,
    No,
    ValenceAOComposition,
    AtomPairInteractions,
    OrbitalFamily,
    BondingRole,
    MulticentreFamily,
    DelocalisedPiFamily,
    MemberMOs,
    ParticipatingAtoms,
    ParticipatingElectrons,
    DonorAcceptorDirection,
    AnalysisMethod,
    MOContribution,
    MOContributionExplanation,
    OutsideMinimalReference,
    AOMetricUnavailable,
    DelocalisedPi,
    Mixed,
    MixedUnd,
    Bonding,
    Antibonding,
    Nonbonding,

    OrientationChannels,
    Topology,
    TopologyPath,
    TopologyCycle,
    TopologyBranchedResonance,
    TopologySpiro,
    TopologyHapticMetal,
    TopologySymmetryDirectSum,
    ChannelAtoms,
    Coherence,
    HideIntermediateFrameworkMOs,
    HideIntermediateFrameworkMOsTooltip,

    LocalLigandField,
    FirstShell,
    CoordinationGeometry,
    CoordinationNumber,
    AngularRMS,
    RadialVariation,
    LocalMolecularGeometries,
    Atom,
    LocalGeometry,
    NotCentredOnThisMO,
    GroupOccupation,
    MetalSPD,
    LigandP,
    SigmaPiChannel,
    MetalLigandOverlap,
    Selection,
    RecoveredFromRawMOBlock,
    WeakFieldTreatment,
    ApproximatelyNonbonding,
    PiInteraction,
    PiSplitting,
    Interaction,
    Splitting,
    SplittingHartree,

    AnnotationDirect,
    AnnotationParsedLabel,
    AnnotationDerived,
    AnnotationHeuristic,
    PiDonorSplitting,
    PiAcceptorSplitting,
    PiCoupledSplitting,
    PiWeakNearNonbonding,

    DiagramSummary,
    ModeValenceCentral,
    ModeDelocalisedPi,
    ModeMulticentre,
    VisibleLevels,
    Occupied,
    Virtual,
    Hidden,
    LocalField,
    Geometry,
    PiPairs,
    ProtectedOverflow,
    SpinPairs,
    UnmatchedVisible,
    RawRecoveredGroups,
    CompactActiveSpaceUnavailable,

    TopologyMultiChannel,
    RingTopologyEvidence,
    ConnectivityModel,
    MayerDistanceModel,
    CovalentDistanceModel,
    RingUnionCentre,
    RingPath,
    AdditionalRingConnection,
    NoChannelRingWitness,
    RingGraphModelExplanation,

    SymmetryOriginalLabel,
    SymmetryExplanation,
    SymmetryLocalTag,
    SymmetryCandidateTag,
    SymmetrySourceTag,
    SymmetryMolecularTag,
    SymmetryMixedTag,
    SymmetryProjectionOrigin,
    SymmetryDimensionOrigin,
    SymmetryPiOrigin,
    SymmetrySpinOrigin,
    SymmetryMolecularOrigin,
    SymmetryTargetMembers,
    SymmetrySourceMembers,
    SymmetryCentreCoverage,
    SymmetryShellPurity,
    SymmetryTargetCoverage,
    SymmetryProjectionMeaning,
    SymmetryGroupUnresolved,
    SymmetryAxes,
    SymmetryDecomposition,
    SymmetryRepresentedRank,
    SymmetryNumericalResidual,
    SymmetryNoSingleLabel,
    CrystalFieldGap,
    CrystalFieldGaps,
    CataloguePrior,
    CatalogueSigmaOnly,
    CatalogueDonor,
    CatalogueAcceptor,
    CatalogueAmbiguous,
    CatalogueEvidenceRelation,
    CatalogueConsistent,
    CatalogueContradicted,
    GapLocalScopeExplanation,
    OrbitalDetails,
    CloseOrbitalDetails,
    OrbitalDetailsHint,
    OrbitalDetailsOutsideDiagram,
    OrbitalDetailsDataScope,
    LevelGroupContainsMOs,
    GroupRepresentativeData,
    MetalLigandGroupAnalysis,
    MetalLigandNotApplicable,
    MetalLigandUnavailable,
    MetalLigandSigmaPiChannel,
    RelatedLocalOrbitals,
    ProjectionRelations,
    CanonicalProjection,
    LocalizedE2,
    AllRelatedOrbitals,
    PiCoupling,
    PiCoupledDimensions,
    CrossFockRange,
    PartitionWeights,
    PiNoRelation,
    OrbitalComposition,
    NAOProjection,
    AOReferenceProjection,
    RetainedNAONorm,
    ReconstructionResidual,
    NotAnalysed,
    InsufficientData,
    IncompatibleData,
    UnsupportedData,
    NotApplicableData,
    BelowReportThreshold,
    NotDetermined,
    BondingMixing,
    AntibondingMixing,
    MixedCoupling,
    SymmetricCoupling,
    LigandToMetal,
    MetalToLigand,
    FullPointGroup,
    AbelianSubgroup,
    AvailableData,
    ComputationFailed,
    TotalDensity,
    SpinDensity,
    SharedSpinDescription,
    MissingSpinDescription,
    MetricDataUnavailable,
    IndependentOverlap,
    LabelSource,
    CalculationInput,
    SupportedInputs,
    SourceMONumber,
    Count,
};

[[nodiscard]] const char* orbital_tr(OrbitalText key, Language language) noexcept;
[[nodiscard]] const char* localised_wavefunction_source(
    WavefunctionSource source, Language language) noexcept;
[[nodiscard]] const char* localised_data_provenance(
    DataProvenance provenance, Language language) noexcept;
[[nodiscard]] const char* localised_annotation_source(
    AnnotationSource source, Language language) noexcept;
[[nodiscard]] const char* localised_bonding_class(
    BondingClass value, Language language) noexcept;
[[nodiscard]] const char* localised_orbital_bonding_role(
    OrbitalBondingRole value, Language language) noexcept;
[[nodiscard]] const char* localised_pi_interaction_kind(
    PiInteractionKind kind, Language language) noexcept;

[[nodiscard]] std::string localised_geometry_name(
    std::string_view geometry_id,
    std::string_view fallback_name,
    Language language);
[[nodiscard]] std::string localised_chemistry_method(
    std::string_view method, Language language);
[[nodiscard]] std::string localised_chemistry_note(
    std::string_view note, Language language);
[[nodiscard]] std::string localised_diagram_selection_summary(
    const MODiagramData& data, Language language);

// All glyphs rendered by this localisation layer. Font construction consumes
// the seed so CJK and accented French strings cannot silently degrade to '?'.
[[nodiscard]] const char* orbital_ui_glyph_seed(Language language) noexcept;

} // namespace cov::ui
