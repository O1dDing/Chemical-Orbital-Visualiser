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
