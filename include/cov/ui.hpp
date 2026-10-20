#pragma once

namespace cov::ui {

enum class Language {
    English = 0,
    ChineseSimplified,
    Japanese,
    French,
    Count,
};

enum class Text {
    AppTitle = 0,
    Tagline,
    LanguageLabel,
    FileSection,
    MoldenPath,
    Load,
    IdleHint,
    WavefunctionSection,
    Atoms,
    Shells,
    BasisFunctions,
    Orbitals,
    ShellConvention,
    ChargeMultiplicity,
    AlphaBetaElectrons,
    SCFStability,
    SpinSquared,
    Converged,
    Failed,
    Stable,
    Unstable,
    FrameTracking,
    AtomMappingCompatibility,
    MatchedSubspaces,
    UnmatchedSubspaces,
    TrackingOptimisation,
    ExactOrNotNeeded,
    ConservativeFallback,
    Compatible,
    Incompatible,
    NoPreviousFrame,
    OrbitalSection,
    MoldenMO,
    InternalIndex,
    Energy,
    Occupation,
    Spin,
    Symmetry,
    RenderingSection,
    Isovalue,
    Grid,
    RecomputeGrid,
    ResetCamera,
    PerformanceSection,
    CUDADevice,
    LastKernel,
    GPUResident,
    InteractionHint,
    IsovalueHint,
    Ready,
    Parsing,
    Loaded,
    GridUpdated,
    Error,
    NoneValue,
    Alpha,
    Beta,
    FontStatus,

    OpenFile,
    RecentFiles,
    CurrentFile,
    OrbitalBrowser,
    Search,
    Filter,
    FilterAuto,
    FilterAll,
    FilterOccupied,
    FilterVirtual,
    FilterCore,
    FilterValence,
    HighVirtualWindow,
    DegeneracyTolerance,
    GroupedLabels,
    RawNumbering,
    DegenerateSet,
    EnergyUnit,
    HOMO,
    LUMO,
    HOMOMinus1,
    LUMOPlus1,
    EnergyDiagram,
    AroundSelected,
    GenerateMODiagram,
    ExportBundle,
    Exported,
    ExportFailed,
    MoleculeStyle,
    MediumBallStick,
    StickDelocalisation,
    AtomSize,
    BondSize,
    MoleculeOpacity,
    OrbitalOpacity,
    ShowHydrogens,
    ShowCoordinationContacts,
    ShowMulticentreSupport,
    ShowPolyhedralCageSupport,
    ShowWeakInteractions,
    WeakInteractionsHint,
    DelocalisationHeuristic,
    SALCUnavailable,
    SymmetryGrouped,
    SimpleDiagram,
    MachineMetadata,
    VisibleOrbitals,
    RawMO,
    Region,
    Core,
    Valence,
    Virtual,
    ReasonableWindow,
    OpenDialogUnsupported,
    CopyMetadata,
    NoOrbitals,
    NonlinearEnergyScale,
    EnergyScale,
    LinearEnergyScale,
    NonlinearFocus,
    OrbitalFamily,
    BondingClassLabel,
    ExactEnergy,
    MulticentreBond,
    DelocalisedPiSystem,
    ClassificationSource,
    DegenerateMembers,
    Count,
};

enum class Tone {
    Neutral,
    Accent,
    Success,
    Danger,
};

[[nodiscard]] const char* tr(Text key, Language language) noexcept;
[[nodiscard]] const char* language_name(Language language) noexcept;
// Text rendered outside the main localisation table (for example the orbital
// material and selected-MO chemistry panels).  Font construction consumes
// this seed so those glyphs cannot silently fall back to '?'.
[[nodiscard]] const char* supplemental_glyph_seed(Language language) noexcept;
[[nodiscard]] const char* scientific_glyph_seed() noexcept;

void apply_theme(float scale = 1.0f);
bool configure_fonts(float pixel_size = 17.0f);
[[nodiscard]] const char* font_status() noexcept;
[[nodiscard]] const char* font_status(Language language);

void section_title(const char* label);
void begin_card(const char* id, float height);
void end_card();
void status_badge(const char* label, Tone tone);

} // namespace cov::ui
