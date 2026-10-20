#include "cov/ui.hpp"
#include "cov/orbital_ui_text.hpp"

namespace cov::ui {

const char* tr_legacy(Text key, Language language) noexcept;

const char* tr(Text key, Language language) noexcept {
    switch (key) {
        case Text::MoldenPath:
            return orbital_tr(OrbitalText::CalculationInput,language);
        case Text::IdleHint:
            return orbital_tr(OrbitalText::SupportedInputs,language);
        case Text::MoldenMO:
            return orbital_tr(OrbitalText::SourceMONumber,language);
        case Text::RawMO:
            switch (language) {
                case Language::ChineseSimplified: return "源 MO";
                case Language::Japanese: return "入力 MO";
                case Language::French: return "MO source";
                default: return "Source MO";
            }
        case Text::ClassificationSource:
            return orbital_tr(OrbitalText::LabelSource,language);
        default:
            return tr_legacy(key,language);
    }
}

} // namespace cov::ui
