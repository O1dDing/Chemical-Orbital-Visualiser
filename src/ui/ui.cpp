#include "cov/ui.hpp"
#include "cov/orbital_ui_text.hpp"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <initializer_list>
#include <string>

namespace cov::ui {
namespace {

constexpr std::size_t kTextCount = static_cast<std::size_t>(Text::Count);

struct LocalisedString {
    const char* en;
    const char* zh;
    const char* ja;
    const char* fr;
};

constexpr std::array<LocalisedString, kTextCount> kStrings{{
    {"Chemical Orbital Visualiser", "Chemical Orbital Visualiser", "Chemical Orbital Visualiser", "Chemical Orbital Visualiser"},
    {"Orbital energies, occupations and connections", "轨道能量、电子占据与轨道联系", "軌道エネルギー・占有数・軌道間のつながり", "Énergies, occupations et liens entre orbitales"},
    {"Language", "语言", "言語", "Langue"},
    {"File", "文件", "ファイル", "Fichier"},
    {"Molden file", "Molden 文件", "Molden ファイル", "Fichier Molden"},
    {"Load", "加载", "読み込む", "Charger"},
    {"Drop a .molden file anywhere on the viewport, or enter a path.", "可将 .molden 文件拖到视口任意位置，或直接输入路径。", ".molden ファイルをビューポートへドロップするか、パスを入力してください。", "Déposez un fichier .molden dans la vue ou saisissez son chemin."},
    {"Wavefunction", "波函数", "波動関数", "Fonction d’onde"},
    {"Atoms", "原子数", "原子数", "Atomes"},
    {"Shells", "壳层数", "シェル数", "Couches"},
    {"Basis functions", "基函数数", "基底関数数", "Fonctions de base"},
    {"Orbitals", "轨道数", "軌道数", "Orbitales"},
    {"Shell convention", "壳层约定", "シェル規約", "Convention des couches"},
    {"Charge / multiplicity", "电荷 / 多重度", "電荷 / 多重度", "Charge / multiplicité"},
    {"Alpha / beta electrons", "Alpha / Beta 电子数", "Alpha / Beta 電子数", "Électrons alpha / bêta"},
    {"SCF / stability", "SCF / 稳定性", "SCF / 安定性", "SCF / stabilité"},
    {"<S²> before / after", "<S²> 消除前 / 后", "<S²> 消去前 / 後", "<S²> avant / après"},
    {"converged", "已收敛", "収束", "convergé"},
    {"failed", "失败", "失敗", "échec"},
    {"stable", "稳定", "安定", "stable"},
    {"unstable", "不稳定", "不安定", "instable"},
    {"Frame continuity", "几何帧连续性", "フレーム連続性", "Continuité des géométries"},
    {"Atom mapping", "原子映射", "原子対応", "Correspondance des atomes"},
    {"Matched orbital groups", "已匹配轨道组", "対応した軌道グループ", "Groupes orbitaux appariés"},
    {"Unmatched previous / current", "未匹配（前帧 / 当前帧）", "未対応（前 / 現在）", "Non appariés précédent / actuel"},
    {"Composite matching", "复合子空间匹配", "複合部分空間の対応", "Appariement des sous-espaces"},
    {"Exact / not needed", "精确 / 无需优化", "厳密 / 不要", "Exact / non requis"},
    {"Conservative fallback", "保守回退", "保守的フォールバック", "Repli conservateur"},
    {"Compatible", "兼容", "互換", "Compatible"},
    {"Incompatible", "不兼容", "非互換", "Incompatible"},
    {"No previous frame", "尚无前一帧", "前のフレームなし", "Aucune géométrie précédente"},
    {"Orbital", "轨道", "軌道", "Orbitale"},
    {"Molden MO (1-based)", "Molden MO（从 1 开始）", "Molden MO（1 始まり）", "MO Molden (base 1)"},
    {"Internal index", "内部索引", "内部インデックス", "Indice interne"},
    {"Energy", "能量", "エネルギー", "Énergie"},
    {"Occupation", "占据数", "占有数", "Occupation"},
    {"Spin", "自旋", "スピン", "Spin"},
    {"Symmetry", "对称性", "対称性", "Symétrie"},
    {"Rendering", "渲染", "レンダリング", "Rendu"},
    {"Isovalue", "等值面", "等値面", "Isovaleur"},
    {"Grid", "网格", "グリッド", "Grille"},
    {"Recompute grid", "重新计算网格", "グリッドを再計算", "Recalculer la grille"},
    {"Reset camera", "重置相机", "カメラをリセット", "Réinitialiser la caméra"},
    {"Performance", "性能", "パフォーマンス", "Performances"},
    {"Compute device", "计算设备", "計算デバイス", "Périphérique de calcul"},
    {"Last grid calculation", "最近网格计算", "直近のグリッド計算", "Dernier calcul de grille"},
    {"Grid ready", "网格就绪", "グリッド準備完了", "Grille prête"},
    {"Left-drag to orbit · mouse wheel to zoom", "按住鼠标左键旋转 · 滚轮缩放", "左ドラッグで回転 · ホイールでズーム", "Glisser gauche : rotation · molette : zoom"},
    {"Changing isovalue updates the display without recalculating the orbital grid.", "调整等值面会更新显示，无需重新计算轨道网格。", "等値面の変更は表示を更新し、軌道グリッドを再計算しません。", "Changer l’isovaleur actualise l’affichage sans recalculer la grille orbitale."},
    {"Experimental MVP · scientific validation in progress", "实验性 MVP · 科学数值验证仍在进行", "実験的 MVP · 科学的検証を継続中", "MVP expérimental · validation scientifique en cours"},
    {"Ready", "就绪", "準備完了", "Prêt"},
    {"Parsing", "正在解析", "解析中", "Analyse"},
    {"Loaded", "已加载", "読み込み完了", "Chargé"},
    {"Grid updated", "网格已更新", "グリッド更新完了", "Grille mise à jour"},
    {"Error", "错误", "エラー", "Erreur"},
    {"(none)", "（无）", "（なし）", "(aucun)"},
    {"Alpha", "Alpha", "Alpha", "Alpha"},
    {"Beta", "Beta", "Beta", "Beta"},
    {"UI fonts", "界面字体", "UI フォント", "Polices UI"},

    {"Open File…", "打开文件…", "ファイルを開く…", "Ouvrir un fichier…"},
    {"Recent files", "最近文件", "最近のファイル", "Fichiers récents"},
    {"Current file", "当前文件", "現在のファイル", "Fichier actuel"},
    {"Orbital browser", "轨道浏览器", "軌道ブラウザ", "Explorateur d’orbitales"},
    {"Search", "搜索", "検索", "Rechercher"},
    {"Filter", "筛选", "フィルター", "Filtre"},
    {"Auto · reasonable", "自动 · 合理范围", "自動 · 妥当範囲", "Auto · plage raisonnable"},
    {"All", "全部", "すべて", "Toutes"},
    {"Occupied", "已占据", "占有", "Occupées"},
    {"Virtual", "虚轨道", "仮想", "Virtuelles"},
    {"Core", "内层", "内殻", "Cœur"},
    {"Valence", "价层", "価電子", "Valence"},
    {"Virtual window", "虚轨道窗口", "仮想軌道ウィンドウ", "Fenêtre virtuelle"},
    {"Degeneracy tolerance", "简并阈值", "縮退判定しきい値", "Tolérance de dégénérescence"},
    {"Grouped labels", "分组标签", "グループ表示", "Étiquettes groupées"},
    {"Raw numbering", "原始编号", "元の番号", "Numérotation brute"},
    {"Degenerate set", "简并组", "縮退組", "Groupe dégénéré"},
    {"Energy unit", "能量单位", "エネルギー単位", "Unité d’énergie"},
    {"HOMO", "HOMO", "HOMO", "HOMO"},
    {"LUMO", "LUMO", "LUMO", "LUMO"},
    {"HOMO-1", "HOMO-1", "HOMO-1", "HOMO-1"},
    {"LUMO+1", "LUMO+1", "LUMO+1", "LUMO+1"},
    {"Valence MO diagram", "价电子层 MO 图", "価電子層 MO 図", "Diagramme MO de valence"},
    {"Valence diagram span", "价电子 MO 图范围", "価電子 MO 図の表示範囲", "Étendue du diagramme MO de valence"},
    {"Generate MO diagram", "生成 MO 图", "MO 図を生成", "Générer le diagramme MO"},
    {"Export diagram + metadata", "导出图与元数据", "図とメタデータを書き出す", "Exporter diagramme + métadonnées"},
    {"Exported", "已导出", "書き出し完了", "Exporté"},
    {"Export failed", "导出失败", "書き出し失敗", "Échec de l’export"},
    {"Molecule style", "分子样式", "分子表示", "Style moléculaire"},
    {"Enhanced ball-and-stick", "增强球棍模型", "強調ボール＆スティック", "Boules et bâtonnets renforcés"},
    {"Stick + delocalisation", "棍 + 离域虚线", "結合 + 非局在化破線", "Bâtons + délocalisation"},
    {"Atom size", "原子大小", "原子サイズ", "Taille des atomes"},
    {"Bond size", "键粗细", "結合の太さ", "Épaisseur des liaisons"},
    {"Molecule opacity", "分子透明度", "分子の不透明度", "Opacité de la molécule"},
    {"Orbital opacity", "轨道透明度", "軌道の不透明度", "Opacité de l’orbitale"},
    {"Show hydrogens", "显示氢原子", "水素を表示", "Afficher les hydrogènes"},
    {"Show coordination contacts", "显示配位连接", "配位結合を表示", "Afficher les contacts de coordination"},
    {"Show multicentre support", "显示多中心连接", "多中心支持を表示", "Afficher le support multicentrique"},
    {"Show polyhedral cage support", "显示多面体笼骨架支撑", "多面体ケージ骨格を表示", "Afficher le support de cage polyédrique"},
    {"Show weak interactions", "显示弱相互作用", "弱い相互作用を表示", "Afficher les interactions faibles"},
    {"Hydrogen-bond, non-covalent and ionic contacts only; ambiguous contacts stay hidden.", "仅显示氢键、非共价和离子接触；歧义接触仍保持隐藏。", "水素結合・非共有結合・イオン接触のみ。曖昧な接触は表示しません。", "Contacts hydrogène, non covalents et ioniques uniquement ; les contacts ambigus restent masqués."},
    {"Dashed bonds use a conservative delocalisation heuristic.", "虚线键使用保守的离域启发式判断。", "破線結合は保守的な非局在化ヒューリスティックです。", "Les liaisons en pointillés utilisent une heuristique prudente de délocalisation."},
    {"Central valence layout", "中央价电子层布局", "中央価電子層レイアウト", "Disposition centrale de valence"},
    {"Valence-grouped levels", "价电子层分组能级", "価電子層のグループ準位", "Niveaux groupés de valence"},
    {"Valence MO diagram", "价电子层 MO 图", "価電子層 MO 図", "Diagramme MO de valence"},
    {"Machine metadata", "机读元数据", "機械可読メタデータ", "Métadonnées machine"},
    {"Visible orbitals", "可见轨道", "表示軌道", "Orbitales visibles"},
    {"Raw MO", "原始 MO", "元の MO", "MO brute"},
    {"Region", "区域", "領域", "Région"},
    {"Core", "内层", "内殻", "Cœur"},
    {"Valence", "价层", "価電子", "Valence"},
    {"Virtual", "虚轨道", "仮想", "Virtuelle"},
    {"Reasonable energy window", "合理能量窗口", "妥当なエネルギー範囲", "Fenêtre d’énergie raisonnable"},
    {"Native Open File is unavailable on this platform.", "当前平台不支持原生“打开文件”。", "このプラットフォームではネイティブのファイル選択を利用できません。", "La boîte de dialogue native n’est pas disponible sur cette plateforme."},
    {"Copy metadata", "复制元数据", "メタデータをコピー", "Copier les métadonnées"},
    {"No orbitals", "无轨道", "軌道がありません", "Aucune orbitale"},
    {"Adaptive nonlinear energy scale (log-gap v3)", "自适应非线性能量轴（log-gap v3）", "適応型非線形エネルギー軸（log-gap v3）", "Échelle d’énergie non linéaire adaptative (log-gap v3)"},
    {"Energy scale", "能量轴", "エネルギー軸", "Échelle d’énergie"},
    {"Linear", "线性", "線形", "Linéaire"},
    {"Adaptive nonlinear", "自适应非线性", "適応型非線形", "Non linéaire adaptative"},
    {"Orbital family", "轨道类型", "軌道タイプ", "Famille orbitale"},
    {"Bonding class", "成键类别", "結合分類", "Classe de liaison"},
    {"Exact energy", "精确能量", "正確なエネルギー", "Énergie exacte"},
    {"Multicentre bond", "多中心键", "多中心結合", "Liaison multicentrique"},
    {"Delocalised π system", "离域 π 体系", "非局在化 π 系", "Système π délocalisé"},
    {"Classification source", "分类来源", "分類の出典", "Source de classification"},
    {"Confidence", "置信度", "信頼度", "Confiance"},
    {"Degenerate members", "简并成员", "縮退メンバー", "Membres dégénérés"},
}};

// These strings deliberately mirror text that is rendered directly by
// main.cpp, ui_text_dispatch.cpp and orbital_ui_dispatch.cpp instead of being
// routed through kStrings.  Building the CJK range from kStrings alone used to
// omit characters in labels such as “轨道材质” and “柔和自动打光”, which made
// Dear ImGui display '?' even though the operating-system CJK font was loaded.
constexpr const char* kSupplementalChinese =
    "轨道材质 标准 玻璃 表面模式 实体 线框 实体 + 线框 柔和自动打光 "
    "波函数文件（FCHK 优先；Molden 兼容） "
    "可拖入 .fchk/.fch/.chk 或兼容的 .molden 文件，也可直接输入路径。 "
    "源文件 MO（从 1 开始） 源 MO "
    "所选 MO 化学性质 化学价层轨道组 是 否 价层 AO 组成 "
    "原子对相互作用 轨道类型 成键性质 多中心 / 离域族 "
    "成员 MO 参与原子 参与电子 "
    "给体 / 受体方向 分析方法 MO 贡献 "
    "MO 贡献来自重叠布居；Mayer 为总密度原子对指数。 "
    "UND / 最小价层参考之外 CUDA 设备 CUDA设备";

constexpr const char* kSupplementalJapanese =
    "軌道マテリアル 標準 ガラス 表示モード ソリッド ワイヤー "
    "ソリッド + ワイヤー ソフト自動照明 "
    "波動関数ファイル（FCHK 優先・Molden 互換） "
    ".fchk/.fch/.chk または互換 .molden ファイルをドロップするか、"
    "パスを入力してください。 入力 MO（1 始まり） 入力 MO "
    "選択 MO の化学的性質 化学原子価軌道群 はい いいえ "
    "原子価 AO 組成 原子対相互作用 軌道型 結合性 "
    "構成 MO 参加原子 参加電子 "
    "多中心 / 非局在化族 供与体 / 受容体 解析法 MO 寄与 "
    "MO 寄与は重なり密度由来、Mayer は全密度の原子対指数。 "
    "UND / 最小原子価参照外 CUDA デバイス";

// Keep all symbols produced by MO labels/annotations in the primary font.
// Π⁵₆ is included as an exact sequence as well as through the complete digit
// sets, which protects both the large-pi family label and future N-centre
// families from atlas-range regressions.
constexpr const char* kScientificGlyphs =
    "● · – — − ± × → ← ↔ ↑ ↓ "
    "α β σ π δ φ Σ Π Δ Φ Γ Π⁵₆ "
    "⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ ⁺ ⁻ "
    "₀ ₁ ₂ ₃ ₄ ₅ ₆ ₇ ₈ ₉ ₊ ₋ ′ ″";

std::string g_font_status = "Dear ImGui default";

std::string first_existing(std::initializer_list<const char*> candidates) {
    std::error_code ec;
    for (const char* candidate : candidates) {
        if (candidate && std::filesystem::exists(candidate, ec) && !ec) return candidate;
        ec.clear();
    }
    return {};
}

const char* file_name_or_default(const std::string& path, const char* fallback) {
    if (path.empty()) return fallback;
    static thread_local std::string name;
    name = std::filesystem::path(path).filename().string();
    return name.c_str();
}

void merge_font_if_available(const std::string& path,
                             const float pixel_size,
                             const ImWchar* ranges,
                             bool& loaded) {
    if (path.empty()) return;
    ImFontConfig cfg{};
    cfg.MergeMode = true;
    cfg.PixelSnapH = true;
    cfg.OversampleH = 1;
    cfg.OversampleV = 1;
    if (ImGui::GetIO().Fonts->AddFontFromFileTTF(path.c_str(), pixel_size, &cfg, ranges)) loaded = true;
}

const char* localised(const LocalisedString& value, const Language language) noexcept {
    switch (language) {
        case Language::ChineseSimplified: return value.zh;
        case Language::Japanese: return value.ja;
        case Language::French: return value.fr;
        default: return value.en;
    }
}

} // namespace

const char* tr(const Text key, const Language language) noexcept {
    const auto k = static_cast<std::size_t>(key);
    if (k >= kTextCount) return "";
    return localised(kStrings[k], language);
}

const char* language_name(const Language language) noexcept {
    switch (language) {
        case Language::ChineseSimplified: return "简体中文";
        case Language::Japanese: return "日本語";
        case Language::French: return "Français";
        default: return "English";
    }
}

const char* supplemental_glyph_seed(const Language language) noexcept {
    static const std::string english=orbital_ui_glyph_seed(Language::English);
    static const std::string chinese=std::string(kSupplementalChinese)+" "+
        orbital_ui_glyph_seed(Language::ChineseSimplified);
    static const std::string japanese=std::string(kSupplementalJapanese)+" "+
        orbital_ui_glyph_seed(Language::Japanese);
    static const std::string french=orbital_ui_glyph_seed(Language::French);
    switch (language) {
        case Language::ChineseSimplified: return chinese.c_str();
        case Language::Japanese: return japanese.c_str();
        case Language::French: return french.c_str();
        default: return english.c_str();
    }
}

const char* scientific_glyph_seed() noexcept { return kScientificGlyphs; }

void apply_theme(const float scale) {
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowPadding = ImVec2(14.0f, 14.0f);
    style.FramePadding = ImVec2(10.0f, 6.0f);
    style.CellPadding = ImVec2(8.0f, 5.0f);
    style.ItemSpacing = ImVec2(9.0f, 8.0f);
    style.ItemInnerSpacing = ImVec2(7.0f, 5.0f);
    style.IndentSpacing = 18.0f;
    style.ScrollbarSize = 11.0f;
    style.GrabMinSize = 10.0f;
    style.WindowRounding = 10.0f;
    style.ChildRounding = 8.0f;
    style.FrameRounding = 7.0f;
    style.PopupRounding = 8.0f;
    style.ScrollbarRounding = 10.0f;
    style.GrabRounding = 6.0f;
    style.TabRounding = 6.0f;
    style.WindowBorderSize = 1.0f;
    style.ChildBorderSize = 1.0f;
    style.FrameBorderSize = 0.0f;
    style.PopupBorderSize = 1.0f;

    ImVec4* c = style.Colors;
    c[ImGuiCol_Text] = ImVec4(0.90f, 0.93f, 0.97f, 1.00f);
    c[ImGuiCol_TextDisabled] = ImVec4(0.53f, 0.59f, 0.67f, 1.00f);
    c[ImGuiCol_WindowBg] = ImVec4(0.055f, 0.073f, 0.102f, 0.965f);
    c[ImGuiCol_ChildBg] = ImVec4(0.082f, 0.108f, 0.145f, 0.96f);
    c[ImGuiCol_PopupBg] = ImVec4(0.070f, 0.091f, 0.122f, 0.985f);
    c[ImGuiCol_Border] = ImVec4(0.18f, 0.23f, 0.30f, 0.95f);
    c[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_FrameBg] = ImVec4(0.112f, 0.145f, 0.195f, 1.00f);
    c[ImGuiCol_FrameBgHovered] = ImVec4(0.145f, 0.195f, 0.270f, 1.00f);
    c[ImGuiCol_FrameBgActive] = ImVec4(0.175f, 0.230f, 0.315f, 1.00f);
    c[ImGuiCol_TitleBg] = c[ImGuiCol_WindowBg];
    c[ImGuiCol_TitleBgActive] = c[ImGuiCol_WindowBg];
    c[ImGuiCol_TitleBgCollapsed] = c[ImGuiCol_WindowBg];
    c[ImGuiCol_MenuBarBg] = ImVec4(0.070f, 0.091f, 0.122f, 1.00f);
    c[ImGuiCol_ScrollbarBg] = ImVec4(0.045f, 0.060f, 0.083f, 0.80f);
    c[ImGuiCol_ScrollbarGrab] = ImVec4(0.21f, 0.27f, 0.35f, 0.95f);
    c[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.30f, 0.39f, 0.50f, 1.00f);
    c[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.36f, 0.46f, 0.59f, 1.00f);
    c[ImGuiCol_CheckMark] = ImVec4(0.36f, 0.56f, 0.98f, 1.00f);
    c[ImGuiCol_SliderGrab] = ImVec4(0.36f, 0.56f, 0.98f, 0.88f);
    c[ImGuiCol_SliderGrabActive] = ImVec4(0.48f, 0.67f, 1.00f, 1.00f);
    c[ImGuiCol_Button] = ImVec4(0.18f, 0.32f, 0.58f, 0.82f);
    c[ImGuiCol_ButtonHovered] = ImVec4(0.26f, 0.43f, 0.78f, 1.00f);
    c[ImGuiCol_ButtonActive] = ImVec4(0.22f, 0.37f, 0.68f, 1.00f);
    c[ImGuiCol_Header] = ImVec4(0.18f, 0.32f, 0.58f, 0.62f);
    c[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.43f, 0.78f, 0.86f);
    c[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.43f, 0.78f, 1.00f);
    c[ImGuiCol_Separator] = ImVec4(0.18f, 0.23f, 0.30f, 0.80f);
    c[ImGuiCol_SeparatorHovered] = ImVec4(0.36f, 0.56f, 0.98f, 0.75f);
    c[ImGuiCol_SeparatorActive] = ImVec4(0.36f, 0.56f, 0.98f, 1.00f);
    c[ImGuiCol_ResizeGrip] = ImVec4(0.36f, 0.56f, 0.98f, 0.20f);
    c[ImGuiCol_ResizeGripHovered] = ImVec4(0.36f, 0.56f, 0.98f, 0.65f);
    c[ImGuiCol_ResizeGripActive] = ImVec4(0.36f, 0.56f, 0.98f, 0.90f);
    c[ImGuiCol_Tab] = ImVec4(0.10f, 0.14f, 0.20f, 1.00f);
    c[ImGuiCol_TabHovered] = ImVec4(0.26f, 0.43f, 0.78f, 0.85f);
    c[ImGuiCol_TabActive] = ImVec4(0.18f, 0.32f, 0.58f, 0.95f);
    c[ImGuiCol_TabUnfocused] = ImVec4(0.08f, 0.11f, 0.15f, 1.00f);
    c[ImGuiCol_TabUnfocusedActive] = ImVec4(0.13f, 0.20f, 0.31f, 1.00f);
    c[ImGuiCol_TableHeaderBg] = ImVec4(0.10f, 0.14f, 0.19f, 1.00f);
    c[ImGuiCol_TableBorderStrong] = ImVec4(0.18f, 0.23f, 0.30f, 0.85f);
    c[ImGuiCol_TableBorderLight] = ImVec4(0.15f, 0.19f, 0.25f, 0.65f);
    c[ImGuiCol_TextSelectedBg] = ImVec4(0.36f, 0.56f, 0.98f, 0.32f);
    c[ImGuiCol_DragDropTarget] = ImVec4(0.31f, 0.82f, 0.65f, 0.95f);
    c[ImGuiCol_NavHighlight] = ImVec4(0.36f, 0.56f, 0.98f, 0.90f);
    c[ImGuiCol_ModalWindowDimBg] = ImVec4(0.02f, 0.03f, 0.05f, 0.70f);
    style.ScaleAllSizes(std::clamp(scale, 0.85f, 2.0f));
}

bool configure_fonts(const float pixel_size) {
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();
    io.Fonts->Flags |= ImFontAtlasFlags_NoPowerOfTwoHeight;
#ifdef _WIN32
    // Windows 11: prefer Segoe UI Variable Text for cleaner high-DPI UI; keep
    // classic Segoe UI as the Windows 10 fallback. CJK glyphs are merged from
    // the OS-provided UI faces, so the downloadable cov.exe ships no font files.
    const std::string base = first_existing({"C:/Windows/Fonts/SegUIVar.ttf", "C:/Windows/Fonts/segoeui.ttf", "C:/Windows/Fonts/arial.ttf"});
    const std::string chinese = first_existing({"C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/msyh.ttf", "C:/Windows/Fonts/simhei.ttf"});
    const std::string japanese = first_existing({"C:/Windows/Fonts/YuGothM.ttc", "C:/Windows/Fonts/YuGothR.ttc", "C:/Windows/Fonts/meiryo.ttc", "C:/Windows/Fonts/msgothic.ttc"});
#elif defined(__APPLE__)
    const std::string base = first_existing({"/System/Library/Fonts/SFNS.ttf", "/System/Library/Fonts/Supplemental/Arial.ttf"});
    const std::string chinese = first_existing({"/System/Library/Fonts/PingFang.ttc", "/Library/Fonts/NotoSansCJK-Regular.ttc"});
    const std::string japanese = first_existing({"/Library/Fonts/NotoSansCJK-Regular.ttc", "/System/Library/Fonts/AppleGothic.ttf"});
#else
    const std::string base = first_existing({"/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf"});
    const std::string chinese = first_existing({"/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", "/usr/share/fonts/opentype/noto/NotoSansCJKsc-Regular.otf"});
    const std::string japanese = first_existing({"/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", "/usr/share/fonts/opentype/noto/NotoSansCJKjp-Regular.otf"});
#endif

    ImFontGlyphRangesBuilder latin_builder;
    latin_builder.AddRanges(io.Fonts->GetGlyphRangesDefault());
    for (const auto& row : kStrings) { latin_builder.AddText(row.en); latin_builder.AddText(row.fr); }
    latin_builder.AddText(language_name(Language::French));
    latin_builder.AddText(supplemental_glyph_seed(Language::English));
    latin_builder.AddText(supplemental_glyph_seed(Language::French));
    latin_builder.AddText(scientific_glyph_seed());
    ImVector<ImWchar> latin_ranges;
    latin_builder.BuildRanges(&latin_ranges);

    ImFont* primary = nullptr;
    if (!base.empty()) {
        ImFontConfig cfg{};
        cfg.OversampleH = 3;
        cfg.OversampleV = 1;
        primary = io.Fonts->AddFontFromFileTTF(base.c_str(), pixel_size, &cfg, latin_ranges.Data);
    }
    if (!primary) primary = io.Fonts->AddFontDefault();

    ImFontGlyphRangesBuilder zh_builder;
    ImFontGlyphRangesBuilder ja_builder;
    for (const auto& row : kStrings) { zh_builder.AddText(row.zh); ja_builder.AddText(row.ja); }
    zh_builder.AddText(language_name(Language::ChineseSimplified));
    ja_builder.AddText(language_name(Language::Japanese));
    zh_builder.AddText(supplemental_glyph_seed(Language::ChineseSimplified));
    ja_builder.AddText(supplemental_glyph_seed(Language::Japanese));
    zh_builder.AddText(scientific_glyph_seed());
    ja_builder.AddText(scientific_glyph_seed());
    ImVector<ImWchar> zh_ranges;
    ImVector<ImWchar> ja_ranges;
    zh_builder.BuildRanges(&zh_ranges);
    ja_builder.BuildRanges(&ja_ranges);

    bool zh_loaded = false;
    bool ja_loaded = false;
    merge_font_if_available(chinese, pixel_size, zh_ranges.Data, zh_loaded);
    merge_font_if_available(japanese, pixel_size, ja_ranges.Data, ja_loaded);
    const bool built = io.Fonts->Build();
    g_font_status = file_name_or_default(base, "ImGui default");
    g_font_status += " + ";
    g_font_status += zh_loaded ? file_name_or_default(chinese, "CJK") : "ZH fallback missing";
    g_font_status += " + ";
    g_font_status += ja_loaded ? file_name_or_default(japanese, "CJK") : "JA fallback missing";
    return built;
}

const char* font_status() noexcept { return g_font_status.c_str(); }

void section_title(const char* label) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.58f, 0.68f, 0.82f, 1.0f));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
}

void begin_card(const char* id, const float height) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(12.0f, 11.0f));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.082f, 0.108f, 0.145f, 0.94f));
    ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.18f, 0.23f, 0.30f, 0.90f));
    ImGui::BeginChild(id, ImVec2(0.0f, height), true, ImGuiWindowFlags_None);
    ImGui::PushTextWrapPos(0.0f);
}

void end_card() {
    ImGui::PopTextWrapPos();
    ImGui::EndChild();
    ImGui::PopStyleColor(2);
    ImGui::PopStyleVar();
}

void status_badge(const char* label, const Tone tone) {
    ImVec4 colour{};
    switch (tone) {
        case Tone::Success: colour = ImVec4(0.31f, 0.82f, 0.65f, 1.0f); break;
        case Tone::Danger: colour = ImVec4(0.94f, 0.42f, 0.47f, 1.0f); break;
        case Tone::Accent: colour = ImVec4(0.42f, 0.64f, 1.00f, 1.0f); break;
        default: colour = ImVec4(0.58f, 0.68f, 0.82f, 1.0f); break;
    }
    ImGui::TextColored(colour, "●  %s", label);
}

} // namespace cov::ui
