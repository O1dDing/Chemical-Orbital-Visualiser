#include "cov/ui.hpp"
#include "cov/nbo_ui.hpp"
#include "cov/nbo_aomo_ui.hpp"
#include "cov/validation.hpp"
#include "cov/orbital_ui_text.hpp"
#include "cov/orbital_inspection_ui.hpp"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <initializer_list>
#include <string>

#ifdef __APPLE__
#include <CoreGraphics/CoreGraphics.h>
#include <CoreText/CoreText.h>
#include <cmath>
#include <cstdio>
#include <memory>
#include <type_traits>
#include <vector>
#endif

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
    {"Orbital energies, occupations and connections", "轨道能量、电子占据与轨道联系", "軌道エネルギー・占有数・軌道間の関係", "Énergies, occupations et relations orbitalaires"},
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
    {"Simplified matching", "已使用简化匹配", "簡略化した対応付け", "Appariement simplifié"},
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
    {"Auto", "自动", "自動", "Auto"},
    {"All", "全部", "すべて", "Toutes"},
    {"Occupied", "已占据", "占有", "Occupées"},
    {"Virtual", "虚轨道", "仮想", "Virtuelles"},
    {"Core", "内层", "内殻", "Cœur"},
    {"Valence", "价层", "価電子", "Valence"},
    {"Virtual window", "虚轨道窗口", "仮想軌道ウィンドウ", "Fenêtre virtuelle"},
    {"Degeneracy tolerance", "简并阈值", "縮退判定しきい値", "Tolérance de dégénérescence"},
    {"Grouped labels", "分组标签", "グループ表示", "Étiquettes groupées"},
    {"Raw numbering", "原始编号", "元の番号", "Numérotation brute"},
    {"Energy-group size", "能级组成员数", "エネルギー群の成分数", "Taille du groupe de niveaux"},
    {"Energy unit", "能量单位", "エネルギー単位", "Unité d’énergie"},
    {"HOMO", "HOMO", "HOMO", "HOMO"},
    {"LUMO", "LUMO", "LUMO", "LUMO"},
    {"HOMO-1", "HOMO-1", "HOMO-1", "HOMO-1"},
    {"LUMO+1", "LUMO+1", "LUMO+1", "LUMO+1"},
    {"Valence MO diagram", "价电子层 MO 图", "価電子層 MO 図", "Diagramme MO de valence"},
    {"Valence diagram span", "价电子 MO 图范围", "価電子 MO 図の表示範囲", "Étendue du diagramme MO de valence"},
    {"Generate MO diagram", "生成 MO 图", "MO 図を生成", "Générer le diagramme MO"},
    {"Export images", "导出图片", "画像を書き出す", "Exporter les images"},
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
    {"Delocalised bonds: dashed lines.", "离域键：虚线。", "非局在化結合：破線。", "Liaisons délocalisées : pointillés."},
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
    {"Nonlinear energy axis", "非线性能量轴", "非線形エネルギー軸", "Axe d’énergie non linéaire"},
    {"Energy scale", "能量轴", "エネルギー軸", "Échelle d’énergie"},
    {"Linear", "线性", "線形", "Linéaire"},
    {"Nonlinear", "非线性", "非線形", "Non linéaire"},
    {"Orbital family", "轨道类型", "軌道タイプ", "Famille orbitale"},
    {"Bonding class", "成键类别", "結合分類", "Classe de liaison"},
    {"Exact energy", "精确能量", "正確なエネルギー", "Énergie exacte"},
    {"Multicentre bond", "多中心键", "多中心結合", "Liaison multicentrique"},
    {"Delocalised π system", "离域 π 体系", "非局在化 π 系", "Système π délocalisé"},
    {"Classification source", "分类来源", "分類の出典", "Source de classification"},
    {"Energy-group members", "能级组成员", "エネルギー群の成分", "Membres du groupe de niveaux"},
}};

// These strings deliberately mirror text that is rendered directly by
// main.cpp, ui_text_dispatch.cpp and orbital_ui_dispatch.cpp instead of being
// routed through kStrings.  Building the CJK range from kStrings alone used to
// omit characters in labels such as “轨道材质” and “柔和自动打光”, which made
// Dear ImGui display '?' even though the operating-system CJK font was loaded.
constexpr const char* kSupplementalChinese =
    "无法定位载入数据当前计算不匹配找到多个文件读取失败详细错误字体缺少中文日文字形分子轨道能级已识别的连接整数键级未确定 "
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
    "互換性のないデータです "
    "データが見つかりません現在の計算と一致しません複数読み込めませんファイル読込失敗エラー詳細中国語日本語のフォントがありません分子軌道のエネルギー整数結合次数が未確定 "
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

// Exact non-ASCII characters used by the four-language integration controls
// in main.cpp, nbo_aomo_ui.cpp and nbo_ui.cpp. The range builder deduplicates
// these source-derived characters instead of loading a full CJK range.
constexpr const char* kIntegrationChineseGlyphs =
    "²–…→−、。一三上下不与且两严个中为主义互些交仅仍代件会位低体作使保候值元先入全关内再出击分划"
    "则删别到前力加动勾化印原及取变叠只可右号合同后和器图在场均型域基堡声壳处复多央失始子存完定实导尾"
    "局层居展属左已布带幅平并度开弱归当待微德恢情成或截所手打拖择指按据接控描放数整文断新方旋无明是显"
    "暗有未杂权来构析架查标核格检正此段母比没注测浏消清源滚灰点片独瓣用电留白百的相看真着确示离称移空"
    "立符等简算类系素紫累红级纳线组结绘续维绿缩缺置而联能自色节荷蓝藏行表要视览角计证该详说请调负轨轮"
    "输这连述适选透逐道部配里重量金键间阈降除随隐集零面项题验骨高（），：；";
constexpr const char* kIntegrationJapaneseGlyphs =
    "²–…→−、。あいえかがきくげこさしすせただちつてでとなにねのはびぶべまみむもらりるれわをんァア"
    "イクグスセタッテデトドピフブプホメラリルロン・ー一上下不中主乗二互交付以位体作係保個候値元入全典"
    "内出分列別利削割力加動化区印原厳去可右号各合含和図在基場変外大央子字存定実密小局展属左布幅底度従"
    "復微成所手択拠拡持指振描損操数整新明書未析査格検構機欠次正残殻注消淡混済準濃灰点独用画異白的目直"
    "相着確示移積空立符等算系素紫累細結続緑線縮群能色荷行表補複見規覧親角解計証詳認説調負赤軌追退透造"
    "連運道選部配重量金開間関除際隠集零電青非面項骨（）：；";
constexpr const char* kIntegrationLatinGlyphs = "²Éèéê–’";

// Direct labels introduced in orbital_ui_v2.cpp and ui_text_dispatch.cpp
// bypass the orbital browser's own localisation table.
constexpr const char* kOrbitalDiagramChineseGlyphs =
    "MO 图设置 正则 MO 能量参考图 "
    "波函数文件或计算目录（FCHK 优先；自动关联 NBO） "
    "可同时拖入波函数与 NBO 文件，或拖入计算目录；也可直接输入文件或目录路径。兼容 .fchk/.fch/.chk 和 .molden。 "
    "源文件 MO（从 1 开始） 源 MO";
constexpr const char* kOrbitalDiagramJapaneseGlyphs =
    "MO 図の設定 正準 MO エネルギー参照図 "
    "波動関数ファイルまたは計算フォルダー（FCHK 優先・NBO 自動関連付け） "
    "波動関数と NBO ファイルをまとめて、または計算フォルダーをドロップできます。ファイルやフォルダーのパス入力も可能です。.fchk/.fch/.chk・.molden に対応。 "
    "入力 MO（1 始まり） 入力 MO";
constexpr const char* kOrbitalDiagramLatinGlyphs =
    "Réglages du diagramme OM Référence énergétique des OM canoniques "
    "Fichier de fonction d’onde ou dossier de calcul (FCHK prioritaire ; association NBO automatique) "
    "Déposez ensemble les fichiers de fonction d’onde et NBO, ou un dossier de calcul ; vous pouvez aussi saisir leur chemin. Formats .fchk/.fch/.chk et .molden compatibles.";

// Keep all symbols produced by MO labels/annotations in the primary font.
// Π⁵₆ is included as an exact sequence as well as through the complete digit
// sets, which protects both the large-pi family label and future N-centre
// families from atlas-range regressions.
constexpr const char* kScientificGlyphs =
    "● · – — − ± × ≈ → ← ↔ ↑ ↓ "
    "α β σ π δ φ η Σ Π Δ Φ Γ Π⁵₆ ∞ "
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

#ifdef __APPLE__
struct CFReleaseOwned {
    void operator()(const void* value) const { if (value) CFRelease(value); }
};
template<class T> using CFOwned = std::unique_ptr<std::remove_pointer_t<T>, CFReleaseOwned>;

std::string core_text_string(CFStringRef value) {
    if (!value) return {};
    const auto capacity = CFStringGetMaximumSizeForEncoding(CFStringGetLength(value), kCFStringEncodingUTF8) + 1;
    std::string text(static_cast<std::size_t>(capacity), '\0');
    if (!CFStringGetCString(value, text.data(), capacity, kCFStringEncodingUTF8)) return {};
    text.resize(std::char_traits<char>::length(text.c_str()));
    return text;
}

struct SystemGlyphBitmap {
    int rectangle = -1;
    int width = 0;
    int height = 0;
    std::vector<unsigned char> pixels;
};

// CoreText selects the actual installed font and collection face. Rendering
// through the system also handles Apple's variable/nonstandard glyph tables,
// which a readable PingFang.ttc plus stb_truetype cannot reliably represent.
bool add_system_glyphs(ImFont* primary, const float pixel_size, const ImWchar* ranges,
                       CFStringRef language, CFStringRef preferred_name,
                       ImFontGlyphRangesBuilder& pending,
                       std::vector<SystemGlyphBitmap>& bitmaps, std::string& source_name) {
    auto* atlas = primary->ContainerAtlas;
    CFOwned<CTFontRef> preferred(CTFontCreateWithName(preferred_name, pixel_size, nullptr));
    if (!preferred) return false;
    bool complete = true;
    bool reported = false;
    for (const ImWchar* range = ranges; range[0]; range += 2) {
        for (unsigned int codepoint = range[0]; codepoint <= range[1]; ++codepoint) {
            if (codepoint <= 0x20 || primary->FindGlyphNoFallback(static_cast<ImWchar>(codepoint)) || pending.GetBit(codepoint)) continue;
            const UniChar character = static_cast<UniChar>(codepoint);
            CFOwned<CFStringRef> text(CFStringCreateWithCharacters(kCFAllocatorDefault, &character, 1));
            CFOwned<CTFontRef> font(CTFontCreateForStringWithLanguage(preferred.get(), text.get(), CFRangeMake(0, 1), language));
            CGGlyph glyph = 0;
            CFOwned<CFStringRef> postscript(font ? CTFontCopyPostScriptName(font.get()) : nullptr);
            const auto name = core_text_string(postscript.get());
            if (!font || name.find("LastResort") != std::string::npos ||
                !CTFontGetGlyphsForCharacters(font.get(), &character, &glyph, 1) || glyph == 0) {
                std::fprintf(stderr, "System font missing %s U+%04X\n", core_text_string(language).c_str(), codepoint);
                complete = false;
                continue;
            }
            CGRect bounds{};
            CGSize advance{};
            CTFontGetBoundingRectsForGlyphs(font.get(), kCTFontOrientationHorizontal, &glyph, &bounds, 1);
            CTFontGetAdvancesForGlyphs(font.get(), kCTFontOrientationHorizontal, &glyph, &advance, 1);
            const bool whitespace = codepoint == 0xA0 || (codepoint >= 0x2000 && codepoint <= 0x200A) ||
                                    codepoint == 0x202F || codepoint == 0x205F || codepoint == 0x3000;
            const int left = static_cast<int>(std::floor(CGRectGetMinX(bounds))) - 1;
            const int bottom = static_cast<int>(std::floor(CGRectGetMinY(bounds))) - 1;
            const int top = static_cast<int>(std::ceil(CGRectGetMaxY(bounds))) + 1;
            const int right = static_cast<int>(std::ceil(CGRectGetMaxX(bounds))) + 1;
            SystemGlyphBitmap bitmap;
            bitmap.width = std::max(1, right - left);
            bitmap.height = std::max(1, top - bottom);
            bitmap.pixels.resize(static_cast<std::size_t>(bitmap.width) * bitmap.height, 0);
            CFOwned<CGColorSpaceRef> gray(CGColorSpaceCreateDeviceGray());
            CFOwned<CGContextRef> context(CGBitmapContextCreate(bitmap.pixels.data(), bitmap.width, bitmap.height,
                                                              8, bitmap.width, gray.get(), kCGImageAlphaNone));
            if (!context) { complete = false; continue; }
            CGContextSetAllowsAntialiasing(context.get(), true);
            CGContextSetShouldAntialias(context.get(), true);
            CGContextSetShouldSmoothFonts(context.get(), false);
            CGContextSetTextDrawingMode(context.get(), kCGTextFill);
            CGContextSetGrayFillColor(context.get(), 1.0, 1.0);
            const CGPoint origin = CGPointMake(-left, -bottom);
            CTFontDrawGlyphs(font.get(), &glyph, &origin, 1, context.get());
            CGContextFlush(context.get());
            if (!whitespace && std::none_of(bitmap.pixels.begin(), bitmap.pixels.end(), [](unsigned char pixel) { return pixel != 0; })) {
                std::fprintf(stderr, "System font produced no pixels for %s U+%04X\n", name.c_str(), codepoint);
                complete = false;
                continue;
            }
            bitmap.rectangle = atlas->AddCustomRectFontGlyph(primary, static_cast<ImWchar>(codepoint),
                bitmap.width, bitmap.height, static_cast<float>(advance.width),
                ImVec2(static_cast<float>(left), std::round(primary->Ascent) - static_cast<float>(top)));
            pending.SetBit(codepoint);
            bitmaps.push_back(std::move(bitmap));
            if ((!reported && codepoint >= 0x3000) || codepoint == 0x6742) {
                source_name = name;
                CFOwned<CFTypeRef> attribute(CTFontCopyAttribute(font.get(), kCTFontURLAttribute));
                CFOwned<CFStringRef> path(attribute && CFGetTypeID(attribute.get()) == CFURLGetTypeID()
                    ? CFURLCopyFileSystemPath(static_cast<CFURLRef>(attribute.get()), kCFURLPOSIXPathStyle) : nullptr);
                std::fprintf(stderr, "COV system font %s U+%04X: %s [%s]\n", core_text_string(language).c_str(), codepoint, name.c_str(), core_text_string(path.get()).c_str());
                reported = true;
            }
        }
    }
    return complete;
}

bool build_system_glyph_atlas(ImFont* primary, const float pixel_size,
                              const ImWchar* chinese, const ImWchar* japanese,
                              bool& zh_loaded, bool& ja_loaded,
                              std::string& chinese_name, std::string& japanese_name) {
    // A newly added font receives its ContainerAtlas during its first build.
    auto* atlas = ImGui::GetIO().Fonts;
    if (!atlas->Build()) return false;
    ImFontGlyphRangesBuilder pending;
    std::vector<SystemGlyphBitmap> bitmaps;
    zh_loaded = add_system_glyphs(primary, pixel_size, chinese, CFSTR("zh-Hans"), CFSTR("PingFangSC-Regular"), pending, bitmaps, chinese_name);
    ja_loaded = add_system_glyphs(primary, pixel_size, japanese, CFSTR("ja"), CFSTR("HiraginoSans-W3"), pending, bitmaps, japanese_name);
    if (!atlas->Build()) return false;
    // Custom rectangles own atlas locations; local bitmaps are copied into the
    // atlas before its RGBA upload is generated. No native pointers survive.
    for (const auto& bitmap : bitmaps) {
        const auto* rectangle = atlas->GetCustomRectByIndex(bitmap.rectangle);
        if (!rectangle || !rectangle->IsPacked()) return false;
        for (int row = 0; row < bitmap.height; ++row) {
            std::copy_n(bitmap.pixels.data() + row * bitmap.width, bitmap.width,
                        atlas->TexPixelsAlpha8 + (rectangle->Y + row) * atlas->TexWidth + rectangle->X);
        }
    }
    return zh_loaded && ja_loaded;
}
#endif

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
        kIntegrationChineseGlyphs+" "+kOrbitalDiagramChineseGlyphs+" "+
        orbital_ui_glyph_seed(Language::ChineseSimplified);
    static const std::string japanese=std::string(kSupplementalJapanese)+" "+
        kIntegrationJapaneseGlyphs+" "+kOrbitalDiagramJapaneseGlyphs+" "+
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
    std::string chinese = "macOS Chinese";
    std::string japanese = "macOS Japanese";
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
    const auto nbo_en=nbo_glyph_seed(Language::English);
    const auto nbo_fr=nbo_glyph_seed(Language::French);
    latin_builder.AddText(nbo_en.c_str());
    latin_builder.AddText(nbo_fr.c_str());
    const auto aomo_en=nbo_aomo_glyph_seed(Language::English);
    const auto aomo_fr=nbo_aomo_glyph_seed(Language::French);
    latin_builder.AddText(aomo_en.c_str());
    latin_builder.AddText(aomo_fr.c_str());
    latin_builder.AddText(scientific_glyph_seed());
    const auto inspection_seed=inspection_glyph_seed();
    latin_builder.AddText(inspection_seed.c_str());
    latin_builder.AddText(kIntegrationLatinGlyphs);
    latin_builder.AddText(kOrbitalDiagramLatinGlyphs);
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
    const auto nbo_zh=nbo_glyph_seed(Language::ChineseSimplified);
    const auto nbo_ja=nbo_glyph_seed(Language::Japanese);
    zh_builder.AddText(nbo_zh.c_str());
    ja_builder.AddText(nbo_ja.c_str());
    const auto aomo_zh=nbo_aomo_glyph_seed(Language::ChineseSimplified);
    const auto aomo_ja=nbo_aomo_glyph_seed(Language::Japanese);
    zh_builder.AddText(aomo_zh.c_str());
    ja_builder.AddText(aomo_ja.c_str());
    zh_builder.AddText(scientific_glyph_seed());
    zh_builder.AddText(inspection_seed.c_str());
    ja_builder.AddText(inspection_seed.c_str());
    ja_builder.AddText(scientific_glyph_seed());
    ImVector<ImWchar> zh_ranges;
    ImVector<ImWchar> ja_ranges;
    zh_builder.BuildRanges(&zh_ranges);
    ja_builder.BuildRanges(&ja_ranges);

    bool zh_loaded = false;
    bool ja_loaded = false;
#ifdef __APPLE__
    const bool built = build_system_glyph_atlas(primary, pixel_size, zh_ranges.Data, ja_ranges.Data,
                                               zh_loaded, ja_loaded, chinese, japanese);
#else
    merge_font_if_available(chinese, pixel_size, zh_ranges.Data, zh_loaded);
    merge_font_if_available(japanese, pixel_size, ja_ranges.Data, ja_loaded);
    const bool built = io.Fonts->Build();
#endif
    g_font_status = file_name_or_default(base, "ImGui default");
    g_font_status += " + ";
    g_font_status += zh_loaded ? file_name_or_default(chinese, "CJK") : "ZH fallback missing";
    g_font_status += " + ";
    g_font_status += ja_loaded ? file_name_or_default(japanese, "CJK") : "JA fallback missing";
    return built;
}

const char* font_status() noexcept { return g_font_status.c_str(); }
const char* font_status(Language language) {
    static std::string display;
    display=g_font_status;
    const char* zh=language==Language::ChineseSimplified?"缺少中文字体":
        language==Language::Japanese?"中国語フォントなし":
        language==Language::French?"Police chinoise absente":"Chinese font unavailable";
    const char* ja=language==Language::ChineseSimplified?"缺少日文字体":
        language==Language::Japanese?"日本語フォントなし":
        language==Language::French?"Police japonaise absente":"Japanese font unavailable";
    for(const auto& pair:{std::pair{"ZH fallback missing",zh},std::pair{"JA fallback missing",ja}}) {
        if(const auto pos=display.find(pair.first);pos!=std::string::npos)
            display.replace(pos,std::char_traits<char>::length(pair.first),pair.second);
    }
    return display.c_str();
}

void section_title(const char* label) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.58f, 0.68f, 0.82f, 1.0f));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
}

void begin_card(const char* id, const float height) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(12.0f, 11.0f));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.082f, 0.108f, 0.145f, 0.94f));
    ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.18f, 0.23f, 0.30f, 0.90f));
    if(height<=0)
        ImGui::BeginChild(id, ImVec2(0,0),
            ImGuiChildFlags_Border|ImGuiChildFlags_AutoResizeY|ImGuiChildFlags_AlwaysAutoResize,
            ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoScrollWithMouse);
    else
        ImGui::BeginChild(id, ImVec2(0.0f, height), true, ImGuiWindowFlags_None);
    if(height<=0)cov::validation::field(std::string("layout.card.")+id+".scroll_max_y",
        std::to_string(ImGui::GetScrollMaxY()));
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
