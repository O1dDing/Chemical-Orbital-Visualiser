#include "cov/orbital_ui_text.hpp"

#include <array>
#include <sstream>
#include <string>
#include <string_view>

namespace cov::ui {
namespace {

struct LocalisedString {
    const char* en;
    const char* zh;
    const char* ja;
    const char* fr;
};

constexpr auto kStrings = std::to_array<LocalisedString>({
    {"symmetry", "对称性", "対称性", "symétrie"},
    {"occupation", "占据数", "占有数", "occupation"},
    {"density", "密度", "密度", "densité"},
    {"overlap", "重叠", "重なり", "recouvrement"},
    {"bond order", "键级", "結合次数", "ordre de liaison"},
    {"point group", "点群", "点群", "groupe ponctuel"},
    {"Gaussian LOG/OUT enrichment attached", "已附加 Gaussian LOG/OUT 增补数据", "Gaussian LOG/OUT 補足データを適用済み", "Enrichissement Gaussian LOG/OUT associé"},
    {"producer", "来源数据", "生成元データ", "producteur"},
    {"derived", "推导数据", "導出データ", "dérivé"},
    {"unavailable", "不可用", "利用不可", "indisponible"},
    {"unknown", "未知", "不明", "inconnue"},

    {"Selected MO chemistry", "所选 MO 化学性质", "選択 MO の化学的性質", "Caractère chimique de l’OM sélectionnée"},
    {"Chemical-valence manifold", "化学价层轨道组", "化学原子価軌道空間", "Espace orbitalaire de valence chimique"},
    {"yes", "是", "はい", "oui"},
    {"no", "否", "いいえ", "non"},
    {"Valence AO composition", "价层 AO 组成", "原子価 AO 組成", "Composition AO de valence"},
    {"Atom-pair interactions", "原子对相互作用", "原子対相互作用", "Interactions par paire atomique"},
    {"Orbital family", "轨道类型", "軌道型", "Famille orbitale"},
    {"Bonding role", "成键性质", "結合性", "Caractère liant"},
    {"Multicentre family", "多中心族", "多中心族", "Famille multicentrique"},
    {"Delocalised pi family", "离域 π 族", "非局在化 π 族", "Famille π délocalisée"},
    {"Member MOs", "成员 MO", "構成 MO", "OM membres"},
    {"Participating atoms", "参与原子", "参加原子", "Atomes participants"},
    {"Participating electrons", "参与电子", "参加電子", "Électrons participants"},
    {"Donor / acceptor direction", "给体 / 受体方向", "供与体 / 受容体方向", "Direction donneur / accepteur"},
    {"Analysis method", "分析方法", "解析法", "Méthode d’analyse"},
    {"MO contribution", "MO 贡献", "MO 寄与", "Contribution OM"},
    {"MO contribution is overlap-population based; Mayer is the total density-level pair index.", "MO 贡献来自重叠布居；Mayer 是总密度层级的原子对指数。", "MO 寄与は重なり密度に基づき、Mayer は全密度レベルの原子対指数です。", "La contribution OM repose sur la population de recouvrement ; Mayer est l’indice de paire au niveau de la densité totale."},
    {"Outside the minimal valence reference", "最小价层参考之外", "最小原子価参照外", "Hors de la référence minimale de valence"},
    {"AO overlap unavailable", "AO 重叠不可用", "AO の重なりを利用できません", "Recouvrement AO indisponible"},
    {"delocalised-pi", "离域 π", "非局在化 π", "π délocalisée"},
    {"mixed", "混合", "混合", "mixte"},
    {"undetermined", "未确定", "未確定", "indéterminé"},
    {"bonding", "成键", "結合性", "liante"},
    {"antibonding", "反键", "反結合性", "antiliante"},
    {"nonbonding", "非键", "非結合性", "non liante"},

    {"Orientation channels", "取向通道", "配向チャネル", "Canaux d’orientation"},
    {"Topology", "拓扑", "トポロジー", "Topologie"},
    {"path-delocalised pi", "链状 π 离域", "鎖状 π 非局在化", "délocalisation π en chaîne"},
    {"cyclic-delocalised pi", "环状 π 离域", "環状 π 非局在化", "délocalisation π cyclique"},
    {"branched-resonance pi network", "支化共振 π 网络", "分岐共鳴 π ネットワーク", "réseau π à résonance ramifiée"},
    {"pi network with a skeletal spiro union", "含骨架螺连接的 π 网络", "骨格スピロ接合を持つ π ネットワーク", "réseau π avec jonction spiro du squelette"},
    {"haptic metal-pi", "多触点金属–π", "多点金属–π", "contact haptique métal–π"},
    {"symmetry-equivalent pi direct sum", "对称等价 π 子体系直和", "対称等価 π 部分系の直和", "somme directe de sous-systèmes π équivalents"},
    {"atoms", "原子", "原子", "atomes"},
    {"coherence", "相干度", "コヒーレンス", "cohérence"},
    {"Hide intermediate framework MOs", "隐藏中间框架轨道", "中間骨格軌道を非表示", "Masquer les OM intermédiaires"},
    {"Hide secondary ligand-centred, ligand-internal and polarisation MOs; keep the principal d levels, necessary sigma bonding/antibonding representatives, the selected MO and donor/acceptor pairs.", "隐藏次级配体中心、配体内部及极化型轨道；保留主要 d 能级、必要的 σ 成键/反键代表、当前轨道和给体/受体配对。", "副次的な配位子中心・配位子内部・分極軌道を隠し、主要 d 準位、必要な σ 結合/反結合代表、選択軌道と供与体/受容体対を保持します。", "Masque les OM secondaires centrées sur les ligands, internes aux ligands ou de polarisation ; conserve les niveaux d principaux, les représentants σ liants/antiliants nécessaires, l’OM sélectionnée et les paires donneur/accepteur."},

    {"local ligand field", "局部配体场", "局所配位子場", "champ de ligands local"},
    {"first shell", "第一配位层", "第一配位圏", "première sphère"},
    {"coordination geometry", "配位几何", "配位幾何", "géométrie de coordination"},
    {"coordination number", "配位数", "配位数", "nombre de coordination"},
    {"angular RMS (rad)", "角度 RMS（rad）", "角度 RMS（rad）", "RMS angulaire (rad)"},
    {"radial variation", "径向变异", "動径変動", "variation radiale"},
    {"local molecular geometries", "局部分子几何数", "局所分子幾何数", "géométries moléculaires locales"},
    {"atom", "原子", "原子", "atome"},
    {"local geometry", "局部几何", "局所幾何", "géométrie locale"},
    {"not centred on this MO", "该 MO 未中心化于这些原子", "この MO は該当中心上にありません", "non centrée sur cette OM"},
    {"group occupation", "轨道组占据数", "軌道群占有数", "occupation du groupe"},
    {"metal s / p / d", "金属 s / p / d", "金属 s / p / d", "métal s / p / d"},
    {"ligand p", "配体 p", "配位子 p", "ligand p"},
    {"sigma / pi channel", "σ / π 通道", "σ / π チャネル", "canal σ / π"},
    {"M–L overlap", "M–L 重叠", "M–L 重なり", "recouvrement M–L"},
    {"selection", "选取方式", "選択", "sélection"},
    {"recovered from complete raw MO block", "从完整原始 MO 数据块恢复", "完全な生 MO ブロックから復元", "récupérée depuis le bloc complet des OM brutes"},
    {"weak-field treatment", "弱场处理", "弱場処理", "traitement de champ faible"},
    {"approximately nonbonding", "近似非键", "近似的に非結合性", "approximativement non liante"},
    {"pi interaction", "π 相互作用", "π 相互作用", "interaction π"},
    {"MO group energy difference", "MO 组能量差", "MO 群のエネルギー差", "écart d’énergie des groupes MO"},
    {"interaction", "相互作用", "相互作用", "interaction"},
    {"splitting", "劈裂", "分裂", "séparation"},
    {"splitting (Ha)", "劈裂（Ha）", "分裂（Ha）", "séparation (Ha)"},

    {"direct", "直接数据", "直接データ", "directe"},
    {"parsed label", "标签解析", "ラベル解析", "étiquette interprétée"},
    {"derived", "推导", "導出", "dérivée"},
    {"heuristic", "启发式", "ヒューリスティック", "heuristique"},
    {"pi donation", "π 给电子", "π 供与", "don π"},
    {"pi back-donation", "π 回馈", "π 逆供与", "rétrodonation π"},
    {"pi coupling", "π 耦合", "π 相互作用", "couplage π"},
    {"weak-field split; approximately nonbonding", "弱场劈裂；近似非键", "弱場分裂；近似的に非結合性", "séparation de champ faible ; approximativement non liante"},

    {"MO diagram summary", "MO 图摘要", "MO 図概要", "Résumé du diagramme OM"},
    {"ligand-field valence groups", "配体场价层组", "配位子場原子価群", "groupes de valence du champ de ligands"},
    {"delocalised pi family groups", "离域 π 族组", "非局在化 π 族群", "groupes de la famille π délocalisée"},
    {"multicentre active-space groups", "多中心活性空间组", "多中心活性空間群", "groupes de l’espace actif multicentrique"},
    {"visible levels", "可见能级", "表示準位", "niveaux visibles"},
    {"occupied", "占据", "占有", "occupés"},
    {"virtual", "虚轨道", "仮想", "virtuels"},
    {"hidden", "隐藏", "非表示", "masqués"},
    {"local field", "局部场", "局所場", "champ local"},
    {"geometry", "几何", "幾何", "géométrie"},
    {"pi pairs", "π 配对", "π 対", "paires π"},
    {"protected overflow", "受保护行溢出", "保護行の超過", "dépassement protégé"},
    {"spin pairs", "自旋配对", "スピン対", "paires de spin"},
    {"unmatched visible", "可见未配对", "表示中の未対応", "visibles non appariées"},
    {"raw-MO recovered groups", "原始 MO 恢复组", "生 MO 復元群", "groupes récupérés des OM brutes"},
    {"compact active space unavailable", "紧凑活性空间不可用", "コンパクト活性空間を利用できません", "espace actif compact indisponible"},
    {"multiple pi orientation channels", "多方向 π 通道", "複数方向の π チャネル", "canaux π à orientations multiples"},
    {"Skeletal rings", "骨架环", "骨格環", "Cycles du squelette"},
    {"Connectivity model", "连接模型", "結合グラフのモデル", "Modèle de connectivité"},
    {"Mayer order and covalent-radius distance", "Mayer 键级与共价半径距离", "Mayer 結合次数と共有結合半径距離", "indice de Mayer et distance des rayons covalents"},
    {"Covalent-radius distance", "共价半径距离", "共有結合半径距離", "distance des rayons covalents"},
    {"Common ring atom", "两环共用原子", "両環の共有原子", "Atome commun aux deux cycles"},
    {"Ring path", "环路径", "環経路", "Chemin cyclique"},
    {"Additional connection between these rings", "两环另有连接", "両環間の追加接続", "Connexion supplémentaire entre ces cycles"},
    {"No spiro junction assigned to these π channels.", "这些 π 通道未分配螺连接。", "これらの π チャネルにスピロ接合の割り当てはありません。", "Aucune jonction spiro attribuée à ces canaux π."},
    {"Skeletal rings and cyclic π channels are recorded separately.", "骨架环与环状 π 通道分别记录。", "骨格環と環状 π チャネルは別々に記録。", "Cycles du squelette et canaux π cycliques enregistrés séparément."},
    {"Full-orbital label", "完整轨道标签", "全軌道のラベル", "Étiquette de l’OM entière"},
    {"Symmetry explanation", "对称性解释", "対称性の解釈", "Interprétation de symétrie"},
    {"local", "局部", "局所", "local"},
    {"candidate", "候选", "候補", "candidat"},
    {"source", "来源", "出力", "source"},
    {"whole MO", "整体", "全体", "OM entière"},
    {"mixed", "混合", "混合", "mixte"},
    {"Local angular projection", "局部角向投影", "局所角成分の射影", "Projection angulaire locale"},
    {"Candidate from shell dimension", "壳层维数候选", "殻の次元からの候補", "Candidat issu de la dimension de la couche"},
    {"Candidate from a pi partner", "π 伙伴候选", "π パートナーからの候補", "Candidat issu d’un partenaire π"},
    {"Candidate from the other spin", "另一自旋的候选", "他スピンからの候補", "Candidat issu de l’autre spin"},
    {"Molecular symmetry operations", "整体分子对称操作", "分子全体の対称操作", "Opérations de symétrie moléculaire"},
    {"MO subspace members", "MO 子空间成员", "MO 部分空間の要素", "Membres du sous-espace des OM"},
    {"Candidate source members", "候选来源成员", "候補の出典軌道", "Membres à l’origine du candidat"},
    {"Centre coverage of target subspace", "中心对目标子空间的覆盖率", "対象部分空間の中心被覆率", "Couverture du sous-espace par le centre"},
    {"Purity within the local shell", "局部壳层内纯度", "局所殻内の純度", "Pureté dans la couche locale"},
    {"Labelled local component / full target", "所标局部分量 / 完整目标", "ラベル付き局所成分 / 全対象", "Composante locale étiquetée / cible entière"},
    {"Subspace projection fractions, without occupation weighting.", "子空间投影比例，未按占据加权。", "占有数で重み付けしない部分空間射影比率。", "Fractions de projection de sous-espace, sans pondération par l’occupation."},
    {"Source label point group unresolved.", "来源标签对应的点群未确定。", "出力ラベルに対応する点群は未確定。", "Groupe ponctuel de l’étiquette source indéterminé."},
    {"Local reference axes in input coordinates", "输入坐标中的局部参考轴", "入力座標での局所基準軸", "Axes de référence locaux dans les coordonnées d’entrée"},

    {"Angular subspace / full target", "角向子空间 / 完整目标", "角部分空間 / 全対象", "Sous-espace angulaire / cible entière"},
    {"Represented spin-orbital rank", "已表示的自旋轨道秩", "表現されたスピン軌道のランク", "Rang représenté des spin-orbitales"},
    {"Maximum numerical residual", "最大数值残差", "最大数値残差", "Résidu numérique maximal"},
    {"No single local irrep assigned.", "未分配单一局部不可约表示。", "単一の局所既約表現の割り当てなし。", "Aucune représentation locale unique attribuée."},
    {"Crystal-field gap", "晶场能隙", "結晶場ギャップ", "Écart du champ cristallin"},
    {"crystal-field gaps", "晶场能隙数", "結晶場ギャップ数", "écarts du champ cristallin"},
    {"Expected ligand type", "配体类型预期", "予想される配位子の種類", "Type de ligand attendu"},
    {"sigma only", "仅 σ", "σ のみ", "σ uniquement"},
    {"pi donor", "π 给体", "π 供与体", "donneur π"},
    {"pi acceptor", "π 受体", "π 受容体", "accepteur π"},
    {"ambiguous", "不确定", "曖昧", "ambigu"},
    {"Orbital analysis / expectation", "轨道结果与预期", "軌道解析と予想", "Analyse orbitale et résultat attendu"},
    {"Matches expectation", "与预期一致", "予想と一致", "Conforme au résultat attendu"},
    {"Differs from expectation", "与预期不一致", "予想と異なる", "Différent du résultat attendu"},
    {"Local metal–ligand contribution; both level groups retained.", "局部金属–配体贡献；保留两端能级组。", "局所金属–配位子成分。両側の準位群を保持。", "Contribution locale métal–ligand ; les deux groupes de niveaux sont conservés."},
    {"Orbital details", "轨道详情", "軌道の詳細", "Détails de l’orbitale"},
    {"Close details", "关闭详情", "詳細を閉じる", "Fermer les détails"},
    {"Orbital details", "轨道详情", "軌道の詳細", "Détails de l’orbitale"},
    {"Selected orbital outside the current diagram.", "所选轨道不在当前能级图中。", "選択軌道は現在の準位図の範囲外です。", "Orbitale sélectionnée hors du diagramme actuel."},
    {"Data scope: the selected orbital and its level group in the current view.", "数据范围：当前视图中的所选轨道及其能级组。", "データの範囲：現在の表示で選択した軌道とその準位群。", "Portée des données : l’orbitale sélectionnée et son groupe de niveaux dans la vue actuelle."},
    {"Level group contains %zu MOs", "能级组含 %zu 个 MO", "準位群には %zu 個の MO が含まれます", "Le groupe de niveaux contient %zu OM"},
    {"Level-group representative: %s", "能级组代表：%s", "準位群の代表：%s", "Représentant du groupe de niveaux : %s"},
    {"Local metal–ligand group analysis", "局部金属–配体组分析", "局所金属–配位子群の解析", "Analyse locale métal–ligand du groupe"},
    {"This diagram’s local metal–ligand model is not applicable.", "本图的局部金属–配体模型不适用。", "この図の局所金属–配位子モデルは適用対象外です。", "Le modèle local métal–ligand de ce diagramme ne s’applique pas."},
    {"Local metal–ligand group data are unavailable.", "局部金属–配体组数据不可用。", "局所金属–配位子群のデータを利用できません。", "Les données locales métal–ligand du groupe sont indisponibles."},
    {"M–L sigma / pi (resolved channels)", "M–L σ / π（已解析通道）", "M–L σ / π（分解済みチャネル）", "M–L σ / π (canaux résolus)"},
    {"Related localized orbitals", "相关局域轨道", "関連する局在軌道", "Orbitales localisées associées"},
    {"Projected relations: %zu", "投影关系：%zu", "射影関係：%zu", "Relations projetées : %zu"},
    {"Canonical MO projection", "正则 MO 投影", "正準 MO への射影", "Projection sur l’OM canonique"},
    {"Localized E(2)", "局域 E(2)", "局在 E(2)", "E(2) localisée"},
    {"All related orbitals", "全部相关轨道", "関連する全軌道", "Toutes les orbitales associées"},
    {"π coupling", "π 耦合", "π カップリング", "Couplage π"},
    {"Coupled dimensions: %zu", "耦合维数：%zu", "結合次元：%zu", "Dimensions couplées : %zu"},
    {"Cross-Fock range", "交叉 Fock 范围", "交差 Fock 範囲", "Plage Fock croisée"},
    {"Metal / ligand projection", "金属 / 配体投影", "金属 / 配位子への射影", "Projection métal / ligand"},
    {"No π coupling assigned to this orbital", "该轨道未分配 π 耦合", "この軌道に π カップリングの割り当てはありません", "Aucun couplage π attribué à cette orbitale"},
    {"Orbital composition", "轨道组成", "軌道組成", "Composition orbitale"},
    {"NAO projection", "NAO 投影", "NAO 射影", "Projection NAO"},
    {"AO-reference projection", "AO 参考投影", "AO 参照射影", "Projection sur la référence AO"},
    {"Retained NAO norm", "保留 NAO 范数", "保持された NAO ノルム", "Norme NAO retenue"},
    {"AO reconstruction residual", "AO 重建残差", "AO 再構成残差", "Résidu de reconstruction AO"},
    {"Not analysed", "未分析", "未解析", "Non analysé"},
    {"Missing data", "数据不足", "データ不足", "Données manquantes"},
    {"Incompatible data", "数据不匹配", "データ不整合", "Données incompatibles"},
    {"Unsupported data", "暂不支持此数据", "未対応のデータ", "Données non prises en charge"},
    {"Not applicable", "不适用", "適用対象外", "Sans objet"},
    {"Not reported above the output threshold", "输出阈值以上未报告", "出力しきい値以上の報告なし", "Non rapporté au-dessus du seuil de sortie"},
    {"Undetermined", "未确定", "未確定", "Indéterminé"},
    {"Bonding mixing", "成键混合", "結合性混合", "Mélange liant"},
    {"Antibonding mixing", "反键混合", "反結合性混合", "Mélange antiliant"},
    {"Mixed coupling", "混合耦合", "混合カップリング", "Couplage mixte"},
    {"Symmetric coupling", "对称耦合", "対称カップリング", "Couplage symétrique"},
    {"Ligand → metal", "配体 → 金属", "配位子 → 金属", "Ligand → métal"},
    {"Metal → ligand", "金属 → 配体", "金属 → 配位子", "Métal → ligand"},
    {"Gaussian full point group", "Gaussian 完整点群", "Gaussian の全点群", "Groupe ponctuel complet Gaussian"},
    {"Gaussian Abelian subgroup", "Gaussian 阿贝尔子群", "Gaussian の可換部分群", "Sous-groupe abélien Gaussian"},
    {"Available", "可用", "利用可能", "Disponible"},
    {"Calculation failed", "计算失败", "計算失敗", "Échec du calcul"},
    {"Total density", "总密度", "全密度", "Densité totale"},
    {"Spin density", "自旋密度", "スピン密度", "Densité de spin"},
    {"Shared integer occupations: paired double occupations, alpha single occupations.", "共享整数占据：双占据成对，单占据为 α。", "共通の整数占有：二重占有は対、単占有は α。", "Occupations entières partagées : occupations doubles appariées, occupations simples alpha."},
    {"Spin partition and full Mayer bond orders unavailable.", "自旋分配与完整 Mayer 键级不可用。", "スピン分配と完全な Mayer 結合次数を利用できません。", "Répartition du spin et indices de liaison Mayer complets indisponibles."},
    {"Bonding analysis and derived orbital labels unavailable.", "成键分析与推导轨道标签不可用。", "結合解析と導出軌道ラベルを利用できません。", "Analyse des liaisons et étiquettes orbitalaires dérivées indisponibles."},
    {"Input overlap differs from the basis; using basis integrals.", "输入重叠与基组不一致；使用基组积分。", "入力の重なりは基底と不一致のため、基底積分を使用。", "Recouvrement d’entrée incompatible avec la base ; utilisation des intégrales de base."},
    {"Label source", "标签来源", "ラベルの由来", "Origine de l’étiquette"},
    {"Calculation file or folder", "计算文件或文件夹", "計算ファイルまたはフォルダー", "Fichier ou dossier de calcul"},
    {"Supports FCHK, FCH, CHK, Molden and NBO data.", "支持 FCHK、FCH、CHK、Molden 和 NBO 数据。", "FCHK、FCH、CHK、Molden、NBO データに対応。", "Formats FCHK, FCH, CHK, Molden et NBO pris en charge."},
    {"Source MO", "原文件 MO 编号", "入力 MO 番号", "Numéro OM source"},
});

static_assert(kStrings.size() == static_cast<std::size_t>(OrbitalText::Count));

constexpr auto kGeometryNames = std::to_array<std::pair<std::string_view, LocalisedString>>({
    {"L-2", {"Linear", "线形", "直線形", "linéaire"}},
    {"A-2", {"Angular", "折线形", "折れ線形", "coudée"}},
    {"TP-3", {"Trigonal planar", "平面三角形", "三角平面形", "trigonale plane"}},
    {"TPY-3", {"Trigonal pyramidal", "三角锥形", "三角錐形", "pyramidale trigonale"}},
    {"TS-3", {"T-shaped", "T 形", "T 字形", "en T"}},
    {"T-4", {"Tetrahedral", "四面体形", "四面体形", "tétraédrique"}},
    {"SP-4", {"Square planar", "平面正方形", "正方形平面形", "plan carré"}},
    {"SS-4", {"Seesaw", "跷跷板形", "シーソー形", "en bascule"}},
    {"vTBPY-4", {"Trigonal-pyramidal (axially vacant trigonal bipyramid)", "三角锥形（轴向空位三角双锥）", "三角錐形（軸位空孔三方両錐）", "pyramidale trigonale (bipyramide trigonale à lacune axiale)"}},
    {"TBPY-5", {"Trigonal bipyramidal", "三角双锥形", "三方両錐形", "bipyramidale trigonale"}},
    {"SPY-5", {"Square pyramidal", "四方锥形", "四角錐形", "pyramidale carrée"}},
    {"OC-6", {"Octahedral", "八面体形", "八面体形", "octaédrique"}},
    {"TPR-6", {"Trigonal prismatic", "三棱柱形", "三角柱形", "prismatique trigonale"}},
    {"PBPY-7", {"Pentagonal bipyramidal", "五角双锥形", "五方両錐形", "bipyramidale pentagonale"}},
    {"COC-7", {"Capped octahedral", "加帽八面体形", "一冠八面体形", "octaèdre coiffé"}},
    {"CTPR-7", {"Capped trigonal prismatic", "加帽三棱柱形", "一冠三角柱形", "prisme trigonal coiffé"}},
    {"SAPR-8", {"Square antiprismatic", "四方反棱柱形", "四角反柱形", "antiprismatique carrée"}},
    {"TDD-8", {"Triangular dodecahedral", "三角十二面体形", "三角十二面体形", "dodécaèdre triangulaire"}},
    {"BTPR-8", {"Bicapped trigonal prismatic", "双帽三棱柱形", "二冠三角柱形", "prisme trigonal bicoiffé"}},
    {"CSAPR-9", {"Capped square antiprismatic", "加帽四方反棱柱形", "一冠四角反柱形", "antiprisme carré coiffé"}},
    {"TCTPR-9", {"Tricapped trigonal prismatic", "三帽三棱柱形", "三冠三角柱形", "prisme trigonal tricoiffé"}},
    {"PPR-10", {"Pentagonal prism", "五棱柱形", "五角柱形", "prisme pentagonal"}},
    {"PAPR-10", {"Pentagonal antiprism", "五方反棱柱形", "五角反柱形", "antiprisme pentagonal"}},
    {"BCSAPR-10", {"Bicapped square antiprismatic", "双帽四方反棱柱形", "二冠四角反柱形", "antiprisme carré bicoiffé"}},
    {"SPC-10", {"Sphenocorona", "斯芬诺冠形", "スフェノコロナ形", "sphénocouronne"}},
    {"TD-10", {"Tetradecahedral (2:6:2)", "十四面体形（2:6:2）", "十四面体形（2:6:2）", "tétradécaédrique (2:6:2)"}},
});

struct MachineTranslation {
    std::string_view machine;
    LocalisedString text;
};

constexpr auto kChemistryMethods = std::to_array<MachineTranslation>({
    {"COV FCHK S-metric minimal atomic-reference projection",
     {"Minimal atomic-reference projection (S metric)",
      "最小原子参考投影（S 度量）",
      "最小原子参照射影（S 計量）",
      "Projection sur la référence atomique minimale (métrique S)"}},
    {"legacy AO-reference projection", {"AO-reference projection", "AO 参考投影", "AO 参照射影", "Projection sur la référence AO"}},
    {"validated NAO projection", {"NAO projection", "NAO 投影", "NAO 射影", "Projection NAO"}},
});

constexpr auto kChemistryNotes = std::to_array<MachineTranslation>({
    {"AO overlap metric unavailable",
     {"AO overlap metric unavailable", "AO 重叠度量不可用", "AO 重なり計量を利用できません", "Métrique de recouvrement des AO indisponible"}},
    {"COV minimal atomic reference could not be built",
     {"COV minimal atomic reference could not be built", "无法构建 COV 最小原子参考", "COV の最小原子参照を構築できませんでした", "La référence atomique minimale COV n’a pas pu être construite"}},
    {"AO-to-atom map or chemical-valence rank unavailable",
     {"AO-to-atom map or chemical-valence rank unavailable", "AO–原子映射或化学价层秩不可用", "AO–原子対応または化学原子価ランクを利用できません", "Correspondance AO–atome ou rang de valence chimique indisponible"}},
    {"No stable atom-pair interaction frame; chemistry remains UND",
     {"Atom-pair interaction unresolved", "原子对相互作用未确定", "原子対相互作用は未確定", "Interaction par paire atomique indéterminée"}},
    {"Outside the selected minimal chemical-valence canonical manifold",
     {"Outside the selected minimal chemical-valence canonical manifold", "位于所选最小化学价层正则轨道空间之外", "選択された最小化学原子価正準軌道空間の外側です", "Hors de l’espace canonique minimal de valence chimique sélectionné"}},
});

const char* localised(const LocalisedString& value, const Language language) noexcept {
    switch (language) {
        case Language::ChineseSimplified: return value.zh;
        case Language::Japanese: return value.ja;
        case Language::French: return value.fr;
        default: return value.en;
    }
}

template <std::size_t N>
std::string translate_machine_value(
    const std::string_view machine,
    const std::array<MachineTranslation, N>& values,
    const Language language) {
    for (const auto& value : values) {
        if (value.machine == machine) return localised(value.text, language);
    }
    return std::string(machine);
}

const char* mode_text(const MODiagramMode mode, const Language language) noexcept {
    switch (mode) {
        case MODiagramMode::DelocalisedPiFamilyOnly:
            return orbital_tr(OrbitalText::ModeDelocalisedPi, language);
        case MODiagramMode::MulticentreActiveSpaceOnly:
            return orbital_tr(OrbitalText::ModeMulticentre, language);
        default:
            return orbital_tr(OrbitalText::ModeValenceCentral, language);
    }
}

} // namespace

const char* orbital_tr(const OrbitalText key, const Language language) noexcept {
    const auto index = static_cast<std::size_t>(key);
    if (index >= kStrings.size()) return "";
    return localised(kStrings[index], language);
}

const char* localised_wavefunction_source(
    const WavefunctionSource source, const Language language) noexcept {
    switch (source) {
        case WavefunctionSource::Fchk: return "FCHK";
        case WavefunctionSource::Molden: return "Molden";
        default: return orbital_tr(OrbitalText::UnknownSource, language);
    }
}

const char* localised_data_provenance(
    const DataProvenance provenance, const Language language) noexcept {
    switch (provenance) {
        case DataProvenance::Producer:
            return orbital_tr(OrbitalText::Producer, language);
        case DataProvenance::Derived:
            return orbital_tr(OrbitalText::Derived, language);
        default:
            return orbital_tr(OrbitalText::Unavailable, language);
    }
}

const char* localised_annotation_source(
    const AnnotationSource source, const Language language) noexcept {
    switch (source) {
        case AnnotationSource::Direct:
            return orbital_tr(OrbitalText::AnnotationDirect, language);
        case AnnotationSource::ParsedLabel:
            return orbital_tr(OrbitalText::AnnotationParsedLabel, language);
        case AnnotationSource::Derived:
            return orbital_tr(OrbitalText::AnnotationDerived, language);
        case AnnotationSource::Heuristic:
            return orbital_tr(OrbitalText::AnnotationHeuristic, language);
        default:
            return orbital_tr(OrbitalText::Unavailable, language);
    }
}

const char* localised_bonding_class(
    const BondingClass value, const Language language) noexcept {
    switch (value) {
        case BondingClass::Bonding:
            return orbital_tr(OrbitalText::Bonding, language);
        case BondingClass::Nonbonding:
            return orbital_tr(OrbitalText::Nonbonding, language);
        case BondingClass::Antibonding:
            return orbital_tr(OrbitalText::Antibonding, language);
        case BondingClass::Mixed:
            switch(language) {
                case Language::ChineseSimplified:return "混合";
                case Language::Japanese:return "混合";
                case Language::French:return "Mixte";
                default:return "Mixed";
            }
        default:
            return "N/A";
    }
}

const char* localised_orbital_bonding_role(
    const OrbitalBondingRole value, const Language language) noexcept {
    switch (value) {
        case OrbitalBondingRole::Bonding:
            return orbital_tr(OrbitalText::Bonding, language);
        case OrbitalBondingRole::Antibonding:
            return orbital_tr(OrbitalText::Antibonding, language);
        case OrbitalBondingRole::Nonbonding:
            return orbital_tr(OrbitalText::Nonbonding, language);
        case OrbitalBondingRole::NotApplicable:
            return "N/A";
        default:
            return orbital_tr(OrbitalText::Mixed, language);
    }
}

const char* localised_pi_interaction_kind(
    const PiInteractionKind kind, const Language language) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor:
            return orbital_tr(OrbitalText::PiDonorSplitting, language);
        case PiInteractionKind::Acceptor:
            return orbital_tr(OrbitalText::PiAcceptorSplitting, language);
        case PiInteractionKind::WeakNearNonbonding:
            return orbital_tr(OrbitalText::PiWeakNearNonbonding, language);
        default:
            return orbital_tr(OrbitalText::PiCoupledSplitting, language);
    }
}

std::string localised_geometry_name(
    const std::string_view geometry_id,
    const std::string_view fallback_name,
    const Language language) {
    for (const auto& [id, name] : kGeometryNames) {
        if (id == geometry_id) return localised(name, language);
    }
    return std::string(fallback_name);
}

std::string localised_chemistry_method(
    const std::string_view method, const Language language) {
    return translate_machine_value(method, kChemistryMethods, language);
}

std::string localised_chemistry_note(
    const std::string_view note, const Language language) {
    return translate_machine_value(note, kChemistryNotes, language);
}

std::string localised_diagram_selection_summary(
    const MODiagramData& data, const Language language) {
    std::ostringstream out;
    if (data.selection.summary.starts_with("compact active space unavailable")) {
        out << orbital_tr(OrbitalText::CompactActiveSpaceUnavailable, language)
            << " · ";
    }
    out << orbital_tr(OrbitalText::DiagramSummary, language)
        << " · " << mode_text(data.mode, language)
        << ": " << orbital_tr(OrbitalText::VisibleLevels, language)
        << '=' << (data.selection.counts_are_final?data.selection.final_member_count:data.levels.size())
        << '/' << data.metadata.size()
        << "; " << orbital_tr(OrbitalText::Occupied, language)
        << '=' << data.selection.valence_occupied_count
        << "; " << orbital_tr(OrbitalText::Virtual, language)
        << '=' << data.selection.frontier_virtual_count
        << "; " << orbital_tr(OrbitalText::Hidden, language)
        << '=' << data.selection.hidden_count;

    if (!data.ligand_field_point_group.empty()) {
        out << "; " << orbital_tr(OrbitalText::LocalField, language)
            << '=' << data.ligand_field_point_group;
    }
    if (!data.ligand_field_geometry_id.empty()) {
        out << "; " << orbital_tr(OrbitalText::Geometry, language)
            << '=' << data.ligand_field_geometry_id
            << "; CN=" << data.ligand_field_coordination_number;
    }
    out << "; " << orbital_tr(OrbitalText::PiPairs, language)
        << '=' << data.pi_interactions.size()
        << "; " << orbital_tr(OrbitalText::CrystalFieldGaps, language)
        << '=' << data.crystal_field_gaps.size()
        << "; " << orbital_tr(OrbitalText::ProtectedOverflow, language)
        << '=' << data.selection.protected_overflow_count;
    if (data.spin_counterpart_pair_count > 0u) {
        out << "; " << orbital_tr(OrbitalText::SpinPairs, language)
            << '=' << data.spin_counterpart_pair_count
            << "; " << orbital_tr(OrbitalText::UnmatchedVisible, language)
            << '=' << data.spin_counterpart_unmatched_visible;
    }
    std::size_t raw_recovered = 0u;
    for (const auto& level : data.levels) {
        if (level.raw_data_fallback) ++raw_recovered;
    }
    if (raw_recovered > 0u) {
        out << "; " << orbital_tr(OrbitalText::RawRecoveredGroups, language)
            << '=' << raw_recovered;
    }
    return out.str();
}

const char* orbital_ui_glyph_seed(const Language language) noexcept {
    static const std::array<std::string, static_cast<std::size_t>(Language::Count)> seeds = [] {
        std::array<std::string, static_cast<std::size_t>(Language::Count)> result;
        for (std::size_t language_index = 0; language_index < result.size(); ++language_index) {
            const auto current = static_cast<Language>(language_index);
            auto& seed = result[language_index];
            const auto append = [&seed](const std::string_view text) {
                if (!seed.empty()) seed.push_back(' ');
                seed.append(text);
            };
            for (const auto& value : kStrings) append(localised(value, current));
            for (const auto& [id, value] : kGeometryNames) {
                append(id);
                append(localised(value, current));
            }
            for (const auto& value : kChemistryMethods) append(localised(value.text, current));
            for (const auto& value : kChemistryNotes) append(localised(value.text, current));
        }
        return result;
    }();
    const auto index = static_cast<std::size_t>(language);
    return index < seeds.size() ? seeds[index].c_str() : "";
}

} // namespace cov::ui
