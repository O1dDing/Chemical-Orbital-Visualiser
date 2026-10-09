# Explorer les liens AO–MO et NBO

[English](AOMO_NBO.md) · [简体中文](AOMO_NBO.zh-CN.md) · [日本語](AOMO_NBO.ja.md) · **Français**

Chemical Orbital Visualiser (COV) réunit les énergies, les occupations et la composition des orbitales dans une même vue. Dans la préversion v0.4.0, vous pouvez ouvrir un fichier Gaussian FCHK existant et ses résultats NBO correspondants, puis accéder aux orbitales et contributions liées à une MO, un atome ou une liaison. Examinez l’orbitale sélectionnée en 3D pour voir sa forme.

## Quels fichiers sont nécessaires ?

COV lit des résultats existants ; il ne lance pas Gaussian ni NBO. Glissez les fichiers ensemble ou leur dossier dans la fenêtre, ou saisissez un chemin. Si le dossier contient plusieurs calculs, choisissez celui à ouvrir. Les réglages avancés des chemins permettent d’ouvrir des fichiers stockés à différents endroits.

| Ce que vous souhaitez voir | Fichiers à fournir |
| --- | --- |
| Énergies et occupations des MO, diagramme de niveaux d’énergie, liens par les coefficients AO–MO de Gaussian et MO en 3D | Fichier Gaussian FCHK/FCH ; Molden permet les vues MO de base. |
| Charges NPA, indices de Wiberg et valeurs E(2) imprimés | FCHK et rapport NBO contenant ces sections. |
| Formes des NAO, coloration selon la charge NPA et les indices de Wiberg | Fichiers de base et AONAO, avec la densité correspondante. |
| Composition des MO en NAO | Fichiers de base, AONAO et NAOMO ; NAONBO fournit une autre transformation. |
| Formes des NBO et liens entre MO canoniques et NBO | Fichiers de base et AONBO ; les liens avec les MO utilisent les coefficients canoniques correspondants du FCHK. |
| Formes et composantes des NHO | Fichiers de base, AONAO, AONHO et NAONHO ; les liens NHO–NBO utilisent aussi NHONBO et AONBO. |
| Formes des NLMO, NBO principale et partie résiduelle | Fichiers de base, AONBO, AONLMO, NBONLMO et entrées correspondantes du rapport. Les liens avec les MO utilisent les coefficients canoniques correspondants du FCHK. |
| Formes des PNAO | Fichiers de base, AONAO et AOPNAO. |
| Formes des SALC | Fichiers de base et AONAO, avec une géométrie et un espace orbital sélectionné permettant la combinaison de symétrie. |
| Énergies des NAO et SALC calculées à partir de l’opérateur | Orbitales correspondantes et données de Fock, avec les composantes de spin nécessaires au calcul. |
| Orbitales donneur–accepteur des interactions E(2) | Fichiers de base, AONBO et entrées du rapport identifiant le donneur et l’accepteur. |

Les fichiers de base sont le FCHK, le rapport NBO et l’archive `.47` correspondants. Le FCHK et les données NBO doivent décrire la même fonction d’onde. Le rapport et les matrices orbitalaires doivent provenir de la même analyse NBO ; une nouvelle analyse GenNBO possède son propre rapport et ses propres matrices. L’archive fournit la structure, la base, le recouvrement, la densité et les coefficients des MO canoniques. La coloration selon le spin nécessite aussi la densité résolue en spin et les données du rapport correspondantes.

NBOMO et NLMOMO fournissent des transformations supplémentaires lorsqu’ils sont disponibles. COV peut aussi calculer les liens avec les MO à partir des coefficients des orbitales locales, du recouvrement et des coefficients des MO canoniques ; ces deux fichiers ne sont donc pas toujours nécessaires. Si les coefficients canoniques d’un spin manquent dans le FCHK, les liens avec ces MO sont indisponibles.

Les données de Fock ne sont pas nécessaires pour simplement lire les valeurs E(2) imprimées ou afficher les orbitales donneur et accepteur correspondantes. Lorsqu’elles sont disponibles, COV peut aussi comparer le couplage et la différence d’énergie imprimés avec la matrice. Les énergies des NAO/NBO imprimées dans le rapport peuvent rester disponibles sans matrice de Fock dans l’archive. Les formes des SALC et leurs énergies calculées ont des besoins distincts.

Les numéros des fichiers de matrices dépendent des réglages de sortie NBO ; COV identifie leur contenu grâce aux en-têtes des matrices. Notre modèle utilise les noms `canonical.fchk` et `analysis.nbo`, mais ces noms ne sont pas obligatoires. Un fichier `.covnbopkg` indique les fichiers à ouvrir ; vous pouvez aussi ouvrir le dossier ou sélectionner les fichiers.

Pour la séquence de calcul, les options de sortie et les fichiers manquants, consultez [Préparer les fichiers de calcul](NBO_ONE_JOB.fr.md). L’association NBO actuelle utilise les données Gaussian FCHK/FCH ; Molden reste utilisable pour les vues MO ordinaires. Ces vues restent disponibles si les données NBO manquent.

## Suivre une orbitale dans le diagramme

Choisissez une MO dans le navigateur ou le diagramme de niveaux d’énergie. Le diagramme se concentre sur les niveaux de valence et les niveaux inoccupés proches ; **Toutes** dans le navigateur affiche la liste complète des orbitales importées. Cliquez sur une MO pour mettre en évidence ses liens disponibles avec les orbitales atomiques Gaussian (AO), les orbitales atomiques naturelles (NAO) et les groupes d’atomes. Cliquez sur une AO, une NAO ou un lien pour examiner une contribution ou trouver les MO associées. Vous pouvez sélectionner plusieurs termes et afficher leur somme partielle signée ou les superposer en 3D. Replier les groupes simplifie le diagramme sans retirer leurs membres.

Les NAO forment une représentation orthogonale : les carrés de leurs coefficients peuvent donc décrire des poids dans cette représentation. Les carrés des coefficients des AO Gaussian d’origine, généralement non orthogonales, ne sont pas des populations atomiques. Lorsque les énergies des NAO ou SALC sont disponibles, elles peuvent utiliser le même axe numérique que les MO. Ce sont des valeurs moyennes de l’opérateur dans l’environnement moléculaire ; les valeurs centrales sont les énergies des MO canoniques. La disposition latérale illustrative organise les orbitales sans utiliser leurs énergies. Un groupe d’orbitales est une SALC uniquement si les données disponibles justifient cette interprétation de symétrie.

Une base locale peut contenir moins d’orbitales que la base AO Gaussian. Le diagramme affiche alors les contributions disponibles et la partie non couverte, sans ramener les contributions à 100 %.

**Vue d’ensemble** propose une vue compacte. **Analyse détaillée** ajoute des détails, et **Toutes les orbitales** inclut les orbitales de cœur et de Rydberg. Ces préréglages modifient les éléments affichés.

## Partir de la molécule

Cliquez sur un atome ou une liaison dans la vue 3D pour ouvrir les valeurs et les liens orbitaux associés. Lorsque le rapport et les matrices correspondantes sont disponibles, les atomes peuvent être colorés selon leur charge NPA ou leur population de spin ; les liaisons peuvent afficher les indices de Wiberg, les orbitales liantes et antiliantes, ainsi que les relations de coordination ou multicentriques. Un indice de Wiberg est une valeur continue, pas un ordre de liaison entier.

Les vues NHO peuvent montrer les lobes directionnels des orbitales et leurs composantes de moment cinétique. Sélectionnez une interaction E(2) pour voir ensemble ses orbitales donneur et accepteur ; la valeur E(2) est une estimation perturbative, pas une énergie de liaison ou de réaction. Les vues NLMO peuvent séparer l’orbitale complète, sa composante NBO principale et la partie résiduelle lorsque les données nécessaires sont disponibles.

Une faible contribution peut être masquée par le seuil d’isosurface actuel. Utilisez **Adapter** pour ajuster le seuil d’affichage ou **Voir les liaisons** pour rendre la molécule plus visible. Gardez le même seuil lorsque vous comparez la taille de différentes orbitales.

## Exporter

**Exporter les images** enregistre les figures PNG et SVG. **Données d’analyse (avancé) → Exporter les données d’analyse** enregistre séparément les fichiers JSON/CSV et les analyses associées. L’exploration interactive en 3D reste dans COV.
