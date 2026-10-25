# Chemical Orbital Visualiser (COV)

[English](README.md) · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · **Français**

Explorez les énergies et les occupations des orbitales, ainsi que les liens entre orbitales atomiques et moléculaires.

COV présente les énergies et les occupations des orbitales dans des diagrammes de niveaux d’énergie interactifs. Vous pouvez exporter les figures et les données. Il lit les fichiers Gaussian FCHK/FCH et Molden.

COV v0.4 ajoute l’analyse NBO. Avec la sortie NBO du calcul, vous pouvez voir comment les orbitales atomiques et localisées contribuent aux orbitales moléculaires, suivre leurs liens dans le diagramme et examiner les charges et les informations sur les liaisons.

Sélectionnez une orbitale pour voir sa forme en 3D. COV fonctionne sous Windows, macOS et Linux, avec le calcul sur CPU disponible sur chaque plateforme. L’interface est disponible en anglais, chinois simplifié, japonais et français.

## Fonctionnalités

- **Lire les diagrammes de niveaux d’énergie.** Consultez les énergies et les occupations électroniques, changez d’unité d’énergie et exportez le diagramme en PNG ou SVG.
- **Lire la symétrie et la composition.** Consultez les noms de symétrie et de groupes de la vue actuelle, les énergies communes des orbitales correspondantes des calculs à couche ouverte restreints, et la composition métal/ligand des complexes concernés. La vue simplifiée peut filtrer les orbitales de cœur et les contributions de fond des ligands avec les contributions AO/SALC correspondantes.
- **Examiner les détails du diagramme.** Les groupes dégénérés et les liaisons équivalentes ont un affichage de liaison cohérent ; les détails flottants des orbitales et les commandes d’étiquettes facilitent la lecture.
- **Suivre les liens entre orbitales.** Voyez comment les orbitales atomiques et localisées contribuent aux orbitales moléculaires. Sélectionnez un lien pour examiner sa contribution.
- **Examiner les charges et les liaisons.** Consultez les charges NPA, les populations de spin, les indices de liaison de Wiberg et les interactions donneur–accepteur lorsque les fichiers NBO contiennent ces données.
- **Exporter les figures et les données.** Enregistrez les diagrammes de niveaux d’énergie en PNG ou SVG, et les données correspondantes en CSV ou JSON.
- **Parcourir les orbitales.** Accédez à la HOMO ou à la LUMO, recherchez dans la liste des orbitales et affichez les orbitales occupées, virtuelles, de cœur ou de valence.
- **Voir les orbitales en 3D.** Sélectionnez une orbitale, réglez son isosurface, puis faites pivoter la vue moléculaire ou zoomez.

## Télécharger

| Version | Contenu | Téléchargements |
|---|---|---|
| [Version stable v0.4.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0) | Analyse NBO, composition et liens des orbitales, vues de symétrie et de liaison | [Windows ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0/Chemical-Orbital-Visualiser-v0.4.0-Windows-x64.zip) · [macOS ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0/Chemical-Orbital-Visualiser-v0.4.0-macOS-universal.zip) · [Linux tar.gz](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0/Chemical-Orbital-Visualiser-v0.4.0-Linux-x86_64.tar.gz) |

## Premiers pas

La [bibliothèque de molécules](examples/molecules/README.fr.md) propose 47 variantes de calcul, dont l’eau, le méthane, le benzène et le ferrocène en conformation décalée. Choisissez un exemple dans le catalogue, téléchargez son archive, puis ouvrez le dossier de calcul extrait dans COV.

1. Téléchargez et extrayez le paquet de votre plateforme. Sous Windows, lancez `cov.exe` ; sous macOS, ouvrez `Chemical Orbital Visualiser.app` ; sous Linux, lancez `./cov`.
2. Ouvrez un fichier Gaussian FCHK/FCH ou un fichier Molden compatible, ou glissez-le dans la fenêtre.
3. Parcourez les énergies et les occupations dans la liste des orbitales et le diagramme. Sélectionnez un niveau pour voir son orbitale, ou exportez le diagramme.
4. Ouvrez le dossier du calcul pour charger ensemble les fichiers de fonction d’onde et NBO. Si le dossier contient plusieurs calculs, choisissez celui à ouvrir.

Si `formchk` est installé, COV peut convertir les fichiers Gaussian CHK en FCHK. Les fichiers FCHK/FCH et Molden fournissent les énergies, les occupations et les formes des MO. Les formes des orbitales NBO et l’analyse de leur composition nécessitent aussi le rapport, l’archive `.47` et les matrices orbitalaires correspondants.

## Fichiers d’entrée et configuration requise

- **Fonctions d’onde :** fichiers Gaussian `.fchk` / `.fch`, ou fichiers compatibles `.molden` / `.mol` / `.input`. Les fichiers Gaussian `.chk` peuvent être ouverts avec un utilitaire `formchk` installé.
- **Fichiers NBO :** [Préparer les fichiers de calcul](docs/NBO_ONE_JOB.fr.md) explique comment les produire ; [Utiliser les résultats NBO](docs/AOMO_NBO.fr.md) indique les fichiers nécessaires à chaque vue.
- **Taille des molécules :** jusqu’à 100 atomes par fichier d’entrée.
- **Plateformes :** Windows x64 ; macOS 12 ou ultérieur sur Apple Silicon ou Intel ; Linux x86_64, avec Ubuntu 22.04 comme base du paquet.
- **Affichage et calcul :** OpenGL 2.1 ou ultérieur. Windows inclut CUDA 12.8 avec repli sur le CPU. CUDA nécessite un pilote NVIDIA compatible. macOS propose le calcul sur Metal et CPU ; Linux utilise le CPU. Un GPU NVIDIA n’est pas nécessaire pour ouvrir COV. Les modules GPU facultatifs peuvent être compilés depuis les sources.

## Documentation

- [Bibliothèque de molécules](examples/molecules/README.fr.md) — catalogue et téléchargements
- [Utiliser COV](docs/UI.fr.md)
- [Utiliser les résultats NBO](docs/AOMO_NBO.fr.md)
- [Préparer les fichiers de calcul](docs/NBO_ONE_JOB.fr.md) — inclut le modèle de travail dans l’arborescence des sources
- [Compiler depuis les sources](docs/BUILD.fr.md)
- [Notes de la version stable](docs/releases/v0.4.0.fr.md) · [Archives des préversions v0.4](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre-archive) · [Ancienne version stable v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0)

[Signaler un problème ou proposer une fonction](https://github.com/O1dDing/Chemical-Orbital-Visualiser/issues/new/choose).

## Licence

[Licence Apache 2.0](LICENSE). Les licences des bibliothèques fournies avec l’application figurent dans les [mentions relatives aux composants tiers](THIRD_PARTY_NOTICES.md).
