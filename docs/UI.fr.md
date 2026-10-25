# Utiliser Chemical Orbital Visualiser

[English](UI.md) · [简体中文](UI.zh-CN.md) · [日本語](UI.ja.md) · **Français**

Ce guide décrit COV v0.4.0. Les [instructions de téléchargement](../README.fr.md#télécharger) indiquent comment lancer les paquets distribués.

Chemical Orbital Visualiser (COV) réunit les énergies, les occupations et les relations orbitalaires dans une même vue. Ouvrez un fichier Gaussian FCHK/FCH ou Molden existant, sélectionnez une orbitale dans la liste, comparez les niveaux dans le diagramme et exportez une figure ou un tableau de données. La molécule et l’orbitale sélectionnée peuvent aussi être examinées en 3D.

## Ouvrir un calcul

Utilisez la commande d’ouverture, glissez un fichier dans la fenêtre ou saisissez son chemin. COV lit les résultats de calcul ; il ne lance pas Gaussian. Un fichier Gaussian CHK peut être converti localement si `formchk` est installé.

Si COV ne trouve pas `formchk`, définissez la variable d’environnement `COV_FORMCHK` avec le chemin complet de son exécutable, puis redémarrez COV. Vous pouvez aussi créer un FCHK avec `formchk`, puis ouvrir ce fichier.

Ouvrir le dossier du calcul permet de charger ensemble la fonction d’onde et les fichiers NBO. Si le dossier contient plusieurs calculs, choisissez celui à ouvrir. La limite est de 100 atomes par molécule.

## Parcourir les orbitales et les énergies

Recherchez dans la liste des orbitales, accédez à la HOMO ou à la LUMO, ou filtrez les orbitales de cœur, de valence et virtuelles. Le choix de toutes les orbitales rétablit la liste complète importée ; les filtres ne suppriment aucune orbitale du calcul. Les détails de l’orbitale sélectionnée indiquent son énergie, son occupation, son spin et les noms de symétrie ou les étiquettes chimiques disponibles.

Les unités d’énergie proposées sont Ha, eV, J/mol, kJ/mol, cal/mol et kcal/mol. L’énergie d’origine reste en Hartree ; changer l’unité d’affichage ne modifie pas le calcul.

Les noms tels que `2a₁` sont affichés lorsque la symétrie et le rang dans cette symétrie sont disponibles. Sinon, COV affiche le numéro MO d’origine et le spin. La recherche accepte les noms et les numéros d’origine. Les niveaux proches peuvent être regroupés à leur énergie moyenne ; chaque orbitale conserve sa propre énergie et son identité dans les détails. La tolérance de regroupement commence à `1e-5 Ha` et peut être modifiée.

## Comparer les niveaux dans le diagramme

Le diagramme central montre les niveaux de valence occupés et les niveaux inoccupés proches, avec les occupations électroniques. Les niveaux de cœur profonds et les niveaux virtuels éloignés peuvent être masqués pour faciliter la lecture ; ils restent dans la liste importée. Sélectionner une orbitale pour la voir en 3D ne l’ajoute pas au diagramme compact et ne réordonne pas ce diagramme.

L’axe d’énergie peut être linéaire ou non linéaire. La vue non linéaire écarte les niveaux serrés et indique que l’axe est non linéaire, avec les valeurs d’énergie et leurs unités sur les graduations. Les niveaux regroupés utilisent l’énergie moyenne de leurs membres.

**Exporter les images** enregistre les figures PNG et SVG. **Données d’analyse (avancé) → Exporter les données d’analyse** enregistre séparément les données JSON et CSV. Les figures utilisent les noms d’orbitales affichés et la langue actuelle de l’interface. Les données comprennent les numéros des orbitales, les énergies, les occupations, le spin, les regroupements et les niveaux présents sur la figure.

## Examiner la molécule et l’orbitale

Faites pivoter la scène et réglez la visibilité des atomes, des liaisons, de la molécule et de l’orbitale. La représentation par défaut **Boules et bâtonnets renforcés** met en évidence le squelette moléculaire. Le mode **Bâtons + délocalisation** réduit l’encombrement des atomes. La visibilité des hydrogènes se règle séparément. Les phases positive et négative de la surface orbitale ont des couleurs différentes.

Certains styles de liaison ou de contact dépendent des informations présentes dans le calcul. Une ligne ou une couleur facilite la lecture et ne constitue pas, à elle seule, une mesure d’ordre de liaison.

La fenêtre de détails peut être déplacée, redimensionnée et rouverte. L’interface est disponible en anglais, chinois simplifié, japonais et français.

Pour les vues AO–MO et NBO, consultez [Utiliser les résultats NBO](AOMO_NBO.fr.md).
