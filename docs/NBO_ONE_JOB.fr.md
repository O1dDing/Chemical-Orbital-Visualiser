# Préparer les fichiers de calcul pour COV

[English](NBO_ONE_JOB.md) · [简体中文](NBO_ONE_JOB.zh-CN.md) · [日本語](NBO_ONE_JOB.ja.md) · **Français**

Les fichiers FCHK/FCH ou Molden fournissent les énergies, les occupations et les formes des orbitales moléculaires. Les formes des orbitales NBO et l’analyse de leur composition nécessitent le rapport NBO, l’archive `.47` et les matrices orbitalaires correspondants. Les fichiers nécessaires à chaque vue sont indiqués dans [Utiliser les résultats NBO](AOMO_NBO.fr.md).

## Du calcul à COV

La fonction d’onde est obtenue en premier. Gaussian peut ensuite lancer NBO dans le même travail. `formchk` exporte le fichier checkpoint canonique au format FCHK ; GenNBO peut lire l’archive `.47` enregistrée pour produire un autre rapport NBO et les matrices demandées sans refaire le calcul SCF.

Le FCHK et les données NBO décrivent la même géométrie, la même base et le même état électronique. Une nouvelle analyse GenNBO peut utiliser d’autres options d’analyse ; son rapport et ses matrices doivent donc être utilisés ensemble. Des calculs d’optimisation et de fréquences peuvent précéder cette étape si nécessaire ; ils ne sont pas requis pour simplement afficher une fonction d’onde existante.

Le [modèle de travail unique](../examples/nbo-one-job/cov_one_job.py) lance successivement Gaussian/NBO, `formchk` et GenNBO, puis rassemble leurs sorties pour COV. L’application COV ouvre les résultats.

## Utiliser le modèle Windows

Le modèle s’exécute depuis l’arborescence des sources de COV avec Python 3.11 ou ultérieur, Gaussian 16W A.03 et NBO7 i8. Gaussian et NBO doivent être déjà installés, et `gaunbo6.bat` doit pointer vers cette installation NBO. A.03 utilise `Pop=NBO6Read` pour appeler NBO7. Les autres révisions de Gaussian ont leur propre configuration d’interface ; ce modèle vérifie la présence d’A.03.

Le modèle utilise [tests/validation_process.py](../tests/validation_process.py) dans l’arborescence des sources.

Copiez [water.json](../examples/nbo-one-job/water.json) et définissez la molécule, la géométrie, la méthode, la base, la charge et la multiplicité. L’exemple est un calcul ponctuel sur l’eau à géométrie fixe avec PBE1PBE / 6-31G(d), de charge 0 et de multiplicité 1. `basis_ecp_tail` accepte le bloc explicite de base/ECP si nécessaire.

Exécutez la commande depuis la racine des sources, avec vos chemins d’installation et un nouveau dossier de sortie :

```powershell
python examples/nbo-one-job/cov_one_job.py `
  --recipe examples/nbo-one-job/water.json `
  --output work/water-job `
  --gaussian-bin 'C:/Gaussian/Gaussian 16 W' `
  --nbo-bin 'C:/NBO/bin'
```

Cette commande présente l’entrée et vérifie les chemins sans créer de fichiers ni lancer de calcul. Ajoutez `--run` à la même commande pour le lancer. Le dossier de sortie ne doit pas déjà exister. Le modèle conserve la géométrie fournie et utilise une seule étape de calcul ponctuel ; les calculs d’optimisation et de fréquences nécessitent leurs propres entrées.

Le modèle autorise jusqu’à trois threads Gaussian. L’arbre des processus est limité à trois cœurs et 20 GiB ; les délais maximaux sont de 900 secondes pour Gaussian, 300 pour `formchk` et 600 pour GenNBO. `--memory-limit-gib` et `--gaussian-timeout` permettent d’ajuster les limites documentées ; `--help` liste les options.

## Ouvrir le résultat

Ouvrez `work/water-job/cov-package` dans COV, ou glissez son fichier `drop.covnbopkg` dans COV. Le dossier contient `canonical.fchk`, un rapport `analysis.nbo`, l’archive `.47` et les matrices orbitalaires demandées.

Le rapport, l’archive et les matrices dans `cov-package` proviennent de la même analyse GenNBO. `canonical.fchk` et `analysis.nbo` sont les noms utilisés par le modèle ; COV n’exige pas ces noms exacts. Le fichier `.covnbopkg` est une liste de fichiers à ouvrir, pas une source supplémentaire de données orbitalaires.

Le modèle conserve le checkpoint original et l’archive Gaussian. Les enregistrements d’exécution restent hors du dossier d’entrée. Les travaux qui échouent conservent leurs sorties ; une nouvelle tentative utilise un nouveau dossier.

## Sortie des matrices

Voici les réglages de sortie du modèle. Les numéros peuvent être modifiés dans NBO ; les en-têtes des matrices identifient leur contenu. Le modèle demande un ensemble étendu de fichiers, même si une vue COV donnée peut n’en utiliser qu’une partie.

| Option de sortie | Fichier | Contenu |
| --- | --- | --- |
| `AONBO=W37` | `FILE.37` | Coefficients des NBO dans la base AO |
| `NBOMO=W49` | `FILE.49` | Coefficients des MO dans la base NBO |
| `SAO=W50` | `FILE.50` | Matrice de recouvrement AO |
| `NAOMO=W51` | `FILE.51` | Coefficients des MO dans la base NAO |
| `AONAO=W52` | `FILE.52` | Coefficients des NAO dans la base AO |
| `NAONBO=W53` | `FILE.53` | Coefficients des NBO dans la base NAO |
| `AONHO=W54` | `FILE.54` | Coefficients des NHO dans la base AO |
| `AONLMO=W55` | `FILE.55` | Coefficients des NLMO dans la base AO |
| `AOPNAO=W56` | `FILE.56` | Coefficients des NAO avant orthogonalisation dans la base AO |
| `NAONHO=W57` | `FILE.57` | Coefficients des NHO dans la base NAO |
| `NAONLMO=W58` | `FILE.58` | Coefficients des NLMO dans la base NAO |
| `NHONBO=W59` | `FILE.59` | Coefficients des NBO dans la base NHO |
| `NBONLMO=W60` | `FILE.60` | Coefficients des NLMO dans la base NBO |
| `NLMOMO=W61` | `FILE.61` | Coefficients des MO dans la base NLMO |
| `AOMO=W62` | `FILE.62` | Coefficients des MO canoniques dans la base AO |

`ARCHIVE` écrit le fichier `.47`. `PRINT=3` inclut les sorties NLMO et Wiberg ; le modèle demande aussi `E2PERT=0.0`. Sa liste complète de mots-clés NBO figure dans [cov_one_job.py](../examples/nbo-one-job/cov_one_job.py).

Un fichier `.47` ne contient pas nécessairement de données de Fock. Ce modèle s’arrête si l’archive ne contient pas les sections de recouvrement, de densité, de MO canoniques ou de Fock, car il prépare aussi les entrées pour les vues d’énergie. COV peut utiliser un ensemble de fichiers plus réduit pour les vues qui n’ont pas besoin des données manquantes. Les vues à couche ouverte nécessitent les données de spin correspondantes ; le modèle ne choisit pas de méthode restreinte ou non restreinte pour la molécule.

## Conserver et récupérer les fichiers

L’entrée originale, le journal du calcul et le checkpoint sont utiles pour de futurs exports. Ils ne sont pas tous nécessaires pour ouvrir un ensemble de fichiers déjà exportés. Le FCHK canonique provient d’un checkpoint qui conserve les MO canoniques ; SaveNBOs ou SaveNLMOs peuvent remplacer ces orbitales dans un checkpoint.

| Élément manquant | Comment le récupérer |
| --- | --- |
| FCHK, lorsque le CHK Gaussian original est encore disponible | Exportez avec `formchk`, ou laissez COV appeler le convertisseur installé. |
| Section du rapport ou matrice orbitale, lorsque les données nécessaires sont encore dans `.47` | Relancez GenNBO avec les options de sortie supplémentaires ; utilisez ensemble le rapport et les matrices de cette exécution. Un nouveau calcul SCF n’est généralement pas nécessaire. |
| Données de Fock ou autres données absentes de l’archive | Exportez depuis l’état du calcul original si le logiciel le permet. GenNBO ne peut pas fournir des données absentes de son entrée. |
| Liste de fichiers à ouvrir `.covnbopkg` | Ouvrez le dossier ou sélectionnez directement les fichiers. |
| Fichiers appartenant à différents calculs | Retrouvez la fonction d’onde et les sorties NBO correspondantes. Renommer les fichiers ne change pas leur contenu. |
| Seules les coordonnées restent disponibles | Un nouveau calcul de structure électronique est nécessaire pour obtenir la fonction d’onde. |

COV construit les SALC à partir des orbitales importées et de la symétrie moléculaire ; elles ne nécessitent pas d’autre calcul SCF. Leurs formes et leurs énergies ont des besoins distincts, décrits dans [Utiliser les résultats NBO](AOMO_NBO.fr.md).

La sortie des matrices et la nouvelle analyse à partir de `.47` sont décrites dans le [manuel NBO 7](https://nbo.chem.wisc.edu/nboman.pdf), sections B.2.4–B.2.6 et B.7.
