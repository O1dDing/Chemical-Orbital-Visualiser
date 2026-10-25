# Compiler Chemical Orbital Visualiser

[English](BUILD.md) · [简体中文](BUILD.zh-CN.md) · [日本語](BUILD.ja.md) · **Français**

Ces instructions concernent cette arborescence des sources. Les paquets distribués et leur configuration requise figurent sur la [page de téléchargement](../README.fr.md#télécharger).

## Prérequis

- CMake 3.28 ou ultérieur
- Compilateur C++20, Ninja et Git
- Accès à Internet pour télécharger les versions d’Eigen, GLFW et Dear ImGui fixées par CMake
- Pour l’application graphique : fichiers de développement OpenGL et, sous Linux, dépendances de développement X11 ou Wayland de GLFW. L’exécution nécessite un contexte de compatibilité OpenGL 2.1.

Sous Windows, utilisez une invite de commandes ou un PowerShell de développement Visual Studio 2022 x64 avec la charge de travail de développement Desktop en C++. Sous Debian/Ubuntu, utilisez une chaîne de compilation C++20 et les paquets de développement OpenGL/GLFW correspondants. Sous macOS, utilisez Xcode 16.2 ou ultérieur, avec son compilateur Metal. La compilation depuis les sources sous FreeBSD utilise le moteur CPU et nécessite les paquets de développement `mesa-libs` et `libglvnd`.

## Préréglages

Exécutez les commandes depuis la racine des sources :

```text
cmake --preset cpu-core
cmake --build --preset cpu-core --parallel
ctest --preset cpu-core
```

`cpu-core` compile les bibliothèques principales et les tests, sans interface graphique ni SDK GPU. Pour compiler l’application graphique avec le calcul des grilles sur CPU :

```text
cmake --preset portable
cmake --build --preset portable --parallel
```

Lancez `build/portable/cov` sous Linux/macOS ou `build\portable\cov.exe` sous Windows. Par exemple :

```bash
./build/portable/cov --compute-backend=cpu ./examples/h2.molden
```

Dans Windows PowerShell :

```powershell
.\build\portable\cov.exe --compute-backend=cpu .\examples\h2.molden
```

Les préréglages `cuda`, `metal` et `webgpu` activent le moteur correspondant. `metal` n’apparaît que sous macOS. Chaque préréglage graphique utilise un dossier de compilation distinct. Le préréglage CUDA nécessite CUDA Toolkit 12.8 ; son exécution nécessite un pilote NVIDIA compatible. Les cibles de code GPU par défaut sont `60, 61, 70, 75, 80, 86, 89, 90, 100, 120`, avec le PTX `compute_120`. Utilisez CUDA 12.8 pour cette liste ; CUDA 13 ne compile plus les cibles les plus anciennes. Exemple pour une seule architecture prise en charge :

```bash
cmake --preset cuda -DCMAKE_CUDA_ARCHITECTURES=86-real
cmake --build --preset cuda --parallel
```

Sous macOS, `cmake --preset metal` puis `cmake --build --preset metal` compilent l’application OpenGL 2.1 et son module de calcul Metal. Pour activer WebGPU, extrayez wgpu-native **v29.0.1.1**, puis exécutez `cmake --preset webgpu -DCOV_WGPU_NATIVE_ROOT=/absolute/path/to/wgpu-native` et `cmake --build --preset webgpu`. Sous Windows, la DLL d’exécution WebGPU doit aussi se trouver à côté de `cov.exe` au lancement.

## Modules de calcul AMD et Intel

HIP et SYCL utilisent des compilateurs distincts ; compilez chacun comme une bibliothèque partagée autonome. L’application principale reste compilée normalement en C++. Sous Linux, depuis la racine des sources :

```bash
# AMD HIP/ROCm : choisir une cible gfx prise en charge par le SDK installé.
cmake -S src/compute/hip -B build/hip -DCOV_ROOT="$PWD" \
  -DCMAKE_HIP_COMPILER=/opt/rocm/llvm/bin/clang++ -DCMAKE_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

# Intel oneAPI DPC++ : le module sélectionne les GPU Intel Level Zero.
source /opt/intel/oneapi/setvars.sh --force
cmake -S src/compute/sycl -B build/sycl -DCOV_ROOT="$PWD" -DCMAKE_CXX_COMPILER=icpx
cmake --build build/sycl --parallel
```

Sous Windows, utilisez une invite de développement Visual Studio 2022 x64 après avoir installé le SDK correspondant :

```bat
cmake -S src/compute/hip -B build/hip -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER="C:/Program Files/AMD/ROCm/7.2/bin/hipcc.exe" -DCOV_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
cmake -S src/compute/sycl -B build/sycl -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER=icx-cl
cmake --build build/sycl --parallel
```

Placez `cov_compute_hip.dll` ou `cov_compute_sycl.dll` sous Windows, ou `libcov_compute_hip.so` ou `libcov_compute_sycl.so` sous Linux, à côté de `cov.exe`/`cov`. Vous pouvez aussi définir `COV_COMPUTE_MODULE_DIR` avec le chemin **absolu** du dossier des modules. Sous Windows, placez les DLL d’exécution correspondantes à côté du module ou définissez `COV_COMPUTE_RUNTIME_DIR` avec le chemin absolu du dossier d’exécution du SDK. Sous Linux, les bibliothèques d’exécution HIP ou oneAPI doivent rester accessibles au chargeur système.

`--compute-backend=auto` privilégie le moteur natif correspondant au GPU d’affichage, essaie ensuite les autres moteurs natifs disponibles, puis utilise le CPU. Les choix explicites sont `cpu`, `cuda`, `hip`, `sycl`, `metal` et `webgpu`. `--compute-device=N` sélectionne le GPU d’un moteur par son indice commençant à zéro. WebGPU nécessite une sélection explicite. OpenCL est réservé à un développement ultérieur et n’a pas d’option de compilation dans cette arborescence.

## Fichiers d’entrée et lecture

COV lit les fonctions d’onde Gaussian FCHK/FCH et Molden. La conversion locale CHK nécessite un utilitaire `formchk` installé séparément ; définissez `COV_FORMCHK` avec son chemin si la détection automatique échoue. Un fichier Gaussian `.log` ou `.out` peut compléter les informations de la fonction d’onde correspondante, sans remplacer l’entrée FCHK/FCH.

La limite est de 100 atomes. Les fonctions de base `s/p/d/f/g` cartésiennes et sphériques réelles sont prises en charge. Pour Molden, le nombre de fonctions de base après développement doit correspondre aux coefficients de chaque orbitale moléculaire. Un indice de coefficient hors limites est traité comme une erreur d’entrée ou de convention des couches ; le lecteur ne tente pas de correction par supposition.

## Organisation du rendu

Le CPU lit l’entrée et prépare l’orbitale sélectionnée pour le calcul de grille. Le CPU ou le moteur de calcul choisi évalue une grille tridimensionnelle ; l’application la transfère ensuite dans une texture 3D OpenGL. CUDA peut aussi écrire directement dans la texture. Le rendu trace les isosurfaces positive et négative à partir de cette texture. Changer l’isovaleur modifie l’affichage sans recalculer la grille orbitale.

À chaque point de grille, l’évaluateur calcule :

$$
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r).
$$

Il n’alloue pas de matrice complète de taille « nombre de points de grille × nombre de fonctions de base ». Les résolutions proposées sont 64³, 128³, 256³ et 512³.

Pour préparer les fichiers Gaussian/NBO, consultez le [modèle de calcul en une tâche](NBO_ONE_JOB.fr.md).
