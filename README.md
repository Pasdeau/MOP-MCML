# Project MOP-MCML (English Version)

First modifications by Songlin Li, improvements by Wenzheng Wang

---

## Launching with Visual Studio (VS) 2022

- Issue with version 2019 as the project .vcxproj needs to be opened, which is not possible with version 2022.
- Find the .sln file and open it with VS 2022 (open a project or solution).
- You may need to download the necessary add-ons (15 minutes).
- Depending on the versions, read the message in the bottom left window (Output) and follow the instructions (right-click on "Solution 'project MCML' ..." -> retarget the project).

---

## Parameter Configuration

1. On the second last line of the .mci file, you can set the parameters for the Photodiode in T and R. The 6 parameters are as follows (unit in cm):

   - Position of PD in R: (Rx, Ry); Length of side: R1;
   - Position of PD in T: (Tx, Ty); Length of side: T1.

2. On the last line, define:
   - Type of light source: 1. Point; 2. Gaussian; 3. Flat.
   - Position of the light source: (x, y) in cm.
   - Standard deviation for Gaussian: r/3 (standard value) or length for flat: d.

3. At the end of the simulation, the terminal displays:
   - Reflectance_CHatterjee and Transmittance_CHatterjee: R and T for definition 2;
   - Reflectance_SL and Transmittance SL: R and T for definition 3;
   - MOP: the mean optical path;
   - User time: the simulation time.

4. In the .mco file, you can find:
   - Diffuse reflectance and Transmittance: R and T for definition 1;
   - N_phR and N_phT: the number of photons received by the PDs in R and T.

---

## Programming Notes

- Modifications mainly in mcmlmain.c + mcmlgo.c (for photon transmission), actually everywhere (.h, etc.).
- momlio is used for display.

---

## Visualization

- Open look_mop.m with Matlab, change the name of the .mco file and run;
- get_mop.m allows users to select the output images they are interested in.

---

## Compilation

- Display --> Terminal;
- In Terminal: 'cl' followed by the file names without the .h extension (cl mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c);
- An .exe file is then generated in the same directory with the name of the first file entered in the list.

---

## Additional Explanation

You can use VS to open the .sln file and modify the MCML to be in T or R mode by changing the code:

- Commented code "Recording photon optical paths in R" can make MCML unavailable in R mode;
- Commented code "Recording photon optical paths in T" can make MCML unavailable in T mode.

---

# Projet MOP-MCML (Version Française)

Première modifications faites par Songlin Li, améliorations par Wenzheng Wang

---

## Lancement avec Visual Studio (VS) 2022

- Souci de version avec 2019 car avec VS 2019 projet .vcxproj à ouvrir. Pas possible avec version 2022.
- Trouver le .sln et l'ouvrir avec VS 2022 (ouvrir un projet ou une solution).
- Il est possible de devoir télécharger les add-ons nécessaires (15 minutes).
- Suivant les versions, message à lire dans la fenêtre en bas à gauche (Sortie) et faire ce qui est dit (clic droit sur « Solution 'projet MCML' ... » -> recible le projet).

---

## Configuration des Paramètres

1. Sur la 2ème dernière ligne du fichier .mci, on peut saisir les paramètres du Photodiode en T et en R, les 6 paramètres sont comme le suivant (unité en cm) :

   - Position du PD en R : (Rx, Ry); Longueur du côté : R1;
   - Position du PD en T : (Tx, Ty); Longueur du côté : T1.

2. Sur la dernière ligne, on définit :
   - Type de source lumineuse : 1. Point; 2. Gaussien; 3. Plat.
   - Position de la source lumineuse : (x, y) en cm.
   - Écart-type pour Gaussien : r/3 (valeur standard) ou longueur pour plat : d.

3. En fin de simulation, sur le terminal, on affiche :
   - Réflectance_CHatterjee et Transmittance_CHatterjee : la R et la T pour la définition 2;
   - Reflectance_SL et Transmittance SL : la R et la T pour la définition 3;
   - MOP : le chemin optique moyen;
   - User time : le temps de simulation.

4. On trouve dans le fichier .mco :
   - Diffuse reflectance et Transmittance : la R et la T pour la définition 1;
   - N_phR et N_phT : le nombre de photons reçus par les PD en R et en T.

---

## Remarques de Programmation

- Modifications principalement dans mcmlmain.c + mcmlgo.c (pour l'envoi des photons), partout en fait (.h, etc).
- momlio sert à l'affichage.

---

## Visualisation

- Avec Matlab ouvrir look_mop.m, changer le nom du fichier .mco et exécuter;
- get_mop.m permet aux utilisateurs de sélectionner les images de sortie qui les intéressent.

---

## Compilation

- Affichage --> Terminal;
- Dans Terminal : 'cl' puis le nom des fichiers sans le .h (cl mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c);
- Un fichier .exe est alors généré dans le même répertoire avec le nom du premier fichier rentré dans la liste.

---

## Explication Complémentaire

Vous pouvez utiliser VS pour ouvrir le fichier .sln et modifier le MCML pour qu'il soit en T ou en R en modifiant le code :

- Code commentaire « Enregistrement des chemins optiques des photons en R » est possible de rendre MCML indisponible en mode R;
- Code commentaire « Enregistrement des chemins optiques des photons en T » est possible de rendre MCML indisponible en mode T.