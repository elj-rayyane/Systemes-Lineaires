# ------------------------------------------------------------------
# Résolution d'un système linéaire Ax = b par décomposition LU (méthode de Doolittle)
#
# Principe :
#   1. Décomposer A en L * U  (L = triangulaire inférieure avec 1 sur la diagonale,
#                               U = triangulaire supérieure)
#   2. Résoudre L * y = b  par substitution avant
#   3. Résoudre U * x = y  par substitution arrière
# ------------------------------------------------------------------

# ---- Saisie de la taille ----
n = int(input("Taille du systeme (n) : "))

# Déclaration et initialisation des structures (tableaux à base de listes)
A = [[0.0] * n for _ in range(n)]
b = [0.0] * n
L = [[0.0] * n for _ in range(n)]
U = [[0.0] * n for _ in range(n)]
y = [0.0] * n
x = [0.0] * n

# ---- Saisie de la matrice A ----
print(f"Entrez la matrice A ({n} x {n}) :")
for i in range(n):
    for j in range(n):
        A[i][j] = float(input(f"  A[{i + 1}][{j + 1}] = "))

# ---- Saisie du vecteur second membre b ----
print(f"Entrez le vecteur b ({n} valeurs) :")
for i in range(n):
    b[i] = float(input(f"  b[{i + 1}] = "))

# ------------------------------------------------------------------
# ETAPE 1 : Initialisation de L et U
# L est initialisée à l'identité (diagonale = 1, reste = 0)
# U est initialisée à zéros
# ------------------------------------------------------------------
for i in range(n):
    for j in range(n):
        L[i][j] = 0.0
        U[i][j] = 0.0
    L[i][i] = 1.0

# ------------------------------------------------------------------
# ETAPE 2 : Décomposition LU — algorithme de Doolittle
#
# Pour chaque pivot k (de 0 à n-1) :
#   - On calcule la ligne k de U : U[k][j] pour j allant de k à n-1
#   - On calcule la colonne k de L : L[i][k] pour i allant de k+1 à n-1
# ------------------------------------------------------------------
for k in range(n):

    # Calcul de la ligne k de U
    for j in range(k, n):
        somme = 0.0
        for i in range(k):
            somme += L[k][i] * U[i][j]
        U[k][j] = A[k][j] - somme

    # Calcul de la colonne k de L (en dessous du pivot)
    for i in range(k + 1, n):
        somme = 0.0
        for j in range(k):
            somme += L[i][j] * U[j][k]
        L[i][k] = (A[i][k] - somme) / U[k][k]

# Affichage de L et U
print("\nMatrice L :")
for i in range(n):
    for j in range(n):
        print(f"{L[i][j]:8.4f}", end=" ")
    print()

print("\nMatrice U :")
for i in range(n):
    for j in range(n):
        print(f"{U[i][j]:8.4f}", end=" ")
    print()

# ------------------------------------------------------------------
# ETAPE 3 : Substitution avant — résolution de L * y = b
#
# y[i] = b[i] - somme(L[i][j] * y[j]) pour j de 0 à i-1
# ------------------------------------------------------------------
for i in range(n):
    somme = 0.0
    for j in range(i):
        somme += L[i][j] * y[j]
    y[i] = b[i] - somme

print("\nVecteur y (solution de Ly = b) :")
for i in range(n):
    print(f"  y[{i + 1}] = {y[i]:8.4f}")

# ------------------------------------------------------------------
# ETAPE 4 : Substitution arrière — résolution de U * x = y
#
# x[i] = (y[i] - somme(U[i][j] * x[j])) / U[i][i]  pour j de i+1 à n-1
# ------------------------------------------------------------------
for i in range(n - 1, -1, -1):
    somme = 0.0
    for j in range(i + 1, n):
        somme += U[i][j] * x[j]
    x[i] = (y[i] - somme) / U[i][i]

print("\nSolution x :")
for i in range(n):
    print(f"  x[{i + 1}] = {x[i]:8.4f}")

# --- Vérification ---
print("\nVérification par substitution :")
print("En remplaçant les inconnues par les valeurs trouvées, on obtient :")
for i in range(n):
    somme_verif = 0.0
    ligne = ""
    for j in range(n):                          
        somme_verif += A[i][j] * x[j]
        ligne += f"{A[i][j]:.2f} * {x[j]:.2f}"
        if j < n - 1:
            ligne += " + "
    ligne += f" = {somme_verif:.2f} (cible : {b[i]:.2f})"
    print(ligne) 