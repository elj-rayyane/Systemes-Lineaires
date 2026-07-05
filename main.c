/*
 * Résolution d'un système linéaire Ax = b par décomposition LU (méthode de Doolittle)
 *
 * Principe :
 *   1. Décomposer A en L * U  (L = triangulaire inférieure avec 1 sur la diagonale,
 *                               U = triangulaire supérieure)
 *   2. Résoudre L * y = b  par substitution avant
 *   3. Résoudre U * x = y  par substitution arrière
 *
 * Le programme fonctionne pour un système de taille quelconque (jusqu'à NMAX).
 */

#include <stdio.h>

#define NMAX 10   /* taille maximale autorisée pour le système */

int main(){
    /* ------------------------------------------------------------------
     * Déclaration des tableaux
     *
     * A[i][j]  : matrice carrée du système, n lignes et n colonnes
     * L[i][j]  : matrice triangulaire inférieure (1 sur la diagonale)
     *            L[i][j] = 0 si j > i, L[i][i] = 1
     * U[i][j]  : matrice triangulaire supérieure
     *            U[i][j] = 0 si j < i
     * b[i]     : vecteur second membre, n valeurs
     * y[i]     : vecteur intermédiaire, solution de L*y = b, n valeurs
     * x[i]     : vecteur solution finale, solution de U*x = y, n valeurs
     * ------------------------------------------------------------------ */
    double A[NMAX][NMAX];
    double L[NMAX][NMAX];
    double U[NMAX][NMAX];
    double b[NMAX];
    double y[NMAX];
    double x[NMAX];

    int n;          /* taille effective du système */
    int i, j, k;   /* indices de boucles */
    double somme;   /* accumulateur temporaire pour les produits scalaires */
    double somme_verif; /* accumulateur des coefficients.inconnues dans la vérification */

    /* ---- Saisie de la taille ---- */
    printf("Taille du systeme (n) : ");
    scanf("%d", &n);

    /* ---- Saisie de la matrice A ---- */
    printf("Entrez la matrice A (%d x %d) :\n", n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("  A[%d][%d] = ", i + 1, j + 1);
            scanf("%lf", &A[i][j]);
        }
    }

    /* ---- Saisie du vecteur second membre b ---- */
    printf("Entrez le vecteur b (%d valeurs) :\n", n);
    for (i = 0; i < n; i++) {
        printf("  b[%d] = ", i + 1);
        scanf("%lf", &b[i]);
    }

    /* ------------------------------------------------------------------
     * ETAPE 1 : Initialisation de L et U
     * L est initialisée à l'identité (diagonale = 1, reste = 0)
     * U est initialisée à zéros
     * Boucle sur i et j de 0 à n-1 ; condition de sortie : fin de [0, n-1]
     * ------------------------------------------------------------------ */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;
        }
        L[i][i] = 1.0;
    }

    /* ------------------------------------------------------------------
     * ETAPE 2 : Décomposition LU — algorithme de Doolittle
     *
     * Pour chaque pivot k (de 0 à n-1) :
     *   - On calcule la ligne k de U : U[k][j] pour j allant de k à n-1
     *   - On calcule la colonne k de L : L[i][k] pour i allant de k+1 à n-1
     *
     * Boucle externe sur k : condition de sortie k = n
     * Boucle interne sur j (ligne de U) : condition de sortie j = n
     * Boucle interne sur i (colonne de L) : condition de sortie i = n
     * ------------------------------------------------------------------ */
    for (k = 0; k < n; k++) {

        /* Calcul de la ligne k de U */
        for (j = k; j < n; j++) {
            somme = 0.0;
            for (i = 0; i < k; i++) {
                somme += L[k][i] * U[i][j];
            }
            U[k][j] = A[k][j] - somme;
        }

        /* Calcul de la colonne k de L (en dessous du pivot) */
        for (i = k + 1; i < n; i++) {
            somme = 0.0;
            for (j = 0; j < k; j++) {
                somme += L[i][j] * U[j][k];
            }
            L[i][k] = (A[i][k] - somme) / U[k][k];
        }
    }

    /* Affichage de L et U */
    printf("\nMatrice L :\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%8.4lf ", L[i][j]);
        }
        printf("\n");
    }

    printf("\nMatrice U :\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%8.4lf ", U[i][j]);
        }
        printf("\n");
    }

    /* ------------------------------------------------------------------
     * ETAPE 3 : Substitution avant — résolution de L * y = b
     *
     * y[i] = b[i] - somme(L[i][j] * y[j]) pour j de 0 à i-1
     * (L[i][i] = 1 donc pas de division)
     *
     * Boucle sur i de 0 à n-1 ; condition de sortie : i = n
     * Boucle interne sur j de 0 à i-1 ; condition de sortie : j = i
     * ------------------------------------------------------------------ */
    for (i = 0; i < n; i++) {
        somme = 0.0;
        for (j = 0; j < i; j++) {
            somme += L[i][j] * y[j];
        }
        y[i] = b[i] - somme;
    }

    printf("\nVecteur y (solution de Ly = b) :\n");
    for (i = 0; i < n; i++) {
        printf("  y[%d] = %8.4lf\n", i + 1, y[i]);
    }

    /* ------------------------------------------------------------------
     * ETAPE 4 : Substitution arrière — résolution de U * x = y
     *
     * x[i] = (y[i] - somme(U[i][j] * x[j])) / U[i][i]  pour j de i+1 à n-1
     *
     * Boucle sur i de n-1 à 0 (sens décroissant) ; condition de sortie : i < 0
     * Boucle interne sur j de i+1 à n-1 ; condition de sortie : j = n
     * ------------------------------------------------------------------ */
    for (i = n - 1; i >= 0; i--) {
        somme = 0.0;
        for (j = i + 1; j < n; j++) {
            somme += U[i][j] * x[j];
        }
        x[i] = (y[i] - somme) / U[i][i];
    }

    printf("\nSolution x :\n");
    for (i = 0; i < n; i++) {
        printf("  x[%d] = %8.4lf\n", i + 1, x[i]);
    }


    /*--- Vérification ---*/
    printf("\nVérification par substitution : \n");
    printf("En remplaçant les inconnues par les valeurs trouvées, on obtient : \n");
    for(i = 0; i < n; i++) {
        somme_verif = 0.0;
        for(j = 0; j < n; j++) {
            printf("%.4lf * %.4lf ", A[i][j], x[j]);
            somme_verif += A[i][j] * x[j];
            /* --- affichage esthétiquement cohérent ---*/
            if(j < n-1){
                printf(" + ");
            }
        }
        
        printf(" = %.4lf (cible : %.4lf)\n", somme_verif, b[i]);
    }

    return 0;
}