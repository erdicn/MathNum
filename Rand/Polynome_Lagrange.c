
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 4

int read_dat(double *x, double *f)
{

    FILE *fichier = fopen("données.dat", "rt");
    if (fichier == NULL)
    {
        printf("Problème d'ouverture du ficheir \n");
        exit(0);
    }

    else
    {
        printf("x f\n");
        for (int i = 0; i < MAX; i++)
        {
            fscanf(fichier, "%lf %lf", &x[i], &f[i]);
            printf("%lf %lf\n", x[i], f[i]);
        }
    }
    return 1;
}

void Polynome_Lagrange(double *zk, double *x, double *f, int n, int m, double *Pk)
{
    double L;
    double hh = (x[MAX - 1] - x[0]) / m;
    // printf("%lf\n",h);
    for (int k = 0; k <= m; k++)
    {
        zk[k] = x[0] + k * hh;
    }

    // printf("\n");

    for (int test = 0; test <= m; test++)
    {
        Pk[test] = 0;
        for (int i = 0; i < n; i++)
        {
            L = 1.0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {

                    L = L * (zk[test] - x[j]) / (x[i] - x[j]);
                }
            }
            Pk[test] = Pk[test] + f[i] * L;
        }
    }
}

void write_dat(double *zk, double *Pk, int m)
{
    // cette fonction va créer un fichier .dat et va éccire dans ce fichier les valeurs de zk et Pk

    FILE *fichier2 = NULL;
    fichier2 = fopen("Values.dat", "w");
    if (fichier2 == NULL)
    {
        printf("Problème d'ouverture du ficheir \n");
        exit(0);
    }

    for (int i = 0; i <= m; i++)
    {
        fprintf(fichier2, "%lf %lf\n", zk[i], Pk[i]);
        printf("%lf %lf\n", zk[i], Pk[i]);
    }

    fclose(fichier2);
}

double Integrale_Trapeze_A(double *Pk, float a, float b, int m)
{

    double Ia;
    double h = (b - a) / m;
    printf("%f\n", h);
    for (int i = 0; i <= m - 1; i++)
    {
        Ia = Ia + h * ((Pk[i] + Pk[i + 1]) / 2);
    }
    return Ia;
}
double Integrale_Trapeze_E()
{

    double Ie = 448.0/9.0;//5. * 6. + 5. / 8. * pow(6, 2) + 1. / 72. * pow(6, 3) - 1. / 96. * pow(6, 4) - (5. * (-2.) + 5. / 8. * pow(-2, 2) + 1. / 72. * pow(-2, 3) - 1. / 96. * pow(-2, 4));
    return Ie;
}

void Erreur(float a, float b, double *H, double *Pk, double Ia, double Ie, int n, double *zk, double *f, double *x, double *E)
{

    int m = 80;
    double IaE;
    double L;

    for (int j = 0; j <= 4; j++)
    {
        int a = 0;
        Polynome_Lagrange(zk, x, f, n, m, Pk);

        for (int i = 0; i <= m; i++)
        {
            IaE = IaE + H[j] * ((Pk[i] + Pk[i + 1]) / 2);
        }

        printf("IaE : %lf\n", IaE);
        E[j] = fabs(Ie - IaE);
        IaE = 0;
        m = m * 10;

        printf("H : %lf\n", H[j]);

        printf("erreur : %18.10e\n", E[j]);
    }
}

void write_erreur_dat(double *H, double *E)
{
    // ecrit les valeurs de h et les erreurs associées

    FILE *fichier3 = NULL;
    fichier3 = fopen("Erreur.dat", "w");
    if (fichier3 == NULL)
    {
        printf("Problème d'ouverture du fichier \n");
        exit(0);
    }

    for (int i = 0; i < 5; i++)
    {
        fprintf(fichier3, "%.16lf %.16lf\n", H[i], E[i]);
    }

    fclose(fichier3);
}

int main()

{
    int m = 10;
    // borne de l'intégrale
    int a = -2.0;
    int b = 6.0;

    // faire malloc
    double *x = NULL;
    x = malloc(sizeof(double) * MAX);
    double *f = NULL;
    f = malloc(sizeof(double) * MAX);
    double *zk = NULL;
    zk = malloc(sizeof(double) * m);
    // double* Pk=NULL;
    // Pk=malloc(sizeof(double)*m);
    int n = MAX;
    double Pk[1000000];
    double Ia;
    double Ie;
    double E[5];
    double H[5] = {0.1, 0.01, 0.001, 0.0001, 0.00001}; // contient les valeurs de h.

    // int x[4]={-2,0,4,6};
    // int f[4]={3,5,8,5};
    read_dat(x, f);

    // printf("%d\n",x[3]); //test ; les valeurs de x et f sont biens défini
    Polynome_Lagrange(zk, x, f, n, m, Pk);
    for (int k = 0; k <= m; k++)
    {
        printf("zk%d : %lf\n", k, zk[k]);
    }
    for (int k = 0; k <= m; k++)
    {

        printf("Pk%d : %lf\n", k, Pk[k]);
    }
    write_dat(zk, Pk, m);

    Ia = Integrale_Trapeze_A(Pk, a, b, m);
    printf("La valeur approche de l'Intégrale du Polynome de Lagrange de [-2,6] vaut : %lf\n", Ia);

    Ie = Integrale_Trapeze_E();
    printf("La valeur exacte de l'Intégrale du Polynome de Lagrange de [-2,6] vaut : %lf\n", Ie);

    printf("L'erreur est égale à %lf\n", Ie - Ia);

    // Pk=malloc(sizeof(double)*800000);
    zk = malloc(sizeof(double) * 800000); // on change la taille pour tout mettre
    Erreur(a, b, H, Pk, Ia, Ie, n, zk, f, x, E);

    write_erreur_dat(H, E);

    return 0;
}
