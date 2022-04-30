#include "header.h"
#include "funzioni.h"
#include "es7.h"

/*** variabili globali ***/
int M = 1; //1,2,4
int N = M*pow(2, 3); //numero di particelle
ofstream dati;

int main() {
    srand(4);//default seed = 1
    struct vec r[N], v[N], a[N];

    int caso = 0;
    double r_c = 3;//
    double rho[] = {0.01, 0.7, 1.2};
    double sigma[] = {1, 1.05, 1.05};
    double L = cbrt(N / rho[caso]);

    double t = 0, t1 = 3;
    double dt = 0.1;
    int N_t = (t1 - t) / dt;


    crea_reticolo(r, L);
    r[0].uguale(0);
    crea_reticolo(v, 0);//v_i nulle

    a_LJ(r,a,r_c,L);

    distr_gauss(v, sigma[caso]);
    v_cm_0(v);
    //stampa_stato(r,v,a);

    dati.open("dati.dat");
    double K = 0, V = 0, E = 0, T = 0;
    double K_c = 0;

    for (int i = 0; i < N_t; ++i) {//tempo
        K_c = K_c * (i + 1.0);

        vel_verlet(r, v, a, dt, r_c, L, &K, &V);

        K_c = (K_c + K) / (i + 2.0);
        //stampa_stato(r,v,a);
        T = 2.0 * K_c / (3.0 * N);
        E = K + V;
        dati << t << "\t" << K << "\t" << V << "\t" << E << "\t" << T << endl;
        t += dt;
    }

    dati.close();
    return 0;
}


