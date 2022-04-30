#include "header.h"
#include "funzioni.h"

/*** variabili globali ***/
int N = pow(n, 3); //numero di particelle
ofstream dati;

int main() {
    struct vec r[N], v[N], F[N], *dr[N];
    vec_2D(dr);//crea "vec" bidimensionale

    int caso = 0;
    double r_c = 3.0;
    double rho[] = {0.01, 0.7, 1.2};
    double sigma[] = {1.1, 1.05, 1.05};
    double L = cbrt(N / rho[caso]); //larghezza lato

    double t = 0, t1 = 1;
    double dt = 0.1;
    int N_t = (t1 - t) / dt;
    //crea_reticolo(r, L); //modifica r

    crea_reticolo(r, L);
    crea_reticolo(v, 0);

    //distr_gauss(v, sigma[caso]);
    //vel_media_0(v);

    F_LJ(r, dr, F, L);//calcola forze F[3N]

    // for (int i = 0; i < N; ++i) {
    //     F[i].x=

    // }
    stampa_coord(dr[3]);

    double K = 0, V = 0, T;
    for (int i = 0; i < N_t; ++i) {//tempo
        K = 0;

        //K *= i + 1;

        vel_verlet(t, r, dr, v, F, dt, L, &K, &V);

        //K /= i + 2;
        //LOG(K)

        for (int i = 0; i < N; ++i) {
            double v_mod = v[i].mod();
            K += 0.5 * v_mod * v_mod;
        }

        T = 2.0 * K / (3.0 * N);
        LOG(T)
    }

    return 0;
}


