#include "header.h"
#include "funzioni.h"

/*** variabili globali ***/
//int M = 1;//tipo di reticolo
int n = 2;//numero di particelle per lato
int N = pow(n, 3); //numero di particelle
ofstream dati;

int main() {
    struct vec r[N], v[N], F[N], *dr[N];
    vec_2D(dr);//crea "vec" bidimensionale

    int caso = 0;

    double rho[] = {0.01, 0.8, 1.2};
    double sigma[] = {1.1, 1.05, 1.05};
    double L = cbrt(N / rho[caso]); //larghezza lato

    double t = 0, t1 = 1;
    double dt = 0.1;
    int N_t = (t1 - t) / dt;
    crea_reticolo(r, L); //modifica r
    //r[0].uguale(0.1);

    distr_gauss(v, sigma[caso]);
    vel_media_0(v);

    F_LJ(r, dr, F, L);//calcola forze F[3N]

    double K = 0, V = 0, T;
    for (int i = 0; i < N_t; ++i) {//tempo
        K = 0;
        t += dt;

        //K *= i + 1;

        vel_verlet(t, r, dr, v, F, dt, L, &K, &V);

        //K /= i + 2;
        //LOG(K)

        T = 2.0 * K / (3.0 * N);
        LOG(T)
    }


    return 0;
}


