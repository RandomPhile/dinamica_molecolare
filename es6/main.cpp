#include "header.h"
#include "funzioni.h"

/*** variabili globali ***/
int N = 1; //numero di particelle
ofstream dati;

int main() {
    struct vec r[N], v[N], F[N];

    double t = 0, t1 = 100;
    double dt = 0.1;
    int N_t = (t1 - t) / dt;

    for (int i = 0; i < N; ++i) {
        r[i].uguale(0);
    }

    distr_gauss(v, 1.1);
    //vel_media_0(v);

    F_osc(r, F);//calcola forze

    //stampa_coord(r);

    dati.open("dati.dat");
    double K = 0, V = 0, T;
    for (int i = 0; i < N_t; ++i) {//tempo
        //K *= i + 1;

        vel_verlet(t, r, v, F, dt, &K, &V);

        //K /= i + 2;
        //cout << t << "\t" << K + V << "\t" << K << "\t" << V << endl;
        dati << t << "\t" << K + V << "\t" << K << "\t" << V << endl;
        T = 2.0 * K / (3.0 * N);
        //LOG(T)
        t += dt;
    }
    dati.close();
    return 0;
}


