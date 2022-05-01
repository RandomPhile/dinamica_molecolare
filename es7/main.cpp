#include "header.h"
#include "funzioni.h"
#include "es7.h"

/* DA FARE:
-aggiustare reticoli M=2,4
-capire dimensioni scatola giuste
*/


/*** variabili globali ***/
//CC, BCC, Fcc
int M = 2; //1,2,4
int N = M * pow(2, 3); //numero di particelle
ofstream dati, coord, gnuplot;

int main() {
    int reticolo = 0;
    srand(2);//default seed = 1
    struct vec r[N], v[N], a[N];

    int caso = 0;
    double rho[] = {0.01, 0.8, 1.2};
    double sigma[] = {1.03, 1.12, 1};
    double L = cbrt(N / rho[caso]);
    double r_c = L / 2; //

    double t = 0, t1 = 2;
    double dt = 0.001;
    int N_t = (t1 - t) / dt;

    double pausa = 0.001;
    int N_step = 200;
    if (N_step > N_t) {N_step = N_t;}
    int skip = rint(N_t / N_step); if (skip == 0) {skip = 1;}

    crea_reticolo(r, L);
    //r[0].uguale(1);
    // crea_reticolo(v, 0);//v_i nulle
    crea_reticolo(a, 0);

    a_LJ(r, a, r_c, L);

    distr_gauss(v, sigma[caso]);
    //stampa_stato(r,v,a);

    dati.open("dati.dat");
    coord.open("coordinate.xyz");
    gnuplot.open("gnuplot.dat");

    gnuplot << reticolo << "\n" << N << "\n" << N_step << "\n" << L << "\n" << pausa << "\n" << skip << "\n" << dt << endl;
    coord << N << endl;
    stampa_costanti(M, N, caso, rho, L, t1, dt, N_t, N_step, skip, pausa);

    double K = 0, V = 0, E = 0, T = 0;
    double K_c = 0, V_c = 0;

    for (int i = 0; i < N_t; ++i) {//tempo
        K_c = K_c * (i + 1.0);
        V_c = V_c * (i + 1.0);

        if (reticolo == 1) {
            for (int p = 0; p < N; ++p) {
                coord << "P" << p << "\t" << r[p].x << "\t" << r[p].y << "\t" << r[p].z << "\t" << v[p].x << "\t" << v[p].y << "\t" << v[p].z << endl;
            }
        }
        vel_verlet(r, v, a, dt, r_c, L, &K, &V);
        v_cm_0(v);
        K_c = (K_c + K) / (i + 2.0);
        V_c = (V_c + V) / (i + 2.0);
        //stampa_stato(r,v,a);
        T = 2.0 * K_c / (3.0 * N);
        E = K_c + V_c;
        //LOG(V)

        dati << t << "\t" << K_c << "\t" << V_c << "\t" << E << "\t" << T << endl;
        t += dt;
    }
    dati.close();
    coord.close();
    gnuplot.close();



    return 0;
}

