#include "header.h"
#include "funzioni.h"
#include "es7.h"

/* DA FARE:
-aggiustare reticoli M=2,4
-capire dimensioni scatola giuste
-potenziale esplode per caso = 1,2
*/
/* NOTE:
-se K o V esplodono molto probabilmente dt è troppo grande
-per stimare sigma t.c. T=1.1, uso N grande ()
*/

/*** variabili globali ***/
//CC, BCC, FCC
int M = 1; //1,2,4
int N = M * pow(5, 3); //numero di particelle
ofstream dati, coord, gnuplot;

int main() {
    srand(1);//default seed = 1
    int caso       = 0;//tre valori di pressione diversi
    double dt      = 0.0001;//passo temporale
    double t1      = 0.0002;//durata simulazione

    double sigma[3][3] = {{0.951309, 1, 1},//CC
        {0.951309, 1, 1},//BCC
        {0.951309, 1, 1}//FCC
    };

    int animazione = 1;//0: grafici; 1: particelle
    double pausa   = 0.5;
    int N_step     = 2;
    double scala_v = 0.2;
    //###################################################
    struct vec r[N], v[N], a[N];
    int reticolo   = log2(M);

    double T_req = 1.1;
    double rho[] = {0.01, 0.8, 1.2};

    double L = cbrt(N / rho[caso]);
    double r_c = L / 2;

    double t = 0;
    int N_t = (t1 - t) / dt;
    if (N_step > N_t) {N_step = N_t;}
    int skip = rint(N_t / N_step); if (skip == 0) {skip = 1;}

    crea_reticolo(r, L);
    crea_reticolo(a, 0);
    distr_gauss(v, sigma[reticolo][caso]);

    vec *dr[N]; vec_2D(dr, N);
    a_LJ(r, a, dr, r_c, L);

    dati.open("dati.dat");
    coord.open("coordinate.xyz");
    gnuplot.open("gnuplot.dat");

    gnuplot << animazione << "\n" << N << "\n" << N_step << "\n" << L << "\n" << pausa << "\n" << skip << "\n" << dt << endl;
    coord << N << endl;
    stampa_costanti(M, N, caso, rho, L, t1, dt, N_t, N_step, skip, pausa);

    double E = 0, T = 0, P = 0;
    double K = 0, V = 0, W = 0;
    double K_c = 0, V_c = 0, W_c = 0;

    for (int i = 0; i < N_t; ++i) {//tempo
        K_c = K_c * (i + 1.0);
        V_c = V_c * (i + 1.0);
        W_c = W_c * (i + 1.0);

        if (animazione == 1) {
            for (int p = 0; p < N; ++p) {
                coord << "P" << p << "\t" << r[p].x << "\t" << r[p].y << "\t" << r[p].z << "\t" << v[p].x*scala_v << "\t" << v[p].y*scala_v << "\t" << v[p].z*scala_v << endl;
            }
        }
        vel_verlet(r, v, a, dt, r_c, L, &K, &V, &W);
        v_cm_0(v);
        K_c = (K_c + K) / (i + 2.0);
        V_c = (V_c + V) / (i + 2.0);
        W_c = (W_c + W) / (i + 2.0);

        T = 2.0 * K_c / (3.0 * N);
        P = rho[caso] * (1 + W_c / (3.0 * T_req)); //P su T_req

        E = K + V;

        dati << t << "\t" << K << "\t" << V << "\t" << E << "\t" << T << "\t" << P << endl;
        t += dt;
    }

    cout << "T = " << T << "\t\tNuovo sigma = " << sigma[reticolo][caso]*sqrt(T_req / T) << endl;
    dati.close();
    coord.close();
    gnuplot.close();



    return 0;
}

