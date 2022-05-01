#include "header.h"
#include "funzioni.h"
#include "es7.h"
/* DA FARE:
-aggiustare reticoli M=2,4
-capire dimensioni scatola giuste
*/


/*** variabili globali ***/
//CC, BCC, Fcc
int M = 1; //1,2,4
int N = M * pow(2, 3); //numero di particelle
ofstream dati, coord, vel, gnuplot;

int main() {
    srand(4);//default seed = 1
    struct vec r[N], v[N], a[N];

    int caso = 0;
    double rho[] = {0.01, 0.8, 1.2};
    double sigma[] = {1.03, 1.12, 1};
    double L = cbrt(N / rho[caso]);
    double r_c = L / 2; //

    double pausa = 0.5;
    int N_step = 1;
    
    double t = 0, t1 = 1;
    double dt = 0.001;
    int N_t = (t1 - t) / dt;
    int skip = rint(N_t / N_step);if (skip==0){skip=1;}

    crea_reticolo(r, L);
    //r[0].uguale(1);
    crea_reticolo(v, 0);//v_i nulle
    crea_reticolo(a, 0);

    a_LJ(r, a, r_c, L);

    distr_gauss(v, sigma[caso]);
    // for (int i = 0; i < N; ++i) {
    //     v[i].per(fattore[caso]);
    // }
    //v_cm_0(v);
    //stampa_stato(r,v,a);

    dati.open("dati.dat");
    coord.open("coordinate.xyz");
    vel.open("vel.xyz");
    gnuplot.open("gnuplot.dat");

    gnuplot << N << "\n" << N_step << "\n" << L << "\n" << pausa << "\n" << skip << "\n" << dt << endl;
    
    cout<<"M      = "<<M<<endl;
    cout<<"N      = "<<N<<endl;
    cout<<"caso   = "<<caso<<endl;
    cout<<"rho    = "<<rho[caso]<<endl;
    cout<<"L      = "<<L<<endl;
    cout<<"t1     = "<<t1<<endl;
    cout<<"dt     = "<<dt<<endl;
    cout<<"N_t    = "<<N_t<<endl;
    cout<<"N_step = "<<N_step<<endl;
    cout<<"skip   = "<<skip<<endl;
    cout<<"pausa  = "<<pausa<<endl;
    LOG("\n\n\n");


    coord << N << endl;
    vel << N << endl;

    double K = 0, V = 0, E = 0, T = 0;
    double K_c = 0, V_c = 0;

    for (int i = 0; i < N_t; ++i) {//tempo
        K_c = K_c * (i + 1.0);
        V_c = V_c * (i + 1.0);
        for (int p = 0; p < N; ++p) {
            coord << "P" << p << "\t" << r[p].x << "\t" << r[p].y << "\t" << r[p].z << endl;
            vel << "P" << p << "\t" << v[p].x << "\t" << v[p].y << "\t" << v[p].z << endl;
        }
        // coord << "\n\n" << endl;
        // vel << "\n\n" << endl;
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
    vel.close();
    gnuplot.close();



    return 0;
}


