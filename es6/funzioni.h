#include "header.h"
struct vec {
    double x; double y; double z;

    double mod() {
        return sqrt(x * x + y * y + z * z);
    }
    void uguale(double val) {
        x = val; y = val; z = val;
    }
    void piu(double val) {
        x += val; y += val; z += val;
    }
    void per(double val) {
        x *= val; y *= val; z *= val;
    }
};
void vec_2D (vec *dr[N]) {
    for (int i = 0; i < N; ++i)
        dr[i] = new vec[N];
}
// void stampa_stato(vec *r, vec *v, vec *F){
//     for (int i = 0; i < N; ++i) {
//         cout << "i = " << i << endl;
//         cout << "rx = " << r[i].x << endl;
//         cout << "vx = " << v[i].x << endl;
//         cout << "fx = " << F[i].x << "\n" << endl;

//         cout << "ry = " << r[i].y << endl;
//         cout << "vy = " << v[i].y << endl;
//         cout << "fy = " << F[i].y << "\n" << endl;

//         cout << "rz = " << r[i].z << endl;
//         cout << "vz = " << v[i].z << endl;
//         cout << "fz = " << F[i].z << "\n" << endl;
//         cout << "\n" << endl;
//     }
// }
void stampa_coord(vec *r) {
    for (int i = 0; i < N; ++i) {
        cout << r[i].x << "\t" << r[i].y << "\t" << r[i].z << "\t" << endl;
    }
    LOG("\n\n");
}
void distr_gauss(vec *x, double sigma) {
    for (int i = 0; i < N; ++i) {
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        x[i].x = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        x[i].y = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        x[i].z = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
    }
}
void vel_media_0(vec *v) {
    vec v_cm = {.x = 0, .y = 0, .z = 0};
    for (int i = 0; i < N; ++i) {
        v_cm.x += v[i].x;
        v_cm.y += v[i].y;
        v_cm.z += v[i].z;
    }
    v_cm.per(1 / N);
    for (int i = 0; i < N; ++i) {
        v[i].x -= v_cm.x;
        v[i].y -= v_cm.y;
        v[i].z -= v_cm.z;
    }
}
void F_osc(vec *r, vec *F) {
    for (int i = 0; i < N; ++i) {
        F[i].uguale(0);
    }
    for (int i = 0; i < N; ++i) {//particella i=1..N
        double r_mod = r[i].mod();

        F[i].x = -r[i].x;
        F[i].x = -r[i].x;
        F[i].x = -r[i].x;
    }
}
void vel_verlet(double t, vec *r, vec *v, vec *F, double dt, double *K, double *V) {
    vec F_old[N]; double v_mod; double r_mod;
    *K = 0; *V = 0;

    for (int i = 0; i < N; ++i) {
        //velocità
        v[i].x += 0.5 * dt * (F[i].x);
        v[i].y += 0.5 * dt * (F[i].y);
        v[i].z += 0.5 * dt * (F[i].z);

        //posizione
        r[i].x += v[i].x * dt;
        r[i].y += v[i].y * dt;
        r[i].z += v[i].z * dt;
    }
    //forze
    //F_osc(r, F);
    for (int i = 0; i < N; ++i) {//particella i=1..N
        double r_mod = r[i].mod();

        F[i].x = -r[i].x;
        F[i].x = -r[i].x;
        F[i].x = -r[i].x;
    }
    for (int i = 0; i < N; ++i) {
        //velocità
        v[i].x += 0.5 * dt * F[i].x;
        v[i].y += 0.5 * dt * F[i].y;
        v[i].z += 0.5 * dt * F[i].z;

        //osservabili fisiche
        r_mod = r[i].mod();
        v_mod = v[i].mod();
        *V += 0.5 * r_mod * r_mod;
        *K += 0.5 * v_mod * v_mod;
    }
}
void vel_verlet2(double t, vec *r, vec *v, vec *F, double dt, double *K, double *V) {
    vec F_old[N]; double v_mod; double r_mod;
    *K = 0; *V = 0;

    for (int i = 0; i < N; ++i) {
        //salvo la forza corrente
        F_old[i] = F[i];
        //posizione
        r[i].x += v[i].x * dt + 0.5 * F_old[i].x * dt * dt;
        r[i].y += v[i].y * dt + 0.5 * F_old[i].y * dt * dt;
        r[i].z += v[i].z * dt + 0.5 * F_old[i].z * dt * dt;
    }
    //forze
    F_osc(r, F);
    for (int i = 0; i < N; ++i) {
        //velocità
        v[i].x += 0.5 * dt * (F_old[i].x + F[i].x);
        v[i].y += 0.5 * dt * (F_old[i].y + F[i].y);
        v[i].z += 0.5 * dt * (F_old[i].z + F[i].z);

        //osservabili fisiche
        r_mod = r[i].mod();
        v_mod = v[i].mod();
        *V += 0.5 * r_mod * r_mod;
        *K += 0.5 * v_mod * v_mod;
    }
}