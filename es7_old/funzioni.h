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
void crea_reticolo(vec *r, double L) {
    int n = cbrt(N);
    float L_cella = L / n;
    int cont = 0;//contatore coordinata particella (fino a 3N)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                r[cont].x = i * L_cella;
                r[cont].y = j * L_cella;
                r[cont].z = k * L_cella;
                cont++;;
            }
        }
    }
}
void stampa_stato(vec *r, vec *v, vec *F){
    for (int i = 0; i < N; ++i) {
        cout << "i = " << i << endl;
        cout << "rx = " << r[i].x << endl;
        cout << "vx = " << v[i].x << endl;
        cout << "fx = " << F[i].x << "\n" << endl;

        cout << "ry = " << r[i].y << endl;
        cout << "vy = " << v[i].y << endl;
        cout << "fy = " << F[i].y << "\n" << endl;

        cout << "rz = " << r[i].z << endl;
        cout << "vz = " << v[i].z << endl;
        cout << "fz = " << F[i].z << "\n" << endl;
        cout << "\n" << endl;
    }
}
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
void cerca_vicini(vec *r, vec **dr, double L) {
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            dr[i][j].x = r[i].x - r[j].x;
            dr[i][j].y = r[i].y - r[j].y;
            dr[i][j].z = r[i].z - r[j].z;

            dr[i][j].x -= L * rint(dr[i][j].x / L);
            dr[i][j].y -= L * rint(dr[i][j].y / L);
            dr[i][j].z -= L * rint(dr[i][j].z / L);

            dr[j][i].x = -dr[i][j].x;
            dr[j][i].y = -dr[i][j].y;
            dr[j][i].z = -dr[i][j].z;
        }
    }
}
void F_LJ(vec *r, vec *dr[N], vec *F, double L) {
    for (int i = 0; i < N; ++i) {
        F[i].uguale(0);
    }
    cerca_vicini(r, dr, L);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; j++) {
            double r_mod = dr[i][j].mod();

            if (r_mod < 3.0) {//r_c
                //double cost = 24 * (pow(r_mod, -8) - 2 * pow(r_mod, -14));
                double cost = 24 * 2 * pow(r_mod, -14);

                F[i].x += dr[i][j].x * cost;
                F[i].y += dr[i][j].y * cost;
                F[i].z += dr[i][j].z * cost;

                F[j].x += -dr[i][j].x * cost;
                F[j].y += -dr[i][j].y * cost;
                F[j].z += -dr[i][j].z * cost;
            }
        }
    }
}
void vel_verlet(double t, vec *r, vec *dr[], vec *v, vec *F, double dt, double L, double *K, double *V) {
    vec F_old[N]; double v_mod;

    for (int i = 0; i < N; ++i) {
        v[i].x += 0.5 * dt * F[i].x;
        v[i].y += 0.5 * dt * F[i].y;
        v[i].z += 0.5 * dt * F[i].z;

        r[i].x += v[i].x * dt;
        r[i].y += v[i].y * dt;
        r[i].z += v[i].z * dt;

        //sposto dentro alla scatola
        r[i].x -= L * rint(r[i].x / L);
        r[i].y -= L * rint(r[i].y / L);
        r[i].z -= L * rint(r[i].z / L);
    }
    F_LJ(r, dr, F, L);
    for (int i = 0; i < N; ++i) {
        v[i].x += 0.5 * dt * F[i].x;
        v[i].y += 0.5 * dt * F[i].y;
        v[i].z += 0.5 * dt * F[i].z;

    }
    t+=dt;
}