#include "header.h"
void crea_reticolo(vec *r, double L) {
    int n = cbrt(N/M);
    float a = L / cbrt(N / M);
    int cont = 0;//contatore coordinata particella (fino a 3N)

    double b[][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    if (M==2) {
        LOG(b[1][2])
        b[1][2] = 0.5;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                vec R = {.x = i * a, .y = j * a, .z = k * a};
                for (int l = 0; l < M; ++l) {
                    r[cont].x = R.x + b[l][0] * a;
                    r[cont].y = R.y + b[l][1] * a;
                    r[cont].z = R.z + b[l][2] * a;
                    cont++;
                }
            }
        }
    }

}

double V_LJ(double r, double L) {
    if (r < L / 2 && r != 0) {
        double VL_2 = 4 * (pow(2 / L, 12) - pow(2 / L, 6));
        //double VpL_2 = 24 * (pow(1 / r, 8) - 2 * pow(1 / r, 14));
        return 4 * (pow(1 / r, 12) - pow(1 / r, 6)) - VL_2;
    } else {
        return 0;
    }
}
void a_LJ(vec *r, vec *a, double r_c, double L) {
    vec *dr[N - 1]; vec_2D(dr, N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j != i) {
                dr[i][j].x = r[i].x - r[j].x;
                dr[i][j].y = r[i].y - r[j].y;
                dr[i][j].z = r[i].z - r[j].z;

                dr[i][j].x -= L * rint(dr[i][j].x / L);
                dr[i][j].y -= L * rint(dr[i][j].y / L);
                dr[i][j].z -= L * rint(dr[i][j].z / L);
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        a[i].uguale(0);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j != i) {
                // cout<<"i = "<<i<<"\n"<<endl;
                // stampa_coord(dr[i]);
                double r_mod = dr[i][j].mod();
                if (r_mod < r_c) {
                    double cost = 24 * (pow(1 / r_mod, 8) - 2 * pow(1 / r_mod, 14));
                    a[i].x += dr[i][j].x * cost;
                    a[i].y += dr[i][j].y * cost;
                    a[i].z += dr[i][j].z * cost;
                }
            }
        }
    }
}
void vel_verlet(vec *r, vec *v, vec *a, double dt, double r_c, double L, double *K, double *V) {
    vec a_prev[N];
    for (int i = 0; i < N; ++i) {
        a_prev[i].x = a[i].x;
        a_prev[i].y = a[i].y;
        a_prev[i].z = a[i].z;

        r[i].x += v[i].x * dt + 0.5 * dt * dt * a[i].x;
        r[i].y += v[i].y * dt + 0.5 * dt * dt * a[i].y;
        r[i].z += v[i].z * dt + 0.5 * dt * dt * a[i].z;
    }
    a_LJ(r, a, r_c, L);
    for (int i = 0; i < N; ++i) {
        v[i].x += 0.5 * dt * (a_prev[i].x + a[i].x);
        v[i].y += 0.5 * dt * (a_prev[i].y + a[i].y);
        v[i].z += 0.5 * dt * (a_prev[i].z + a[i].z);
    }

    *K = 0; *V = 0;
    double v_mod, r_mod;
    for (int i = 0; i < N; ++i) {
        v_mod = v[i].mod();
        r_mod = r[i].mod();

        *K += 0.5 * v_mod * v_mod;
        *V += V_LJ(r_mod, L);
    }
}