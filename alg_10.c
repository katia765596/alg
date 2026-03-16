#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
int gauss(int n, const double* A, const double* b, double* x) {
    double* a = (double*)malloc(n * n * sizeof(double));
    double* c = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n * n; ++i) a[i] = A[i];
    for (int i = 0; i < n; ++i) c[i] = b[i];
    for (int k = 0; k < n; ++k) {
        int maxrow = k;
        double maxval = fabs(a[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(a[i * n + k]) > maxval) {
                maxval = fabs(a[i * n + k]);
                maxrow = i;
            }
        }
        if (maxval < 1e-12) {
            free(a); free(c);
            return 1;   
        }
        if (maxrow != k) {
            for (int j = k; j < n; ++j) {
                double tmp = a[k * n + j];
                a[k * n + j] = a[maxrow * n + j];
                a[maxrow * n + j] = tmp;
            }
            double tmp = c[k];
            c[k] = c[maxrow];
            c[maxrow] = tmp;
        }
        double pivot = a[k * n + k];
        for (int j = k; j < n; ++j) a[k * n + j] /= pivot;
        c[k] /= pivot;
        for (int i = k + 1; i < n; ++i) {
            double factor = a[i * n + k];
            for (int j = k; j < n; ++j) {
                a[i * n + j] -= factor * a[k * n + j];
            }
            c[i] -= factor * c[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        x[i] = c[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= a[i * n + j] * x[j];
        }
    }
    free(a); free(c);
    return 0;
}
int newton_nd(
    int n,
    void (*F)(const double*, double*),
    void (*J)(const double*, double*),
    double* x,
    double eps,        
    int max_iter,
    int verbose
) {
    double* f = (double*)malloc(n * sizeof(double));
    double* dx = (double*)malloc(n * sizeof(double));
    double* x_prev = (double*)malloc(n * sizeof(double));
    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) x_prev[i] = x[i];
        F(x, f);
        double* Jmat = (double*)malloc(n * n * sizeof(double));
        J(x, Jmat);
        double* b = (double*)malloc(n * sizeof(double));
        for (int i = 0; i < n; ++i) b[i] = -f[i];
        int singular = gauss(n, Jmat, b, dx);
        free(Jmat); free(b);
        if (singular) {
            if (verbose) printf("   Якобиан вырожден на итерации %d\n", iter);
            free(f); free(dx); free(x_prev);
            return 1;
        }
        double normDx = 0.0;
        for (int i = 0; i < n; ++i) {
            x[i] += dx[i];
            double a = fabs(dx[i]);
            if (a > normDx) normDx = a;
        }
        if (verbose) {
            printf("   iter %2d: x = (", iter);
            for (int i = 0; i < n; ++i) printf("%12.8f ", x[i]);
            printf(")  ||dx||_inf = %8.2e\n", normDx);
        }
        if (normDx < eps) {
            if (verbose) printf("   Сошлось за %d итераций\n", iter + 1);
            free(f); free(dx); free(x_prev);
            return 0;
        }
    }
    if (verbose) printf("   Достигнуто максимальное число итераций\n");
    free(f); free(dx); free(x_prev);
    return 2;
}
void F10a(const double* x, double* f) {
    double x1 = x[0], x2 = x[1], x3 = x[2];
    f[0] = x1 * x1 * x1 + x1 * x1 * x2 - x1 * x3 + 6.0;
    f[1] = exp(x1) + exp(x2) - x3;
    f[2] = x2 * x2 - 2.0 * x1 * x3 - 4.0;
}
void J10a(const double* x, double* J) {
    double x1 = x[0], x2 = x[1], x3 = x[2];
    J[0] = 3.0 * x1 * x1 + 2.0 * x1 * x2 - x3;
    J[1] = x1 * x1;
    J[2] = -x1;
    J[3] = exp(x1);
    J[4] = exp(x2);
    J[5] = -1.0;
    J[6] = -2.0 * x3;
    J[7] = 2.0 * x2;
    J[8] = -2.0 * x1;
}
void F10b(const double* x, double* f) {
    double x1 = x[0], x2 = x[1], x3 = x[2];
    f[0] = 6.0 * x1 - 2.0 * cos(x2 * x3) - 1.0;
    double tmp = x1 * x1 + sin(x3) + 1.06;
    f[1] = 9.0 * x2 + sqrt(tmp) + 0.9; 
    f[2] = 60.0 * x3 + 3.0 * exp(-x1 * x2) + 10.0 * M_PI - 3.0;
}
void J10b(const double* x, double* J) {
    double x1 = x[0], x2 = x[1], x3 = x[2];
    double u = x2 * x3;
    double su = sin(u), cu = cos(u);
    J[0] = 6.0;
    J[1] = 2.0 * su * x3;
    J[2] = 2.0 * su * x2;
    double tmp = x1 * x1 + sin(x3) + 1.06;
    double sqrt_tmp = sqrt(tmp);
    if (sqrt_tmp < 1e-12) sqrt_tmp = 1e-12;
    J[3] = x1 / sqrt_tmp;
    J[4] = 9.0;
    J[5] = cos(x3) / (2.0 * sqrt_tmp);
    double exp_term = exp(-x1 * x2);
    J[6] = -3.0 * x2 * exp_term;
    J[7] = -3.0 * x1 * exp_term;
    J[8] = 60.0;
}
void F10c(const double* x, double* f) {
    double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3];
    f[0] = 4.0 * x1 - x2 + x3 - x1 * x4;
    f[1] = -x1 + 3.0 * x2 - 2.0 * x3 - x2 * x4;
    f[2] = x1 - 2.0 * x2 + 3.0 * x3 - x3 * x4;
    f[3] = x1 * x1 + x2 * x2 + x3 * x3 - 1.0;
}
void J10c(const double* x, double* J) {
    double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3];
    J[0] = 4.0 - x4;
    J[1] = -1.0;
    J[2] = 1.0;
    J[3] = -x1;
    J[4] = -1.0;
    J[5] = 3.0 - x4;
    J[6] = -2.0;
    J[7] = -x2;
    J[8] = 1.0;
    J[9] = -2.0;
    J[10] = 3.0 - x4;
    J[11] = -x3;
    J[12] = 2.0 * x1;
    J[13] = 2.0 * x2;
    J[14] = 2.0 * x3;
    J[15] = 0.0;
}
int main() {
    const double eps = 1e-8;    
    const int max_iter = 100;
    printf("========== Задание 10a ==========\n");
    double starts10a[][3] = {
        {-2.0, 0.0, 1.0},
        {-1.5, 0.5, 0.5},
        {-1.0,-1.0, 2.0},
        { 0.0, 0.0, 0.0}
    };
    int nst10a = sizeof(starts10a) / sizeof(starts10a[0]);
    for (int s = 0; s < nst10a; ++s) {
        double x[3] = { starts10a[s][0], starts10a[s][1], starts10a[s][2] };
        printf("Начальное приближение: (%5.2f, %5.2f, %5.2f)\n", x[0], x[1], x[2]);
        int res = newton_nd(3, F10a, J10a, x, eps, max_iter, 1);
        if (res == 0) {
            double f[3];
            F10a(x, f);
            printf("   Решение: (%12.8f, %12.8f, %12.8f)  F = (%8.2e, %8.2e, %8.2e)\n",
                x[0], x[1], x[2], f[0], f[1], f[2]);
        }
        else {
            printf("   Решение не найдено.\n");
        }
    }
    printf("\n========== Задание 10b ==========\n");
    double starts10b[][3] = {
        {0.0, -0.5, 0.0},
        {0.1, -0.4, 0.05},
        {0.2, -0.3, 0.1},
        {-0.1,-0.5, 0.0}
    };
    int nst10b = sizeof(starts10b) / sizeof(starts10b[0]);
    for (int s = 0; s < nst10b; ++s) {
        double x[3] = { starts10b[s][0], starts10b[s][1], starts10b[s][2] };
        printf("Начальное приближение: (%5.2f, %5.2f, %5.2f)\n", x[0], x[1], x[2]);
        int res = newton_nd(3, F10b, J10b, x, eps, max_iter, 1);
        if (res == 0) {
            double f[3];
            F10b(x, f);
            printf("   Решение: (%12.8f, %12.8f, %12.8f)  F = (%8.2e, %8.2e, %8.2e)\n",
                x[0], x[1], x[2], f[0], f[1], f[2]);
        }
        else {
            printf("   Решение не найдено.\n");
        }
    }
    printf("\n========== Задание 10c ==========\n");
    double starts10c[][4] = {
        {0.5, 0.5, 0.5, 1.0},
        {0.3, 0.4, 0.5, 2.0},
        {1.0, 0.0, 0.0, 3.0},
        {0.0, 1.0, 0.0, 2.0},
        {0.0, 0.0, 1.0, 1.0}
    };
    int nst10c = sizeof(starts10c) / sizeof(starts10c[0]);
    for (int s = 0; s < nst10c; ++s) {
        double x[4] = { starts10c[s][0], starts10c[s][1], starts10c[s][2], starts10c[s][3] };
        printf("Начальное приближение: (%5.2f, %5.2f, %5.2f, %5.2f)\n",
            x[0], x[1], x[2], x[3]);
        int res = newton_nd(4, F10c, J10c, x, eps, max_iter, 1);
        if (res == 0) {
            double f[4];
            F10c(x, f);
            printf("   Решение: (%12.8f, %12.8f, %12.8f, %12.8f)  F = (%8.2e, %8.2e, %8.2e, %8.2e)\n",
                x[0], x[1], x[2], x[3], f[0], f[1], f[2], f[3]);
        }
        else {
            printf("   Решение не найдено.\n");
        }
    }
    return 0;
}