#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void multiply(const double* a, int na, const double* b, int nb, double* res) {
    for (int i = 0; i < na + nb - 1; i++)
        res[i] = 0.0;
    for (int i = 0; i < na; i++)
        for (int j = 0; j < nb; j++)
            res[i + j] += a[i] * b[j];
}
double* invert_direct(const double* f, int n_f, int n) {
    if (n <= 0) return NULL;
    double* g = (double*)calloc(n, sizeof(double));
    if (f[0] == 0.0) {
        fprintf(stderr, "Ошибка: свободный член равен нулю, ряд не обратим.\n");
        free(g);
        return NULL;
    }
    g[0] = 1.0 / f[0];
    for (int k = 1; k < n; k++) {
        double sum = 0.0;
        for (int i = 1; i <= k; i++) {
            if (i < n_f) 
                sum += f[i] * g[k - i];
        }
        g[k] = -sum / f[0];
    }
    return g;
}
double* invert_newton(const double* f, int n_f, int n) {
    if (n <= 0) return NULL;
    if (f[0] == 0.0) {
        fprintf(stderr, "Ошибка: свободный член равен нулю.\n");
        return NULL;
    }
    double* g = (double*)malloc(sizeof(double));
    g[0] = 1.0 / f[0];
    int cur_len = 1; 
    while (cur_len < n) {
        int new_len = (cur_len * 2 < n) ? cur_len * 2 : n; 
        int t_len = new_len; 
        double* t = (double*)calloc(t_len, sizeof(double));
        for (int i = 0; i < t_len; i++) {
            double s = 0.0;
            for (int j = 0; j <= i; j++) {
                int fi = j;
                int gi = i - j;
                if (fi < n_f && gi < cur_len)
                    s += f[fi] * g[gi];
            }
            t[i] = s;
        }
        double* u = (double*)malloc(new_len * sizeof(double));
        for (int i = 0; i < new_len; i++) {
            u[i] = (i == 0) ? (2.0 - t[0]) : (-t[i]);
        }
        double* new_g = (double*)calloc(new_len, sizeof(double));
        for (int i = 0; i < new_len; i++) {
            double s = 0.0;
            for (int j = 0; j <= i; j++) {
                int gj = j;
                int ui = i - j;
                if (gj < cur_len && ui < new_len)
                    s += g[gj] * u[ui];
            }
            new_g[i] = s;
        }
        free(g);
        free(t);
        free(u);
        g = new_g;
        cur_len = new_len;
    }
    return g;
}
void print_series(const char* name, const double* coeff, int n) {
    printf("%s = ", name);
    for (int i = 0; i < n; i++) {
        printf("%+.6f x^%d ", coeff[i], i);
        if (i < n - 1) printf("+ ");
    }
    printf("\n");
}
int main() {
    int M = 10;
    double* f1 = (double*)malloc((M + 1) * sizeof(double));
    for (int k = 0; k <= M; k++) {
        double fact = 1.0;
        for (int i = 1; i <= k; i++) fact *= i;
        f1[k] = 1.0 / fact;
    }
    printf("Пример 1: f(x) = sum_{k=0}^{%d} x^k/k!\n", M);
    print_series("f", f1, M + 1);
    int n = 8;
    double* g_direct = invert_direct(f1, M + 1, n);
    double* g_newton = invert_newton(f1, M + 1, n);
    if (g_direct) {
        printf("\nПрямой метод, первые %d коэффициентов g:\n", n);
        print_series("g", g_direct, n);
        free(g_direct);
    }
    if (g_newton) {
        printf("Метод Ньютона, первые %d коэффициентов g:\n", n);
        print_series("g", g_newton, n);
        free(g_newton);
    }
    double f2[] = { -1.0, -1.0, 1.0 }; 
    int n_f2 = 3;
    printf("\nПример 2: f(x) = x^2 - x - 1\n");
    print_series("f", f2, n_f2);

    n = 10;
    g_direct = invert_direct(f2, n_f2, n);
    g_newton = invert_newton(f2, n_f2, n);

    if (g_direct) {
        printf("\nПрямой метод, первые %d коэффициентов g:\n", n);
        print_series("g", g_direct, n);
        free(g_direct);
    }
    if (g_newton) {
        printf("Метод Ньютона, первые %d коэффициентов g:\n", n);
        print_series("g", g_newton, n);
        free(g_newton);
    }
    free(f1);
    return 0;
}
