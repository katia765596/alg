#include <stdio.h>
#include <math.h>
#include <stdlib.h>
int gauss_2x2(const double* A, const double* b, double* x) {
    double det = A[0] * A[3] - A[1] * A[2];
    if (fabs(det) < 1e-12) return 1;   
    x[0] = (b[0] * A[3] - b[1] * A[1]) / det;
    x[1] = (A[0] * b[1] - A[2] * b[0]) / det;
    return 0;
}
int newton_2d(
    void (*F)(const double*, double*),        
    void (*J)(const double*, double*),      
    double* x,                              
    double eps,                           
    int max_iter,                           
    int verbose                               
) {
    double f[2], x_prev[2], dx[2];
    for (int iter = 0; iter < max_iter; ++iter) {
        x_prev[0] = x[0]; x_prev[1] = x[1];
        F(x, f);
        double normF = hypot(f[0], f[1]);
        double Jmat[4];
        J(x, Jmat);
        double b[2] = { -f[0], -f[1] };
        if (gauss_2x2(Jmat, b, dx)) {
            if (verbose) printf("   Якобиан вырожден на итерации %d\n", iter);
            return 1;
        }
        x[0] += dx[0];
        x[1] += dx[1];
        double normDx = hypot(dx[0], dx[1]);
        if (verbose) {
            printf("   iter %2d: x = (%12.8f, %12.8f)  ||F||=%8.2e  ||dx||=%8.2e\n",
                iter, x[0], x[1], normF, normDx);
        }
        if (normF < eps && normDx < eps) {
            if (verbose) printf("   Сошлось за %d итераций\n", iter + 1);
            return 0;
        }
    }
    if (verbose) printf("   Достигнуто максимальное число итераций\n");
    return 2;
}
double A_param, a_param, b_param;   
void F9a(const double* x, double* f) {
    double X = x[0], Y = x[1];
    f[0] = tan(X * Y + A_param) - X * X;
    f[1] = a_param * X * X + b_param * Y * Y - 1.0;
}
void J9a(const double* x, double* J) {
    double X = x[0], Y = x[1];
    double u = X * Y + A_param;
    double sec2 = 1.0 / (cos(u) * cos(u));
    J[0] = Y * sec2 - 2.0 * X;      
    J[1] = X * sec2;               
    J[2] = 2.0 * a_param * X;         
    J[3] = 2.0 * b_param * Y;      
}
void F9b(const double* x, double* f) {
    double x1 = x[0], x2 = x[1];
    f[0] = x1 * x1 + x2 * x2 - 2.0;
    f[1] = exp(x1 - 1.0) + x2 * x2 * x2 - 2.0;
}
void J9b(const double* x, double* J) {
    double x1 = x[0], x2 = x[1];
    J[0] = 2.0 * x1;
    J[1] = 2.0 * x2;
    J[2] = exp(x1 - 1.0);
    J[3] = 3.0 * x2 * x2;
}
int main() {
    const double eps = 1e-6;
    const int max_iter = 100;
    printf("========== Задание 9a ==========\n");
    struct {
        double A, a, b;
        const char* desc;
    } params9a[] = {
        {0.2, 0.6, 2.0, "A=0.2, a=0.6, b=2.0"},
        {0.4, 0.8, 2.0, "A=0.4, a=0.8, b=2.0"},
        {0.3, 0.2, 3.0, "A=0.3, a=0.2, b=3.0"},
        {0.0, 0.6, 2.0, "A=0.0, a=0.6, b=2.0"}
    };
    double starts9a[][2] = {
        { 0.5,  0.5},
        { 0.5, -0.5},
        {-0.5,  0.5},
        {-0.5, -0.5},
        { 0.2,  0.2},
        { 0.8,  0.3}
    };
    int nst9a = sizeof(starts9a) / sizeof(starts9a[0]);
    for (int p = 0; p < 4; ++p) {
        A_param = params9a[p].A;
        a_param = params9a[p].a;
        b_param = params9a[p].b;
        printf("\n--- %s ---\n", params9a[p].desc);
        for (int s = 0; s < nst9a; ++s) {
            double x[2] = { starts9a[s][0], starts9a[s][1] };
            printf("Начальное приближение: (%5.2f, %5.2f)\n", x[0], x[1]);
            int res = newton_2d(F9a, J9a, x, eps, max_iter, 1);
            if (res == 0) {
                double f[2];
                F9a(x, f);
                printf("   Решение: (%12.8f, %12.8f)  F = (%8.2e, %8.2e)\n",
                    x[0], x[1], f[0], f[1]);
            }
            else {
                printf("   Решение не найдено.\n");
            }
        }
    }
    printf("\n========== Задание 9b ==========\n");
    double starts9b[][2] = {
        {0.5, 0.5},
        {1.2, 0.8},
        {0.8, 1.2},
        {1.0, 1.0}
    };
    int nst9b = sizeof(starts9b) / sizeof(starts9b[0]);
    for (int s = 0; s < nst9b; ++s) {
        double x[2] = { starts9b[s][0], starts9b[s][1] };
        printf("Начальное приближение: (%5.2f, %5.2f)\n", x[0], x[1]);
        int res = newton_2d(F9b, J9b, x, eps, max_iter, 1);
        if (res == 0) {
            double f[2];
            F9b(x, f);
            printf("   Решение: (%12.8f, %12.8f)  F = (%8.2e, %8.2e)\n",
                x[0], x[1], f[0], f[1]);
        }
        else {
            printf("   Решение не найдено.\n");
        }
    }
    return 0;
}