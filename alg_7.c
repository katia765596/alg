#include <stdio.h>
#include <math.h>
#define EPS 1e-12
#define MAX_ITER 1000
#define STEP_INIT 0.5
#define STEP_MIN 0.01
#define STEP_MAX 0.5
#define ROOTS_NEEDED 3
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
double f(double x) {
    return pow(x - 1, 3) * sin(M_PI * x) * (cos(2 * M_PI * x) - 1);
}
double df(double x) {
    double term1 = 3 * pow(x - 1, 2) * sin(M_PI * x) * (cos(2 * M_PI * x) - 1);
    double term2 = pow(x - 1, 3) * M_PI * cos(M_PI * x) * (cos(2 * M_PI * x) - 1);
    double term3 = pow(x - 1, 3) * sin(M_PI * x) * (-2 * M_PI * sin(2 * M_PI * x));
    return term1 + term2 + term3;
}
double d2f(double x) {
    double t1 = 6 * (x - 1) * sin(M_PI * x) * (cos(2 * M_PI * x) - 1);
    double t2 = 6 * pow(x - 1, 2) * M_PI * cos(M_PI * x) * (cos(2 * M_PI * x) - 1);
    double t3 = 3 * pow(x - 1, 2) * sin(M_PI * x) * (-2 * M_PI * sin(2 * M_PI * x));
    double t4 = pow(x - 1, 3) * (-M_PI * M_PI * sin(M_PI * x)) * (cos(2 * M_PI * x) - 1);
    double t5 = pow(x - 1, 3) * M_PI * cos(M_PI * x) * (-2 * M_PI * sin(2 * M_PI * x));
    double t6 = pow(x - 1, 3) * sin(M_PI * x) * (-4 * M_PI * M_PI * cos(2 * M_PI * x));
    return t1 + t2 + t3 + t4 + t5 + t6;
}
int already_found(double roots[], int count, double candidate) {
    for (int i = 0; i < count; i++)
        if (fabs(candidate - roots[i]) < 1e-6) return 1;
    return 0;
}
double newton_auto(double (*f)(double), double (*df)(double), double (*d2f)(double), double x0) {
    double x = x0;
    for (int i = 0; i < MAX_ITER; i++) {
        double fx = f(x);
        double dfx = df(x);
        double d2fx = d2f(x);
        if (fabs(dfx) < 1e-14) break;
        double p = 1.0;
        if (fabs(fx) > EPS) {
            p = fabs(fx * d2fx / (dfx * dfx)) + 1.0;
        }
        double x_next = x - p * fx / dfx;
        if (fabs(x_next - x) < EPS) return x_next;
        x = x_next;
    }
    return x;
}
int main() {
    double roots[ROOTS_NEEDED];
    int count = 0;
    double x = 0.5;
    double step = STEP_INIT;
    while (count < ROOTS_NEEDED && x < 20) {
        double root = newton_auto(f, df, d2f, x);
        if (root > 0 && !already_found(roots, count, root)) {
            roots[count++] = root;
        }
        double fx = fabs(f(x));
        if (fx > 1.0) step = (step / 2 > STEP_MIN) ? step / 2 : STEP_MIN;
        else step = (step * 1.2 < STEP_MAX) ? step * 1.2 : STEP_MAX;
        x += step;
    }
    for (int i = 0; i < count; i++) {
        printf("Root %d: %.12f\n", i + 1, roots[i]);
    }
    return 0;
}