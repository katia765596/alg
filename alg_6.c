#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_ITER 1000
typedef double (*Func)(double);         
typedef double (*FuncDeriv)(double);    
int bisection(Func f, double a, double b, double eps, double* root, int* iterations)
{
    double fa = f(a), fb = f(b);
    if (fa * fb > 0) return 0; 
    int k = 0;
    while ((b - a) / 2 > eps && k < MAX_ITER)
    {
        double c = (a + b) / 2;
        double fc = f(c);
        k++;
        printf("  [Бисекция] итерация %d: x = %.12f\n", k, c);
        if (fabs(fc) < eps) { *root = c; *iterations = k; return 1; }
        if (fa * fc < 0) { b = c; fb = fc; }
        else { a = c; fa = fc; }
    }
    *root = (a + b) / 2;
    *iterations = k;
    return 1;
}
int newton(Func f, FuncDeriv df, double x0, double eps, double* root, int* iterations)
{
    double x = x0;
    int k = 0;
    while (k < MAX_ITER)
    {
        double fx = f(x);
        double dfx = df(x);
        k++;
        printf("  [Ньютон] итерация %d: x = %.12f\n", k, x);
        if (fabs(fx) < eps) { *root = x; *iterations = k; return 1; }
        if (fabs(dfx) < 1e-12) return 0; 
        x = x - fx / dfx;
    }
    *root = x;
    *iterations = k;
    return 1;
}
double f_a(double x) { return sin(x) - 2 * x * x + 0.5; }
double df_a(double x) { return cos(x) - 4 * x; }
double f_b(double x) { return x * x * x - 8; }
double df_b(double x) { return 3 * x * x; }
double f_c(double x) { return x <= -1.0 || x >= 1.0 ? 1e6 : sqrt(1 - x * x) - exp(x) + 0.1; }
double df_c(double x) { return x <= -1.0 || x >= 1.0 ? 0 : (-x) / sqrt(1 - x * x) - exp(x); }
double f_d(double x) { return pow(x, 6) - 5 * pow(x, 3) - 2; }
double df_d(double x) { return 6 * pow(x, 5) - 15 * pow(x, 2); }
double f_e(double x) { return x <= 0 ? 1e6 : log2(x) - 1 / (1 + x * x); }
double df_e(double x) { return x <= 0 ? 0 : 1 / (x * log(2)) + 2 * x / ((1 + x * x) * (1 + x * x)); }
double f_f(double x) { return sin(x) * sin(x) - 1; }
double df_f(double x) { return 2 * sin(x) * cos(x); }
double f_g(double x) { return x <= 0 ? 1e6 : log(x) - 1; }
double df_g(double x) { return x <= 0 ? 0 : 1 / x; }
void find_roots(Func f, FuncDeriv df, double x_start, double x_end, double delta, double eps)
{
    double roots[3];
    int count = 0;
    printf("\nЛокализация корней и сравнение методов:\n");
    printf("--------------------------------------------------------\n");

    for (double x = x_start; x < x_end && count < 3; x += delta)
    {
        double a = x, b = x + delta;
        if (f(a) * f(b) <= 0) 
        {
            double r_b, r_n;
            int iter_b, iter_n;
            printf("\nКорень %d на [%lf,%lf]:\n", count + 1, a, b);
            int conv_b = bisection(f, a, b, eps, &r_b, &iter_b);
            int conv_n = newton(f, df, (a + b) / 2, eps, &r_n, &iter_n);
            if (conv_b && conv_n)
            {
                printf("  Результат:\n");
                printf("    Бисекция: %.12f (итераций %d)\n", r_b, iter_b);
                printf("    Ньютон:    %.12f (итераций %d)\n", r_n, iter_n);
            }
            else
            {
                printf("  Один из методов не сошелся.\n");
            }
            roots[count++] = r_b; 
        }
    }
    if (count == 0) printf("Корни не найдены на заданном интервале.\n");
}
int main()
{
    int choice, n;
    double eps;
    double x_start, x_end, delta;
    printf("Выберите уравнение (1-7):\n");
    printf("1) sin(x)-2x^2+0.5=0\n");
    printf("2) x^3-8=0\n");
    printf("3) sqrt(1-x^2)-exp(x)+0.1=0\n");
    printf("4) x^6-5x^3-2=0\n");
    printf("5) log2(x)-1/(1+x^2)=0\n");
    printf("6) sin(x)^2-1=0\n");
    printf("7) ln(x)-1=0\n");
    scanf("%d", &choice);
    printf("Введите n (точность eps = 10^-n): ");
    scanf("%d", &n);
    eps = pow(10, -n);
    printf("Введите начало интервала локализации x_start: ");
    scanf("%lf", &x_start);
    printf("Введите конец интервала x_end: ");
    scanf("%lf", &x_end);
    printf("Введите шаг локализации delta: ");
    scanf("%lf", &delta);
    switch (choice)
    {
    case 1: find_roots(f_a, df_a, x_start, x_end, delta, eps); break;
    case 2: find_roots(f_b, df_b, x_start, x_end, delta, eps); break;
    case 3: find_roots(f_c, df_c, x_start, x_end, delta, eps); break;
    case 4: find_roots(f_d, df_d, x_start, x_end, delta, eps); break;
    case 5: find_roots(f_e, df_e, x_start, x_end, delta, eps); break;
    case 6: find_roots(f_f, df_f, x_start, x_end, delta, eps); break;
    case 7: find_roots(f_g, df_g, x_start, x_end, delta, eps); break;
    default: printf("Неверный выбор\n"); return 0;
    }
    return 0;
}