#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846
#define MAX_ITER 1000
typedef double (*Func)(double);
typedef double (*Deriv)(double);
double f(double x)
{
    return pow(x - 1, 3) * sin(PI * x) * (cos(2 * PI * x) - 1);
}
double df(double x)
{
    double u = pow(x - 1, 3);
    double du = 3 * pow(x - 1, 2);
    double v = sin(PI * x);
    double dv = PI * cos(PI * x);
    double w = cos(2 * PI * x) - 1;
    double dw = -2 * PI * sin(2 * PI * x);
    return du * v * w + u * dv * w + u * v * dw;
}
int newton(Func f, Deriv df, double x0, double eps, double* root)
{
    double x = x0;
    for (int k = 1; k <= MAX_ITER; k++)
    {
        double d = df(x);
        if (fabs(d) < 1e-12)
            return -1;
        double x_next = x - f(x) / d;
        printf("Newton iter %d  x = %.10f\n", k, x_next);
        if (fabs(x_next - x) < eps)
        {
            *root = x_next;
            return k;
        }
        x = x_next;
    }
    return -1;
}
int modified_newton(Func f, Deriv df, double x0, double sigma, double eps, double* root)
{
    double x = x0;
    for (int k = 1; k <= MAX_ITER; k++)
    {
        double d = df(x);
        if (fabs(d) < 1e-12)
            return -1;
        double x_next = x + sigma * f(x) / d;
        printf("Modified iter %d  x = %.10f\n", k, x_next);
        if (fabs(x_next - x) < eps)
        {
            *root = x_next;
            return k;
        }
        x = x_next;
    }
    return -1;
}
int main()
{
    double roots[3] = { 1.0, 2.0, 3.0 }; 
    int p = 3;
    double sigma = -p;  
    int n_values[4] = { 3,4,5,6 };
    for (int k = 0; k < 4; k++)
    {
        double eps = pow(10, -n_values[k]);
        printf("\n==============================\n");
        printf("eps = 10^-%d\n", n_values[k]);
        printf("==============================\n");
        for (int i = 0; i < 3; i++)
        {
            double x0 = roots[i] + 0.2;
            printf("\nRoot near %.2f\n", roots[i]);
            double r1, r2;
            printf("\nClassic Newton\n");
            int it1 = newton(f, df, x0, eps, &r1);
            printf("\nModified Newton\n");
            int it2 = modified_newton(f, df, x0, sigma, eps, &r2);
            printf("\nResult:\n");
            printf("Newton root = %.10f  iterations = %d\n", r1, it1);
            printf("Modified root = %.10f  iterations = %d\n", r2, it2);
            printf("\n----------------------------\n");
        }
    }
}

