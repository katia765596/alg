#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#define STEP 0.01
typedef double (*Func)(double, double, double, double);
double phi_a(double x, double a, double b, double c)
{
    return c + a * sin(2 * x) + b * cos(2 * x);
}
double phi_b(double x, double a, double b, double c)
{
    return c + a * exp(-b * x * x);
}
double dphi_a(double x, double a, double b, double c)
{
    return 2 * a * cos(2 * x) - 2 * b * sin(2 * x);
}
double dphi_b(double x, double a, double b, double c)
{
    return -2 * a * b * x * exp(-b * x * x);
}
void research(Func dphi, double a, double b, double c, double left, double right)
{
    double max = 0;
    for (double x = left; x <= right; x += STEP)
    {
        double val = fabs(dphi(x, a, b, c));
        if (val > max)
            max = val;
    }
    printf("\nmax |phi'(x)| = %.6f\n", max);
    if (max < 1)
        printf("Метод простой итерации сходится при любом начальном приближении\n");
    else
        printf("Сходимость не гарантирована\n");
}
int main()
{
    double a, b, c;
    int choice;
    printf("Исследование сходимости метода простой итерации\n");
    printf("1) φ(x) = c + a*sin(2x) + b*cos(2x)\n");
    printf("2) φ(x) = c + a*exp(-b*x^2)\n");
    scanf("%d", &choice);
    printf("Введите a: ");
    scanf("%lf", &a);
    printf("Введите b: ");
    scanf("%lf", &b);
    printf("Введите c: ");
    scanf("%lf", &c);
    switch (choice)
    {
    case 1:
        research(dphi_a, a, b, c, -10, 10);
        break;

    case 2:
        research(dphi_b, a, b, c, -10, 10);
        break;

    default:
        printf("Неверный выбор\n");
    }
    return 0;
}