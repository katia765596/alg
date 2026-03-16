#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#define EPS 1e-6
#define MAX_ITER 100
#define MAX_ROOTS 20
typedef double (*PhiFunction)(double);
typedef struct
{
    char name[100];
    PhiFunction phi;
} Equation;
double simple_iteration(PhiFunction phi, double x0, int* iterations)
{
    double x_prev = x0;
    double x_next;
    printf("\n---------------------------------------------\n");
    printf(" k        x_k           |x_k - x_k-1|\n");
    printf("---------------------------------------------\n");
    for (int k = 1; k <= MAX_ITER; k++)
    {
        x_next = phi(x_prev);
        printf("%2d   %12.8f    %12.8f\n",
            k, x_next, fabs(x_next - x_prev));
        if (fabs(x_next - x_prev) < EPS)
        {
            *iterations = k;
            return x_next;
        }
        x_prev = x_next;
    }
    *iterations = MAX_ITER;
    return x_prev;
}
int is_new_root(double root, double roots[], int count)
{
    for (int i = 0; i < count; i++)
        if (fabs(root - roots[i]) < EPS)
            return 0;
    return 1;
}
double phi_a(double x)
{
    return (1 - x * x * x) / (3 * x);
}
double phi_b(double x)
{
    return pow(x, 3);
}
double phi_c(double x)
{
    return (x * x + 2) / 3;
}
void find_roots(Equation eq, double a, double b, double step)
{
    double roots[MAX_ROOTS];
    int root_count = 0;
    printf("\n=================================\n");
    printf("Уравнение: %s\n", eq.name);
    printf("Интервал поиска: [%.2f , %.2f]\n", a, b);
    for (double x = a; x <= b; x += step)
    {
        int iterations;
        double root = simple_iteration(eq.phi, x, &iterations);
        if (is_new_root(root, roots, root_count))
        {
            roots[root_count++] = root;
            printf("\n>>> Найден корень: %.8f\n", root);
            printf(">>> Итераций: %d\n", iterations);
        }
    }
    printf("\nВсего найдено корней: %d\n", root_count);

    for (int i = 0; i < root_count; i++)
        printf("Корень %d = %.8f\n", i + 1, roots[i]);
}
int main()
{
    Equation eq1 = { "x^3 + 3x^2 - 1 = 0", phi_a };
    Equation eq2 = { "x^4 - x^3 = 0", phi_b };
    Equation eq3 = { "x^2 - 3x + 2 = 0", phi_c };
    int choice;
    printf("Метод простой итерации\n");
    printf("1 - x^3 + 3x^2 - 1 = 0\n");
    printf("2 - x^4 - x^3 = 0\n");
    printf("3 - x^2 - 3x + 2 = 0\n");
    scanf("%d", &choice);
    switch (choice)
    {
    case 1:
        find_roots(eq1, -3, 2, 0.5);
        break;

    case 2:
        find_roots(eq2, -2, 2, 0.5);
        break;

    case 3:
        find_roots(eq3, -1, 3, 0.5);
        break;

    default:
        printf("Неверный выбор\n");
    }
    return 0;
}