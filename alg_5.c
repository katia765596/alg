#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 1e-6
#define MAX_ITER 50
typedef double (*PhiFunc)(double, double, double, double);
int simple_iteration(PhiFunc phi, double x0, double p1, double p2, double p3, double* root, int* iterations)
{
    double x_prev = x0;
    double x_next;
    for (int k = 1; k <= MAX_ITER; k++)
    {
        x_next = phi(x_prev, p1, p2, p3);
        if (fabs(x_next - x_prev) < EPS)
        {
            *root = x_next;
            *iterations = k;
            return 1; 
        }
        x_prev = x_next;
    }
    *root = x_prev;
    *iterations = MAX_ITER;
    return 0; 
}
double phi_a(double x, double p1, double p2, double p3) { (void)p1; (void)p2; (void)p3; return 2 * x - 1; }
double phi_b(double x, double p1, double p2, double p3) { (void)p1; (void)p2; (void)p3; return exp(2 * x) - 1; }
double phi_c(double x, double A, double p2, double p3) { (void)p2; (void)p3; return x <= 0 ? 1e6 : A - log(x); }
double phi_d(double x, double alpha, double beta, double p3) { (void)p3; return alpha * exp(-x) + beta * x; }
void investigate(PhiFunc phi, double p1, double p2, double p3, double x_start, double x_end, double step)
{
    printf("\nТаблица сходимости:\n");
    printf("x0      | Результат | Кол-во итераций | Сходимость\n");
    printf("-----------------------------------------------------\n");
    int num_steps = (int)((x_end - x_start) / step + 0.5) + 1;
    if (num_steps <= 0) num_steps = 1;
    char* graph = malloc((num_steps + 1) * sizeof(char));
    if (!graph) { printf("Ошибка выделения памяти!\n"); return; }
    graph[num_steps] = '\0';
    int idx = 0;
    for (int i = 0; i < num_steps; i++)
    {
        double x0 = x_start + i * step;
        double root;
        int iter;
        int converged = simple_iteration(phi, x0, p1, p2, p3, &root, &iter);
        if (converged)
        {
            printf("%7.3f | %10.6f | %15d |   + (сошлось)\n", x0, root, iter);
            graph[idx++] = '+';
        }
        else
        {
            printf("%7.3f | %10s | %15d |   - (не сошлось)\n", x0, "N/A", iter);
            graph[idx++] = '-';
        }
    }
    graph[idx] = '\0';
    printf("\nГоризонтальная визуализация сходимости по x0:\n");
    printf("%s\n", graph);
    printf("( '+' = сошлось, '-' = не сошлось )\n");
    free(graph);
}
int main()
{
    int choice;
    double p1 = 0, p2 = 0, p3 = 0;
    double x_start, x_end, step;
    printf("Исследование сходимости метода простой итерации\n");
    printf("1) x_{n+1} = 2*x_n - 1\n");
    printf("2) x_{n+1} = e^{2*x_n} - 1\n");
    printf("3) x_{n+1} = A - ln(x_n)\n");
    printf("4) x_{n+1} = alpha * e^{-x_n} + beta * x_n\n");
    printf("Выберите вариант (1-4): ");
    scanf("%d", &choice);
    switch (choice)
    {
    case 3: printf("Введите A: "); scanf("%lf", &p1); break;
    case 4:
        printf("Введите alpha: "); scanf("%lf", &p1);
        printf("Введите beta: "); scanf("%lf", &p2);
        break;
    }
    printf("Введите начало интервала x0: "); scanf("%lf", &x_start);
    printf("Введите конец интервала x0: "); scanf("%lf", &x_end);
    printf("Введите шаг перебора x0: "); scanf("%lf", &step);
    switch (choice)
    {
    case 1: investigate(phi_a, 0, 0, 0, x_start, x_end, step); break;
    case 2: investigate(phi_b, 0, 0, 0, x_start, x_end, step); break;
    case 3: investigate(phi_c, p1, 0, 0, x_start, x_end, step); break;
    case 4: investigate(phi_d, p1, p2, 0, x_start, x_end, step); break;
    default: printf("Неверный выбор\n"); return 0;
    }
    return 0;
}