#include <stdio.h>
#define _CRT_SECURE_NO_WARNINGS
void com_div_min(
    double a0, double a1,
    double b0, double b1,
    double* re, double* im)
{
    // карацуба
    double y1 = a0 * b0;
    double y2 = a1 * b1;
    double y3 = (a0 + a1) * (b0 - b1);
    double re_num = y1 + y2;
    double mn_num = y3 - y1 + y2;
    double znm = b0 * b0 + b1 * b1;
    *re = re_num / znm;
    *im = mn_num / znm;
}
//примеееер
int main()
{
    double re, im;
    com_div_min(4, 3, 1, -2, &re, &im);
    printf("res: %.6f + %.6fi\n", re, im);
    return 0;
}