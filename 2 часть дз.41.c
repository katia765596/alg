#include <stdio.h>

void compl_(double a, double b, double c, double d, double* re, double* lm) {
    double p = a * c;       
    double q = b * d;          
    double r = (a + b) * (c + d); 
    *re = p - q;             
    *lm = r - p - q;        
}
int main() {
    double a = 3.0, b = 2.0, c = 1.0, d = 4.0;
    double re, lm;
    compl_(a, b, c, d, &re, &lm);
    printf("(%g + %gi) * (%g + %gi) = %g + %gi\n", a, b, c, d, re, lm);
    return 0;
}