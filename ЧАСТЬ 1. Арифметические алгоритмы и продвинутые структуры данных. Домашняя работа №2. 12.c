#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
/*ЧАСТЬ 1 1),2),3),4),5) ЗАДАНИЕ*/
/* полином */
void polynom(uint64_t a) {
    int flag=1;
    for (int i=63;i>=0;i--) {
        if ((a>>i)&1) {
            if (flag==0) 
                printf(" + ");
            flag=0;
            if (i==0) 
                printf("1");
            else if (i==1) 
                printf("x");
            else 
                printf("x^%d", i);
        }
    }
    if (flag==1) 
        printf("0");
    printf("\n");
}
/* бит запись */
uint64_t bit_polynom(const char* s) {
    uint64_t res=0;
    int i=0;
    if (s==0) 
        return 0;
    while (s[i]) {
        if (s[i]==' ') 
        { i++; continue; }
        if (s[i]=='1') 
        { res|=1;i++; }
        else if (s[i]=='x') {
            i++;
            if (s[i]=='^') {
                i++;
                int p = 0;
                while (isdigit(s[i])) 
                { p=p*10+(s[i]-'0'); i++; }
                res |= (1ULL << p);
            }
            else 
            { res |= (1ULL << 1); }
        }
        else { i++; }
    }
    return res;
}
/* степени */
int degree(uint64_t a) {
    for (int i=63; i>=0; i--)
        if ((a>>i)&1) 
            return i;
    return -1;
}
/* умножение 32 алг умн в столбик */
uint64_t polynom_ymn(uint64_t a, uint64_t b) {
    uint64_t res=0;
    int db = degree(b);
    for (int i=0; i<=db; i++)
        if ((b>>i)&1)
            res ^= (a<<i);
    return res;
}
/* деление для расширенного алгоритма евклидаа */
void polynom_div(uint64_t a, uint64_t b, uint64_t* ch, uint64_t* o) {
    *ch = 0;
    *o = a;
    int db = degree(b);
    while (degree(*o)>=db) {
        int difff=degree(*o)-db;
        *ch^=(1ULL<<difff);
        *o^=(b<<difff);
    }
}
/* умножение c непр полиномом */
uint64_t nepr_mul(uint64_t a, uint64_t b, uint64_t f) {
    uint64_t res = polynom_ymn(a, b);
    int n = degree(f);
    for (int i=degree(res);i>=n;i--)
        if ((res>>i)&1)
            res^=(f<<(i - n));
    return res;
}
/* расширенный алгоритм евклида, возвр обр эл */
uint64_t evkl(uint64_t a, uint64_t f) {
    uint64_t dol=f, del=a;
    uint64_t u=0, v=1;
    while (del!=0) {
        uint64_t ch, o;
        polynom_div(dol, del, &ch, &o);
        uint64_t t =u^polynom_ymn(ch, v);
        dol = del; del = o; u = v; v = t;
    }
    return u;
}
int main() {

        char input[500], str_a[100], str_b[100];

        printf("Enter irreducible polynomials f(x) separated by commas:\n");
        fgets(input, 500, stdin);
        input[strcspn(input, "\n")] = 0;
        printf("Enter element a: ");
        fgets(str_a, 100, stdin);
        str_a[strcspn(str_a, "\n")] = 0; 
        printf("Enter element b: ");
        fgets(str_b, 100, stdin);
        str_b[strcspn(str_b, "\n")] = 0;
        uint64_t a = bit_polynom(str_a);
        uint64_t b = bit_polynom(str_b);
    char* token = strtok(input, ",");
    int count = 1;
    while (token) {
        uint64_t f = bit_polynom(token);
        if (degree(f) <= 0) {
            printf("\nPolynomial %s is invalid (zero or constant). Skipping.\n", token);
            token = strtok(NULL, ",");
            count++;
            continue;
        }
        printf("\n==============================\n");
        printf("Polynomial #%d: ", count++);
        polynom(f);
        printf("Element a: ");
        polynom(a);
        printf("Element b: ");
        polynom(b);
        printf("\nMultiplication a*b (without reduction):\n");
        polynom(polynom_ymn(a, b));
        printf("\nMultiplication in GF(2^n):\n");
        polynom(nepr_mul(a, b, f));
        printf("\nInverse of a in GF(2^n):\n");
        uint64_t inv = evkl(a, f);
        polynom(inv);
        printf("\nCheck: a * a^(-1) mod f = ");
        polynom(nepr_mul(a, inv, f));
        token = strtok(NULL, ",");
    }
    return 0;
}