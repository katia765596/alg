#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_VARS 20
#define MAX_LEN 256
#define EPS 1e-9


#define PINK "\033[1;35m"
#define RESET "\033[0m"

int MOD = -1;                
int g_n_vars = 0;             
int g_order = 0;         
#define ORDER_LEX      0
#define ORDER_GRLEX    1
#define ORDER_GREVLEX  2
#define ORDER_INVLEX   3
#define ORDER_RINVLEX  4


typedef struct {
    double re, im;
} coef_t;

typedef struct {
    coef_t c;
    int power[MAX_VARS];
} Mono;

typedef struct Trie {
    coef_t coef;          
    int is_end;           
    struct Child {
        int power;
        struct Trie* node;
        struct Child* next;
    }*children;
} Trie;

typedef struct {
    Mono* arr;
    int size, cap;
} List;

void list_init(List* L);
void list_add(List* L, Mono m);
void list_free(List* L);
Trie* create_node();
Trie* get_child(Trie* node, int power);
void insert(Trie* root, int* p, int n, coef_t c);
void collect(Trie* node, int* path, int depth, int n, List* L);
void free_trie(Trie* t);
int mono_equal(Mono a, Mono b, int n);
void normalize(List* L, int n);
int degree_mono(Mono m, int n);
int mono_cmp(const void* a, const void* b);          
void sort_monos(List* L, int n);
void print_coef(coef_t c);
void print_poly(Trie* root, char** vars, int n);
Trie* poly_add(Trie* A, Trie* B, int n);
Trie* poly_sub(Trie* A, Trie* B, int n);
Trie* poly_mul(Trie* A, Trie* B, int n);
Trie* poly_mul_mono(Trie* P, coef_t coeff, int* pow, int n);
int mono_divisible(int* a, int* b, int n);
void mono_div(int* a, int* b, int* q, int n);
int leading_mono(List* L, int n);
void poly_divmod(Trie* A, Trie* B, int n, Trie** Q, Trie** R);
Trie* poly_div(Trie* A, Trie* B, int n);
Trie* poly_mod(Trie* A, Trie* B, int n);
Trie* poly_derivative(Trie* P, int idx, int n);
coef_t eval_trie_complex(Trie* root, coef_t* vals, int n);
int homogeneous(Trie* root, int n);
void decompose(Trie* root, char** vars, int n);
void supp(Trie* root, int n);
int poly_equal(Trie* A, Trie* B, int n);
Trie* get_homogeneous_component(Trie* root, int degree, int n);
void normalize_poly(Trie* p, int n);                 
void leading_term(Trie* p, int n, Mono* out);       
coef_t lc(Trie* p, int n);                      
void multideg(Trie* p, int n, int* deg);         
void lm_str(Trie* p, int n, char* buf);            
void lt_str(Trie* p, int n, char* buf);            



void list_init(List* L) {
    L->size = 0;
    L->cap = 8;
    L->arr = malloc(sizeof(Mono) * L->cap);
}

void list_add(List* L, Mono m) {
    if (L->size == L->cap) {
        L->cap *= 2;
        L->arr = realloc(L->arr, sizeof(Mono) * L->cap);
    }
    L->arr[L->size++] = m;
}

void list_free(List* L) {
    free(L->arr);
}

Trie* create_node() {
    Trie* t = calloc(1, sizeof(Trie));
    return t;
}

Trie* get_child(Trie* node, int power) {
    struct Child* cur = node->children;
    while (cur) {
        if (cur->power == power)
            return cur->node;
        cur = cur->next;
    }
    struct Child* new_child = malloc(sizeof(struct Child));
    new_child->power = power;
    new_child->node = create_node();
    new_child->next = node->children;
    node->children = new_child;
    return new_child->node;
}

void insert(Trie* root, int* p, int n, coef_t c) {
    Trie* cur = root;
    for (int i = 0; i < n; i++) {
        cur = get_child(cur, p[i]);
    }
    cur->is_end = 1;
    cur->coef = (coef_t){ cur->coef.re + c.re, cur->coef.im + c.im };
}

void collect(Trie* node, int* path, int depth, int n, List* L) {
    if (!node) return;
    if (depth > MAX_VARS) return;
    if (node->is_end) {
        Mono m;
        m.c = node->coef;
        memset(m.power, 0, sizeof(m.power));
        if (depth > 0)
            memcpy(m.power, path, sizeof(int) * depth);
        list_add(L, m);
    }
    struct Child* child = node->children;
    while (child) {
        path[depth] = child->power;
        collect(child->node, path, depth + 1, n, L);
        child = child->next;
    }
}

void free_trie(Trie* t) {
    if (!t) return;
    struct Child* child = t->children;
    while (child) {
        free_trie(child->node);
        struct Child* next = child->next;
        free(child);
        child = next;
    }
    free(t);
}

int mono_equal(Mono a, Mono b, int n) {
    for (int i = 0; i < n; i++)
        if (a.power[i] != b.power[i]) return 0;
    return fabs(a.c.re - b.c.re) < EPS && fabs(a.c.im - b.c.im) < EPS;
}

void normalize(List* L, int n) {
    for (int i = 0; i < L->size; i++) {
        for (int j = i + 1; j < L->size; j++) {
            if (mono_equal(L->arr[i], L->arr[j], n)) {
                L->arr[i].c.re += L->arr[j].c.re;
                L->arr[i].c.im += L->arr[j].c.im;
                for (int k = j; k < L->size - 1; k++)
                    L->arr[k] = L->arr[k + 1];
                L->size--;
                j--;
            }
        }
    }
    for (int i = 0; i < L->size; i++) {
        if (fabs(L->arr[i].c.re) < EPS && fabs(L->arr[i].c.im) < EPS) {
            for (int k = i; k < L->size - 1; k++)
                L->arr[k] = L->arr[k + 1];
            L->size--;
            i--;
        }
    }
}

int degree_mono(Mono m, int n) {
    int d = 0;
    for (int i = 0; i < n; i++) d += m.power[i];
    return d;
}

int mono_cmp(const void* a, const void* b) {
    const Mono* ma = (const Mono*)a;
    const Mono* mb = (const Mono*)b;
    int n = g_n_vars;

    int da = degree_mono(*ma, n);
    int db = degree_mono(*mb, n);

    switch (g_order) {
    case ORDER_LEX:
        for (int i = 0; i < n; i++) {
            if (ma->power[i] != mb->power[i])
                return mb->power[i] - ma->power[i];   
        }
        return 0;

    case ORDER_GRLEX:
        if (da != db) return db - da;                
        for (int i = 0; i < n; i++) {
            if (ma->power[i] != mb->power[i])
                return mb->power[i] - ma->power[i];
        }
        return 0;

    case ORDER_GREVLEX:
        if (da != db) return db - da;                 
        for (int i = n - 1; i >= 0; i--) {         
            if (ma->power[i] != mb->power[i])
                return ma->power[i] - mb->power[i]; 
        }
        return 0;

    case ORDER_INVLEX:
    
        for (int i = n - 1; i >= 0; i--) {
            if (ma->power[i] != mb->power[i])
                return mb->power[i] - ma->power[i];   
        }
        return 0;

    case ORDER_RINVLEX:

        for (int i = n - 1; i >= 0; i--) {
            if (ma->power[i] != mb->power[i])
                return ma->power[i] - mb->power[i];
        }
        return 0;

    default:
        return 0;
    }
}

void sort_monos(List* L, int n) {
    g_n_vars = n;
    qsort(L->arr, L->size, sizeof(Mono), mono_cmp);
}

void print_coef(coef_t c) {
    int re_int = (int)(c.re + 0.5);
    int im_int = (int)(c.im + 0.5);
    if (fabs(c.re - re_int) < EPS) c.re = re_int;
    if (fabs(c.im - im_int) < EPS) c.im = im_int;
    if (fabs(c.im) < EPS) {
        if (fabs(c.re) < EPS) printf("0");
        else if (fabs(c.re - 1) < EPS) printf("1");
        else if (fabs(c.re + 1) < EPS) printf("-1");
        else printf("%g", c.re);
        return;
    }
    if (fabs(c.re) < EPS) {
        if (fabs(c.im - 1) < EPS) printf("i");
        else if (fabs(c.im + 1) < EPS) printf("-i");
        else printf("%gi", c.im);
        return;
    }
    if (fabs(c.re - 1) > EPS && fabs(c.re + 1) > EPS) printf("%g", c.re);
    else if (fabs(c.re + 1) < EPS) printf("-");
    if (c.im > 0) {
        if (fabs(c.im - 1) < EPS) printf("+i");
        else printf("+%gi", c.im);
    }
    else {
        if (fabs(c.im + 1) < EPS) printf("-i");
        else printf("%gi", c.im);
    }
}

void print_poly(Trie* root, char** vars, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(root, path, 0, n, &L);
    normalize(&L, n);
    sort_monos(&L, n);
    if (L.size == 0) {
        printf("0\n");
        list_free(&L);
        return;
    }
    int first = 1;
    for (int i = 0; i < L.size; i++) {
        coef_t c = L.arr[i].c;
        int has_vars = 0;
        for (int j = 0; j < n; j++) if (L.arr[i].power[j]) { has_vars = 1; break; }
        if (first) {
            if (c.re < -EPS || (fabs(c.re) < EPS && c.im < -EPS)) {
                printf("-");
                c.re = -c.re;
                c.im = -c.im;
            }
        }
        else {
            if (c.re > EPS || (fabs(c.re) < EPS && c.im > EPS)) {
                printf(" + ");
            }
            else {
                printf(" - ");
                c.re = -c.re;
                c.im = -c.im;
            }
        }
        int coef_is_one = (fabs(c.re - 1) < EPS && fabs(c.im) < EPS);
        if (!(has_vars && coef_is_one)) {
            print_coef(c);
        }
        for (int j = 0; j < n; j++) {
            if (L.arr[i].power[j]) {
                printf("%s", vars[j]);
                if (L.arr[i].power[j] > 1) printf("^%d", L.arr[i].power[j]);
            }
        }
        first = 0;
    }
    printf("\n");
    list_free(&L);
}

void normalize_poly(Trie* p, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(p, path, 0, n, &L);
    normalize(&L, n);

    free_trie(p);

    p->children = NULL;
    p->coef = (coef_t){ 0,0 };
    p->is_end = 0;
    for (int i = 0; i < L.size; i++) {
        insert(p, L.arr[i].power, n, L.arr[i].c);
    }
    list_free(&L);
}

void leading_term(Trie* p, int n, Mono* out) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(p, path, 0, n, &L);
    normalize(&L, n);
    sort_monos(&L, n);
    if (L.size > 0) {
        *out = L.arr[0];
    }
    else {
        memset(out, 0, sizeof(Mono));
        out->c = (coef_t){ 0,0 };
    }
    list_free(&L);
}

coef_t lc(Trie* p, int n) {
    Mono m;
    leading_term(p, n, &m);
    return m.c;
}

void multideg(Trie* p, int n, int* deg) {
    Mono m;
    leading_term(p, n, &m);
    memcpy(deg, m.power, sizeof(int) * n);
}

void lm_str(Trie* p, int n, char* buf) {
    Mono m;
    leading_term(p, n, &m);
    if (m.c.re == 0 && m.c.im == 0) {
        strcpy(buf, "0");
        return;
    }
    char tmp[256] = "";
    for (int i = 0; i < n; i++) {
        if (m.power[i]) {
            char v[32];
            sprintf(v, "x%d", i + 1); 
            strcat(tmp, v);
            if (m.power[i] > 1) {
                char exp[16];
                sprintf(exp, "^%d", m.power[i]);
                strcat(tmp, exp);
            }
        }
    }
    if (strlen(tmp) == 0) strcpy(tmp, "1");
    strcpy(buf, tmp);
}

void lt_str(Trie* p, int n, char* buf) {
    Mono m;
    leading_term(p, n, &m);
    if (m.c.re == 0 && m.c.im == 0) {
        strcpy(buf, "0");
        return;
    }
    char coeff[64];
    if (fabs(m.c.im) < EPS) {
        if (fabs(m.c.re - 1) < EPS) coeff[0] = 0;
        else if (fabs(m.c.re + 1) < EPS) sprintf(coeff, "-");
        else sprintf(coeff, "%g", m.c.re);
    }
    else {
        sprintf(coeff, "(");
        print_coef(m.c);
        sprintf(coeff + strlen(coeff), ")");
    }
    char mon[256] = "";
    for (int i = 0; i < n; i++) {
        if (m.power[i]) {
            char v[32];
            sprintf(v, "x%d", i + 1);
            strcat(mon, v);
            if (m.power[i] > 1) {
                char exp[16];
                sprintf(exp, "^%d", m.power[i]);
                strcat(mon, exp);
            }
        }
    }
    if (strlen(mon) == 0) strcpy(mon, "1");
    if (coeff[0] == 0) sprintf(buf, "%s", mon);
    else sprintf(buf, "%s%s", coeff, mon);
}
Trie* poly_add(Trie* A, Trie* B, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(A, path, 0, n, &L);
    collect(B, path, 0, n, &L);
    normalize(&L, n);
    Trie* res = create_node();
    for (int i = 0; i < L.size; i++)
        insert(res, L.arr[i].power, n, L.arr[i].c);
    list_free(&L);
    normalize_poly(res, n);
    return res;
}

Trie* poly_sub(Trie* A, Trie* B, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(A, path, 0, n, &L);
    List LB;
    list_init(&LB);
    collect(B, path, 0, n, &LB);
    for (int i = 0; i < LB.size; i++) {
        LB.arr[i].c.re = -LB.arr[i].c.re;
        LB.arr[i].c.im = -LB.arr[i].c.im;
        list_add(&L, LB.arr[i]);
    }
    normalize(&L, n);
    Trie* res = create_node();
    for (int i = 0; i < L.size; i++)
        insert(res, L.arr[i].power, n, L.arr[i].c);
    list_free(&L);
    list_free(&LB);
    normalize_poly(res, n);
    return res;
}

Trie* poly_mul(Trie* A, Trie* B, int n) {
    List LA, LB;
    list_init(&LA); list_init(&LB);
    int path[MAX_VARS] = { 0 };
    collect(A, path, 0, n, &LA);
    collect(B, path, 0, n, &LB);
    Trie* res = create_node();
    for (int i = 0; i < LA.size; i++) {
        for (int j = 0; j < LB.size; j++) {
            Mono m;
            m.c.re = LA.arr[i].c.re * LB.arr[j].c.re - LA.arr[i].c.im * LB.arr[j].c.im;
            m.c.im = LA.arr[i].c.re * LB.arr[j].c.im + LA.arr[i].c.im * LB.arr[j].c.re;
            for (int k = 0; k < n; k++)
                m.power[k] = LA.arr[i].power[k] + LB.arr[j].power[k];
            insert(res, m.power, n, m.c);
        }
    }
    list_free(&LA); list_free(&LB);
    normalize_poly(res, n);
    return res;
}

Trie* poly_mul_mono(Trie* P, coef_t coeff, int* pow, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(P, path, 0, n, &L);
    Trie* res = create_node();
    for (int i = 0; i < L.size; i++) {
        Mono m;
        m.c.re = L.arr[i].c.re * coeff.re - L.arr[i].c.im * coeff.im;
        m.c.im = L.arr[i].c.re * coeff.im + L.arr[i].c.im * coeff.re;
        for (int j = 0; j < n; j++)
            m.power[j] = L.arr[i].power[j] + pow[j];
        insert(res, m.power, n, m.c);
    }
    list_free(&L);
    normalize_poly(res, n);
    return res;
}

int mono_divisible(int* a, int* b, int n) {
    for (int i = 0; i < n; i++) if (a[i] < b[i]) return 0;
    return 1;
}

void mono_div(int* a, int* b, int* q, int n) {
    for (int i = 0; i < n; i++) q[i] = a[i] - b[i];
}

int leading_mono(List* L, int n) {
    if (L->size == 0) return -1;
    g_n_vars = n;
    int best = 0;
    for (int i = 1; i < L->size; i++) {
        if (mono_cmp(&L->arr[i], &L->arr[best]) < 0) best = i;
    }
    return best;
}

void poly_divmod(Trie* A, Trie* B, int n, Trie** Q, Trie** R) {
    if (MOD != -1) {
        printf("Division not supported when MOD active.\n");
        *Q = create_node();
        *R = create_node();
        return;
    }
    *Q = create_node();
    *R = create_node();
    int path[MAX_VARS] = { 0 };
    List LR;
    list_init(&LR);
    collect(A, path, 0, n, &LR);
    for (int i = 0; i < LR.size; i++)
        insert(*R, LR.arr[i].power, n, LR.arr[i].c);
    list_free(&LR);

    List LB;
    list_init(&LB);
    collect(B, path, 0, n, &LB);
    normalize(&LB, n);
    if (LB.size == 0) {
        printf("Division by zero.\n");
        list_free(&LB);
        return;
    }

    int safety = 0;
    while (safety++ < 10000) {
        List LR2;
        list_init(&LR2);
        collect(*R, path, 0, n, &LR2);
        normalize(&LR2, n);
        if (LR2.size == 0) {
            list_free(&LR2);
            break;
        }
        int leadR = leading_mono(&LR2, n);
        if (leadR == -1) { list_free(&LR2); break; }
        int* leadR_pow = LR2.arr[leadR].power;
        coef_t leadR_coef = LR2.arr[leadR].c;

        int leadB = leading_mono(&LB, n);
        int* leadB_pow = LB.arr[leadB].power;
        coef_t leadB_coef = LB.arr[leadB].c;

        if (!mono_divisible(leadR_pow, leadB_pow, n)) {
            list_free(&LR2);
            break;
        }

        int qpow[MAX_VARS];
        mono_div(leadR_pow, leadB_pow, qpow, n);
        double den = leadB_coef.re * leadB_coef.re + leadB_coef.im * leadB_coef.im;
        coef_t qcoef = {
            (leadR_coef.re * leadB_coef.re + leadR_coef.im * leadB_coef.im) / den,
            (leadR_coef.im * leadB_coef.re - leadR_coef.re * leadB_coef.im) / den
        };
        insert(*Q, qpow, n, qcoef);
        Trie* qB = poly_mul_mono(B, qcoef, qpow, n);
        Trie* newR = poly_sub(*R, qB, n);
        free_trie(*R);
        free_trie(qB);
        *R = newR;
        list_free(&LR2);
        List test;
        list_init(&test);
        collect(*R, path, 0, n, &test);
        if (test.size == 0) {
            list_free(&test);
            break;
        }
        list_free(&test);
    }
    list_free(&LB);
    normalize_poly(*Q, n);
    normalize_poly(*R, n);
}

Trie* poly_div(Trie* A, Trie* B, int n) {
    Trie* Q, * R;
    poly_divmod(A, B, n, &Q, &R);
    free_trie(R);
    return Q;
}

Trie* poly_mod(Trie* A, Trie* B, int n) {
    Trie* Q, * R;
    poly_divmod(A, B, n, &Q, &R);
    free_trie(Q);
    return R;
}

Trie* poly_derivative(Trie* P, int idx, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(P, path, 0, n, &L);
    Trie* res = create_node();
    for (int i = 0; i < L.size; i++) {
        int exp = L.arr[i].power[idx];
        if (exp == 0) continue;
        coef_t newc = { L.arr[i].c.re * exp, L.arr[i].c.im * exp };
        int newpow[MAX_VARS];
        memcpy(newpow, L.arr[i].power, sizeof(int) * n);
        newpow[idx]--;
        insert(res, newpow, n, newc);
    }
    list_free(&L);
    normalize_poly(res, n);
    return res;
}

coef_t eval_trie_complex(Trie* root, coef_t* vals, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(root, path, 0, n, &L);
    coef_t res = { 0,0 };
    for (int i = 0; i < L.size; i++) {
        coef_t term = L.arr[i].c;
        for (int j = 0; j < n; j++) {
            if (L.arr[i].power[j] > 0) {
                coef_t base = vals[j];
                coef_t pow = { 1,0 };
                int e = L.arr[i].power[j];
                while (e--) {
                    pow.re = pow.re * base.re - pow.im * base.im;
                    pow.im = pow.re * base.im + pow.im * base.re;
                }
                term.re = term.re * pow.re - term.im * pow.im;
                term.im = term.re * pow.im + term.im * pow.re;
            }
        }
        res.re += term.re;
        res.im += term.im;
    }
    list_free(&L);
    return res;
}

int homogeneous(Trie* root, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(root, path, 0, n, &L);
    if (L.size == 0) { list_free(&L); return 1; }
    int deg = degree_mono(L.arr[0], n);
    for (int i = 1; i < L.size; i++)
        if (degree_mono(L.arr[i], n) != deg) { list_free(&L); return 0; }
    list_free(&L);
    return 1;
}

void decompose(Trie* root, char** vars, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(root, path, 0, n, &L);
    normalize(&L, n);
    sort_monos(&L, n);
    for (int d = 0; d <= 20; d++) {
        int printed = 0;
        for (int i = 0; i < L.size; i++) {
            if (degree_mono(L.arr[i], n) == d) {
                if (!printed) { printf("deg %d: ", d); printed = 1; }
                print_coef(L.arr[i].c);
                for (int j = 0; j < n; j++) if (L.arr[i].power[j]) {
                    printf("%s", vars[j]);
                    if (L.arr[i].power[j] > 1) printf("^%d", L.arr[i].power[j]);
                }
                printf(" ");
            }
        }
        if (printed) printf("\n");
    }
    list_free(&L);
}

void supp(Trie* root, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(root, path, 0, n, &L);
    for (int i = 0; i < L.size; i++) {
        printf("(");
        for (int j = 0; j < n; j++) {
            printf("%d", L.arr[i].power[j]);
            if (j < n - 1) printf(",");
        }
        printf(")\n");
    }
    list_free(&L);
}

int poly_equal(Trie* A, Trie* B, int n) {
    List LA, LB;
    list_init(&LA); list_init(&LB);
    int path[MAX_VARS] = { 0 };
    collect(A, path, 0, n, &LA);
    collect(B, path, 0, n, &LB);
    normalize(&LA, n); normalize(&LB, n);
    if (LA.size != LB.size) { list_free(&LA); list_free(&LB); return 0; }
    for (int i = 0; i < LA.size; i++) {
        int found = 0;
        for (int j = 0; j < LB.size; j++)
            if (mono_equal(LA.arr[i], LB.arr[j], n)) { found = 1; break; }
        if (!found) { list_free(&LA); list_free(&LB); return 0; }
    }
    list_free(&LA); list_free(&LB);
    return 1;
}

Trie* get_homogeneous_component(Trie* root, int degree, int n) {
    List L;
    list_init(&L);
    int path[MAX_VARS] = { 0 };
    collect(root, path, 0, n, &L);
    Trie* res = create_node();
    for (int i = 0; i < L.size; i++)
        if (degree_mono(L.arr[i], n) == degree)
            insert(res, L.arr[i].power, n, L.arr[i].c);
    list_free(&L);
    normalize_poly(res, n);
    return res;
}
typedef struct {
    char* s;
    int pos;
    char** vars;
    int n_vars;
} ParserState;

Trie* parse_expression(ParserState* p, int n);
Trie* parse_mul(ParserState* p, int n);
Trie* parse_power(ParserState* p, int n);
Trie* parse_atom(ParserState* p, int n);
coef_t parse_complex(ParserState* p);
double parse_number(ParserState* p);
void skip_spaces(ParserState* p);
int is_alpha(char c);
int is_digit(char c);
int end_of_string(ParserState* p);

int end_of_string(ParserState* p) { return p->s[p->pos] == '\0'; }
void skip_spaces(ParserState* p) { while (!end_of_string(p) && p->s[p->pos] == ' ') p->pos++; }
int is_alpha(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }
int is_digit(char c) { return c >= '0' && c <= '9'; }

double parse_number(ParserState* p) {
    double sign = 1;
    if (p->s[p->pos] == '-') { sign = -1; p->pos++; }
    else if (p->s[p->pos] == '+') p->pos++;
    if (!end_of_string(p) && p->s[p->pos] == 'i') { p->pos++; return sign; }
    double n = 0;
    while (!end_of_string(p) && is_digit(p->s[p->pos])) {
        n = n * 10 + (p->s[p->pos] - '0');
        p->pos++;
    }
    if (!end_of_string(p) && p->s[p->pos] == '.') {
        p->pos++;
        double f = 0.1;
        while (!end_of_string(p) && is_digit(p->s[p->pos])) {
            n += f * (p->s[p->pos] - '0');
            f *= 0.1;
            p->pos++;
        }
    }
    if (!end_of_string(p) && p->s[p->pos] == 'i') {
        p->pos++;
        return sign * n;
    }
    if (MOD != -1) {
        double intpart;
        if (fabs(modf(n, &intpart)) > 1e-9) {
            printf("ERROR: non-integer coefficient not allowed when MOD is active\n");
            exit(1);
        }
    }
    return sign * n;
}

coef_t parse_complex(ParserState* p) {
    if (MOD != -1) {
        double val = parse_number(p);
        if (!end_of_string(p) && p->s[p->pos] == 'i') {
            printf("ERROR: complex numbers not allowed when MOD is active\n");
            exit(1);
        }
        return (coef_t) { val, 0 };
    }
    coef_t c = { 0,0 };
    if (!end_of_string(p) && p->s[p->pos] == '(') {
        p->pos++;
        c = parse_complex(p);
        if (!end_of_string(p) && p->s[p->pos] != ')') {
            printf("ERROR: expected ')'\n");
            exit(1);
        }
        p->pos++;
        return c;
    }
    int start = p->pos;
    double re = parse_number(p);
    if (!end_of_string(p) && p->s[p->pos] == 'i') {
        p->pos++;
        c.im = re;
        return c;
    }
    c.re = re;
    if (!end_of_string(p) && (p->s[p->pos] == '+' || p->s[p->pos] == '-')) {
        int sign_im = (p->s[p->pos] == '-') ? -1 : 1;
        p->pos++;
        int start2 = p->pos;
        double im_val = parse_number(p);
        if (!end_of_string(p) && p->s[p->pos] == 'i') {
            p->pos++;
            c.im = sign_im * im_val;
        }
        else {
            p->pos = start;
            c.re = 0; c.im = 0;
            re = parse_number(p);
            if (!end_of_string(p) && p->s[p->pos] == 'i') {
                p->pos++;
                c.im = re;
            }
            else {
                c.re = re;
            }
        }
    }
    return c;
}

Trie* parse_atom(ParserState* p, int n) {
    skip_spaces(p);
    int sign = 1;
    if (!end_of_string(p) && p->s[p->pos] == '-') { sign = -1; p->pos++; skip_spaces(p); }
    else if (!end_of_string(p) && p->s[p->pos] == '+') { p->pos++; skip_spaces(p); }
    Trie* t = NULL;
    if (!end_of_string(p) && p->s[p->pos] == '(') {
        p->pos++;
        t = parse_expression(p, n);
        if (!end_of_string(p) && p->s[p->pos] != ')') { printf("ERROR: expected ')'\n"); exit(1); }
        p->pos++;
    }
    else if (!end_of_string(p) && is_alpha(p->s[p->pos])) {
        char name[32]; int k = 0;
        while (!end_of_string(p) && is_alpha(p->s[p->pos])) {
            if (k >= 31) { printf("ERROR: variable name too long\n"); exit(1); }
            name[k++] = p->s[p->pos++];
        }
        name[k] = 0;
        int idx = -1;
        for (int i = 0; i < n; i++) if (strcmp(p->vars[i], name) == 0) { idx = i; break; }
        if (idx == -1) { printf("ERROR: unknown variable '%s'\n", name); exit(1); }
        int pow[MAX_VARS] = { 0 };
        pow[idx] = 1;
        t = create_node();
        insert(t, pow, n, (coef_t) { 1, 0 });
    }
    else {
        coef_t c = parse_complex(p);
        int pow[MAX_VARS] = { 0 };
        t = create_node();
        insert(t, pow, n, c);
    }
    if (sign == -1) {
        List L;
        list_init(&L);
        int path[MAX_VARS] = { 0 };
        collect(t, path, 0, n, &L);
        Trie* neg = create_node();
        for (int i = 0; i < L.size; i++) {
            coef_t nc = { -L.arr[i].c.re, -L.arr[i].c.im };
            insert(neg, L.arr[i].power, n, nc);
        }
        list_free(&L);
        free_trie(t);
        t = neg;
    }
    return t;
}

Trie* parse_power(ParserState* p, int n) {
    Trie* left = parse_atom(p, n);
    skip_spaces(p);
    if (!end_of_string(p) && p->s[p->pos] == '^') {
        p->pos++;
        skip_spaces(p);
        double exp = parse_number(p);
        if (exp != (int)exp || exp < 0) {
            printf("ERROR: exponent must be non-negative integer\n");
            exit(1);
        }
        int e = (int)exp;
        Trie* result = create_node();
        int zero[MAX_VARS] = { 0 };
        insert(result, zero, n, (coef_t) { 1, 0 });
        for (int i = 0; i < e; i++) {
            Trie* tmp = poly_mul(result, left, n);
            free_trie(result);
            result = tmp;
        }
        free_trie(left);
        return result;
    }
    return left;
}

Trie* parse_mul(ParserState* p, int n) {
    Trie* left = parse_power(p, n);
    while (1) {
        skip_spaces(p);
        if (end_of_string(p)) break;
        if (p->s[p->pos] == '*') {
            p->pos++;
            Trie* right = parse_power(p, n);
            Trie* prod = poly_mul(left, right, n);
            free_trie(left); free_trie(right);
            left = prod;
        }
        else if (is_alpha(p->s[p->pos]) || p->s[p->pos] == '(' || is_digit(p->s[p->pos])) {
            Trie* right = parse_power(p, n);
            Trie* prod = poly_mul(left, right, n);
            free_trie(left); free_trie(right);
            left = prod;
        }
        else break;
    }
    return left;
}

Trie* parse_expression(ParserState* p, int n) {
    Trie* left = parse_mul(p, n);
    while (1) {
        skip_spaces(p);
        if (end_of_string(p)) break;
        char op = p->s[p->pos];
        if (op == '+' || op == '-') {
            p->pos++;
            Trie* right = parse_mul(p, n);
            Trie* res = (op == '+') ? poly_add(left, right, n) : poly_sub(left, right, n);
            free_trie(left); free_trie(right);
            left = res;
        }
        else break;
    }
    return left;
}

Trie* parse_poly_full(char* s, int n, char** vars) {
    ParserState p = { s, 0, vars, n };
    Trie* res = parse_expression(&p, n);
    skip_spaces(&p);
    if (!end_of_string(&p))
        printf("WARNING: extra characters after expression\n");
    return res;
}
void print_menu() {
    printf("\n" PINK "=== MENU ===\n" RESET);
    printf("1. Enter new polynomials f and g\n");
    printf("2. Show f and g\n");
    printf("3. f + g\n");
    printf("4. f - g\n");
    printf("5. f * g\n");
    if (MOD == -1) {
        printf("6. f / g (quotient)\n");
        printf("7. f mod g (remainder)\n");
    }
    else {
        printf("6. f / g (quotient) [disabled when MOD active]\n");
        printf("7. f mod g (remainder) [disabled when MOD active]\n");
    }
    printf("8. Evaluate f at given point (complex)\n");
    printf("9. Check if f is homogeneous\n");
    printf("10. Decompose f into homogeneous components\n");
    printf("11. Show supp(f)\n");
    printf("12. Check if f == g\n");
    printf("13. Derivative of f (by variable)\n");
    printf("14. Change modulus (MOD)\n");
    printf("15. Show all results (f+g, f-g, f*g, ...)\n");
    printf("16. Extract homogeneous component of f (given degree)\n");
    printf("17. Set monomial order\n");
    printf("18. Show lm, lt, lc, multideg for f and g\n");
    printf("0. Exit\n");
    printf("Choose: ");
}

void show_all(Trie* F, Trie* G, char** vars, int n) {
    if (!F || !G) { printf("No polynomials loaded.\n"); return; }
    printf(PINK "\n**********  COMPLETE ANALYSIS  **********\n" RESET);
    printf(PINK "\n> f = " RESET); print_poly(F, vars, n);
    printf(PINK "> g = " RESET); print_poly(G, vars, n);
    Trie* sum = poly_add(F, G, n);
    printf(PINK "\n> f + g = " RESET); print_poly(sum, vars, n);
    free_trie(sum);
    Trie* diff = poly_sub(F, G, n);
    printf(PINK "> f - g = " RESET); print_poly(diff, vars, n);
    free_trie(diff);
    Trie* prod = poly_mul(F, G, n);
    printf(PINK "> f * g = " RESET); print_poly(prod, vars, n);
    free_trie(prod);
    if (MOD == -1) {
        Trie* quot = poly_div(F, G, n);
        printf(PINK "> f / g (quotient) = " RESET); print_poly(quot, vars, n);
        free_trie(quot);
        Trie* rem = poly_mod(F, G, n);
        printf(PINK "> f mod g (remainder) = " RESET); print_poly(rem, vars, n);
        free_trie(rem);
    }
    else {
        printf(PINK "> Division disabled when MOD is active\n" RESET);
    }
    printf(PINK "\n> f is %shomogeneous\n" RESET, homogeneous(F, n) ? "" : "not ");
    printf(PINK "\n> Homogeneous decomposition of f:\n" RESET);
    decompose(F, vars, n);
    printf(PINK "\n> Support of f:\n" RESET);
    supp(F, n);
    printf(PINK "\n> f == g ? %s\n" RESET, poly_equal(F, G, n) ? "YES" : "NO");
    printf(PINK "\n> Derivatives of f:\n" RESET);
    for (int i = 0; i < n; i++) {
        Trie* df = poly_derivative(F, i, n);
        printf("   d/d%s = ", vars[i]);
        print_poly(df, vars, n);
        free_trie(df);
    }
    printf(PINK "\n**********  END  **********\n" RESET);
}

void set_order() {
    printf("Select monomial order:\n");
    printf("0. lex\n");
    printf("1. grlex\n");
    printf("2. grevlex\n");
    printf("3. invlex\n");
    printf("4. rinvlex\n");
    int ord;
    if (scanf("%d", &ord) != 1) {
        printf("Invalid input.\n");
        while (getchar() != '\n');
        return;
    }
    if (ord >= 0 && ord <= 4) g_order = ord;
    else printf("Invalid order.\n");
    printf("Order set to ");
    switch (g_order) {
    case 0: printf("lex\n"); break;
    case 1: printf("grlex\n"); break;
    case 2: printf("grevlex\n"); break;
    case 3: printf("invlex\n"); break;
    case 4: printf("rinvlex\n"); break;
    }
}

void show_lm_lt_lc_multideg(Trie* P, char** vars, int n, const char* name) {
    if (!P) { printf("%s not loaded.\n", name); return; }
    Mono m;
    leading_term(P, n, &m);
    char mon_str[256], term_str[256];
    lm_str(P, n, mon_str);
    lt_str(P, n, term_str);
    int deg[MAX_VARS];
    multideg(P, n, deg);
    printf("\n%s:\n", name);
    printf("  multideg = (");
    for (int i = 0; i < n; i++) printf("%d%s", deg[i], i == n - 1 ? ")" : ",");
    printf("\n");
    printf("  lm = %s\n", mon_str);
    printf("  lc = "); print_coef(lc(P, n)); printf("\n");
    printf("  lt = %s\n", term_str);
}

int main() {
    int n = 0;
    char** vars = NULL;
    Trie* F = NULL;
    Trie* G = NULL;
    int choice;
    char input[MAX_LEN];
    while (1) {
        print_menu();
        if (fgets(input, MAX_LEN, stdin) == NULL) break;
        if (sscanf(input, "%d", &choice) != 1) {
            printf("Invalid input.\n");
            continue;
        }
        if (choice == 0) break;
        switch (choice) {
        case 1: {
            if (F) free_trie(F);
            if (G) free_trie(G);
            if (vars) {
                for (int i = 0; i < n; i++) free(vars[i]);
                free(vars);
            }
            printf("Number of variables: ");
            if (fgets(input, MAX_LEN, stdin) == NULL) break;
            if (sscanf(input, "%d", &n) != 1 || n <= 0 || n > MAX_VARS) {
                printf("ERROR: invalid number of variables (1..%d)\n", MAX_VARS);
                n = 0;
                continue;
            }
            vars = malloc(sizeof(char*) * n);
            for (int i = 0; i < n; i++) {
                vars[i] = malloc(32);
                printf("Variable %d name: ", i + 1);
                if (fgets(input, MAX_LEN, stdin) == NULL) break;
                sscanf(input, "%s", vars[i]);
            }
            printf("MOD (0 to disable): ");
            if (fgets(input, MAX_LEN, stdin) == NULL) break;
            int mod_input;
            if (sscanf(input, "%d", &mod_input) == 1) {
                MOD = (mod_input == 0) ? -1 : mod_input;
            }
            char a[MAX_LEN], b[MAX_LEN];
            printf("f: ");
            if (fgets(a, MAX_LEN, stdin) == NULL) break;
            a[strcspn(a, "\n")] = 0;
            printf("g: ");
            if (fgets(b, MAX_LEN, stdin) == NULL) break;
            b[strcspn(b, "\n")] = 0;
            F = parse_poly_full(a, n, vars);
            G = parse_poly_full(b, n, vars);
            printf("\nPolynomials loaded.\n");
            break;
        }
        case 2: {
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            printf("f: "); print_poly(F, vars, n);
            printf("g: "); print_poly(G, vars, n);
            break;
        }
        case 3: {
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            Trie* sum = poly_add(F, G, n);
            printf("f+g: "); print_poly(sum, vars, n);
            free_trie(sum);
            break;
        }
        case 4: {
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            Trie* diff = poly_sub(F, G, n);
            printf("f-g: "); print_poly(diff, vars, n);
            free_trie(diff);
            break;
        }
        case 5: {
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            Trie* prod = poly_mul(F, G, n);
            printf("f*g: "); print_poly(prod, vars, n);
            free_trie(prod);
            break;
        }
        case 6: {
            if (MOD != -1) { printf("Division disabled when MOD active.\n"); break; }
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            Trie* quot = poly_div(F, G, n);
            printf("f/g (quotient): "); print_poly(quot, vars, n);
            free_trie(quot);
            break;
        }
        case 7: {
            if (MOD != -1) { printf("Division disabled when MOD active.\n"); break; }
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            Trie* rem = poly_mod(F, G, n);
            printf("f mod g (remainder): "); print_poly(rem, vars, n);
            free_trie(rem);
            break;
        }
        case 8: {
            if (!F) { printf("No f loaded.\n"); break; }
            coef_t* vals = malloc(sizeof(coef_t) * n);
            for (int i = 0; i < n; i++) {
                printf("%s = ", vars[i]);
                char buf[64];
                if (fgets(buf, 64, stdin) == NULL) break;
                buf[strcspn(buf, "\n")] = 0;
                ParserState ps = { buf, 0, vars, n };
                vals[i] = parse_complex(&ps);
            }
            coef_t res = eval_trie_complex(F, vals, n);
            printf("f = ");
            print_coef(res);
            printf("\n");
            free(vals);
            break;
        }
        case 9: {
            if (!F) { printf("No f loaded.\n"); break; }
            printf("f is %shomogeneous.\n", homogeneous(F, n) ? "" : "not ");
            break;
        }
        case 10: {
            if (!F) { printf("No f loaded.\n"); break; }
            decompose(F, vars, n);
            break;
        }
        case 11: {
            if (!F) { printf("No f loaded.\n"); break; }
            supp(F, n);
            break;
        }
        case 12: {
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            printf("f == g ? %s\n", poly_equal(F, G, n) ? "YES" : "NO");
            break;
        }
        case 13: {
            if (!F) { printf("No f loaded.\n"); break; }
            printf("Derivative by variable (0..%d): ", n - 1);
            int idx;
            if (fgets(input, MAX_LEN, stdin) == NULL) break;
            if (sscanf(input, "%d", &idx) != 1 || idx < 0 || idx >= n) {
                printf("Invalid index\n");
                break;
            }
            Trie* df = poly_derivative(F, idx, n);
            printf("d/d%s = ", vars[idx]);
            print_poly(df, vars, n);
            free_trie(df);
            break;
        }
        case 14: {
            printf("MOD (0 to disable): ");
            if (fgets(input, MAX_LEN, stdin) == NULL) break;
            int mod_input;
            if (sscanf(input, "%d", &mod_input) == 1) {
                MOD = (mod_input == 0) ? -1 : mod_input;
            }
            break;
        }
        case 15: {
            show_all(F, G, vars, n);
            break;
        }
        case 16: {
            if (!F) { printf("No f loaded.\n"); break; }
            printf("Enter degree: ");
            int deg;
            if (fgets(input, MAX_LEN, stdin) == NULL) break;
            if (sscanf(input, "%d", &deg) != 1) {
                printf("Invalid degree\n");
                break;
            }
            Trie* h = get_homogeneous_component(F, deg, n);
            printf("Homogeneous component of degree %d: ", deg);
            print_poly(h, vars, n);
            free_trie(h);
            break;
        }
        case 17: {
            set_order();
            break;
        }
        case 18: {
            if (!F || !G) { printf("No polynomials loaded.\n"); break; }
            show_lm_lt_lc_multideg(F, vars, n, "f");
            show_lm_lt_lc_multideg(G, vars, n, "g");
            break;
        }
        default:
            printf("Invalid choice.\n");
        }
    }
    if (F) free_trie(F);
    if (G) free_trie(G);
    if (vars) {
        for (int i = 0; i < n; i++) free(vars[i]);
        free(vars);
    }
    return 0;
}