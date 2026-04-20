#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <complex>

class rational {
public:
    using integer = long long;
    rational() : num_(0), den_(1) {}
    rational(integer n) : num_(n), den_(1) {}
    rational(integer n, integer d) : num_(n), den_(d) {
        if (den_ == 0) throw std::invalid_argument("Zero denominator");
        norm();
    }
    integer num() const { return num_; }
    integer den() const { return den_; }
    rational operator-() const { return rational(-num_, den_); }
    rational operator+(const rational& o) const { return rational(num_ * o.den_ + o.num_ * den_, den_ * o.den_); }
    rational operator-(const rational& o) const { return rational(num_ * o.den_ - o.num_ * den_, den_ * o.den_); }
    rational operator*(const rational& o) const { return rational(num_ * o.num_, den_ * o.den_); }
    rational operator/(const rational& o) const {
        if (o.num_ == 0) throw std::invalid_argument("Division by zero");
        return rational(num_ * o.den_, den_ * o.num_);
    }
    bool operator==(const rational& o) const { return num_ == o.num_ && den_ == o.den_; }
    bool operator!=(const rational& o) const { return !(*this == o); }
    bool operator<(const rational& o) const { return num_ * o.den_ < o.num_ * den_; }
    bool operator>(const rational& o) const { return o < *this; }
    bool operator<=(const rational& o) const { return !(o < *this); }
    bool operator>=(const rational& o) const { return !(*this < o); }
    explicit operator double() const { return (double)num_ / den_; }
private:
    integer num_, den_;
    void norm() {
        if (den_ < 0) { num_ = -num_; den_ = -den_; }
        integer g = std::gcd(num_ < 0 ? -num_ : num_, den_);
        num_ /= g; den_ /= g;
    }
};
std::ostream& operator<<(std::ostream& os, const rational& r) {
    if (r.den() == 1) os << r.num();
    else os << r.num() << "/" << r.den();
    return os;
}
rational abs(const rational& r) { return r.num() < 0 ? -r : r; }
template<typename T> struct coeff_traits {
    static bool is_zero(const T& v) { return v == T(0); }
    static T abs(const T& v) { return v < 0 ? -v : v; }
};
template<> struct coeff_traits<double> {
    static bool is_zero(double v) {
        return std::abs(v) <= 1e-12 * std::max(1.0, std::abs(v));
    }
    static double abs(double v) { return std::abs(v); }
};
template<> struct coeff_traits<rational> {
    static bool is_zero(const rational& v) { return v == rational(0); }
    static rational abs(const rational& v) { return ::abs(v); }
};
template<typename T> struct coeff_traits<std::complex<T>> {
    static bool is_zero(const std::complex<T>& v) {
        return coeff_traits<T>::is_zero(v.real()) && coeff_traits<T>::is_zero(v.imag());
    }
    static std::complex<T> abs(const std::complex<T>& v) {
        return std::complex<T>(std::abs(v));
    }
};
enum class ordering_type { lex, grlex, grevlex, invlex, rinvlex };
template<typename coeff>
class poly_ring {
public:
    using degs = std::vector<int>;
    poly_ring(std::vector<std::string> vars, ordering_type ord) : vars_(vars), ord_(ord) {}
    const std::vector<std::string>& vars() const { return vars_; }
    size_t nvars() const { return vars_.size(); }
    ordering_type ord() const { return ord_; }
    int cmp(const degs& a, const degs& b) const {
        switch (ord_) {
        case ordering_type::lex: return cmp_lex(a, b);
        case ordering_type::grlex: return cmp_grlex(a, b);
        case ordering_type::grevlex: return cmp_grevlex(a, b);
        case ordering_type::invlex: return cmp_invlex(a, b);
        case ordering_type::rinvlex: return cmp_rinvlex(a, b);
        }
        return 0;
    }
    static bool divides(const degs& a, const degs& b) {
        for (size_t i = 0; i < a.size(); ++i) if (a[i] < b[i]) return false;
        return true;
    }
    static degs lcm(const degs& a, const degs& b) {
        degs r(a.size());
        for (size_t i = 0; i < a.size(); ++i) r[i] = std::max(a[i], b[i]);
        return r;
    }
    static int total_deg(const degs& d) {
        int s = 0; for (int x : d) s += x; return s;
    }
private:
    std::vector<std::string> vars_;
    ordering_type ord_;
    static int cmp_lex(const degs& a, const degs& b) {
        for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return a[i] - b[i];
        return 0;
    }
    static int cmp_grlex(const degs& a, const degs& b) {
        int sa = total_deg(a), sb = total_deg(b);
        if (sa != sb) return sa - sb;
        return cmp_lex(a, b);
    }
    static int cmp_grevlex(const degs& a, const degs& b) {
        int sa = total_deg(a), sb = total_deg(b);
        if (sa != sb) return sa - sb;
        for (int i = (int)a.size() - 1; i >= 0; --i) if (a[i] != b[i]) return b[i] - a[i];
        return 0;
    }
    static int cmp_invlex(const degs& a, const degs& b) {
        for (int i = (int)a.size() - 1; i >= 0; --i) if (a[i] != b[i]) return a[i] - b[i];
        return 0;
    }
    static int cmp_rinvlex(const degs& a, const degs& b) {
        for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return b[i] - a[i];
        return 0;
    }
};
template<typename coeff>
class polynomial {
public:
    using degs = std::vector<int>;
    using ring_ptr = std::shared_ptr<poly_ring<coeff>>;
    using traits = coeff_traits<coeff>;
private:
    struct node {
        coeff coeff = coeff(0);
        std::map<int, std::unique_ptr<node>> children;
    };
    ring_ptr ring_;
    std::unique_ptr<node> root_;
    void check_ring(const polynomial& o) const { if (ring_ != o.ring_) throw std::invalid_argument("different rings"); }
    static std::unique_ptr<node> clone(const node* src) {
        if (!src) return nullptr;
        auto n = std::make_unique<node>();
        n->coeff = src->coeff;
        for (auto& p : src->children) n->children[p.first] = clone(p.second.get());
        return n;
    }
    void insert(node* node, size_t depth, const degs& degs, coeff c) {
        if (depth == nvars()) {
            node->coeff = node->coeff + c;
            return;
        }
        int d = degs[depth];
        auto& child = node->children[d];
        if (!child) child = std::make_unique<node>();
        insert(child.get(), depth + 1, degs, c);
    }
    static bool get_coeff(const node* node, size_t depth, const degs& degs, coeff& out) {
        if (!node) return false;
        if (depth == degs.size()) { out = node->coeff; return true; }
        auto it = node->children.find(degs[depth]);
        if (it == node->children.end()) return false;
        return get_coeff(it->second.get(), depth + 1, degs, out);
    }
    static void add_trees(node* res, const node* a, const node* b, size_t depth) {
        if (!a && !b) return;
        if (a) res->coeff = res->coeff + a->coeff;
        if (b) res->coeff = res->coeff + b->coeff;
        std::set<int> degs;
        if (a) for (auto& p : a->children) degs.insert(p.first);
        if (b) for (auto& p : b->children) degs.insert(p.first);
        for (int d : degs) {
            const node* ca = nullptr, * cb = nullptr;
            if (a) { auto it = a->children.find(d); if (it != a->children.end()) ca = it->second.get(); }
            if (b) { auto it = b->children.find(d); if (it != b->children.end()) cb = it->second.get(); }
            auto& rc = res->children[d];
            if (!rc) rc = std::make_unique<node>();
            add_trees(rc.get(), ca, cb, depth + 1);
        }
        prune(res, depth);
    }
    static void sub_trees(node* res, const node* a, const node* b, size_t depth) {
        if (!a && !b) return;
        if (a) res->coeff = res->coeff + a->coeff;
        if (b) res->coeff = res->coeff - b->coeff;
        std::set<int> degs;
        if (a) for (auto& p : a->children) degs.insert(p.first);
        if (b) for (auto& p : b->children) degs.insert(p.first);
        for (int d : degs) {
            const node* ca = nullptr, * cb = nullptr;
            if (a) { auto it = a->children.find(d); if (it != a->children.end()) ca = it->second.get(); }
            if (b) { auto it = b->children.find(d); if (it != b->children.end()) cb = it->second.get(); }
            auto& rc = res->children[d];
            if (!rc) rc = std::make_unique<node>();
            sub_trees(rc.get(), ca, cb, depth + 1);
        }
        prune(res, depth);
    }
    static void mul_scalar(node* res, const node* src, size_t depth, coeff s) {
        if (!src) return;
        res->coeff = src->coeff * s;
        for (auto& p : src->children) {
            auto& rc = res->children[p.first];
            if (!rc) rc = std::make_unique<node>();
            mul_scalar(rc.get(), p.second.get(), depth + 1, s);
        }
        prune(res, depth);
    }
    static bool trees_eq(const node* a, const node* b) {
        if (!a && !b) return true;
        if (!a || !b) return false;
        if (!traits::is_zero(a->coeff - b->coeff)) return false;
        if (a->children.size() != b->children.size()) return false;
        for (auto& p : a->children) {
            auto it = b->children.find(p.first);
            if (it == b->children.end()) return false;
            if (!trees_eq(p.second.get(), it->second.get())) return false;
        }
        return true;
    }
    static bool is_zero_node(const node* node) {
        if (!node) return true;
        if (!traits::is_zero(node->coeff)) return false;
        for (auto& p : node->children) if (!is_zero_node(p.second.get())) return false;
        return true;
    }
    static void prune(node* node, size_t depth) {
        if (!node) return;
        for (auto it = node->children.begin(); it != node->children.end(); ) {
            prune(it->second.get(), depth + 1);
            if (is_zero_node(it->second.get()))
                it = node->children.erase(it);
            else
                ++it;
        }
    }
    void trav(const node* node, size_t depth, degs& cur,
        const std::function<void(const degs&, coeff)>& f) const {
        if (!node) return;
        if (!traits::is_zero(node->coeff)) f(cur, node->coeff);
        if (depth == nvars()) return;
        for (auto& p : node->children) {
            cur[depth] = p.first;
            trav(p.second.get(), depth + 1, cur, f);
        }
        cur[depth] = 0;
    }
    void shift(const degs& s) {
        polynomial sh(ring_);
        traverse([&](const degs& d, coeff c) {
            degs nd = d;
            for (size_t i = 0; i < nd.size(); ++i) nd[i] += s[i];
            sh.add_monomial(nd, c);
            });
        *this = std::move(sh);
    }
    degs lead_deg() const {
        degs best(nvars(), 0);
        bool first = true;
        traverse([&](const degs& d, coeff c) {
            if (traits::is_zero(c)) return;
            if (first) { best = d; first = false; }
            else if (ring_->cmp(d, best) > 0) best = d;
            });
        return best;
    }
public:
    explicit polynomial(ring_ptr ring) : ring_(ring), root_(std::make_unique<node>()) {}
    polynomial(const polynomial& o) : ring_(o.ring_), root_(clone(o.root_.get())) {}
    polynomial(polynomial&&) = default;
    polynomial& operator=(const polynomial& o) {
        if (this != &o) { ring_ = o.ring_; root_ = clone(o.root_.get()); }
        return *this;
    }
    const poly_ring<coeff>& ring() const { return *ring_; }
    size_t nvars() const { return ring_->nvars(); }
    void add_monomial(const degs& d, coeff c) {
        if (traits::is_zero(c)) return;
        insert(root_.get(), 0, d, c);
        prune(root_.get(), 0);
    }
    polynomial operator+(const polynomial& o) const {
        check_ring(o);
        polynomial r(ring_);
        add_trees(r.root_.get(), root_.get(), o.root_.get(), 0);
        prune(r.root_.get(), 0);
        return r;
    }
    polynomial operator-(const polynomial& o) const {
        check_ring(o);
        polynomial r(ring_);
        sub_trees(r.root_.get(), root_.get(), o.root_.get(), 0);
        prune(r.root_.get(), 0);
        return r;
    }
    polynomial operator*(const polynomial& o) const {
        check_ring(o);
        polynomial r(ring_);
        traverse([&](const degs& d, coeff c) {
            polynomial t = o * c;
            t.shift(d);
            r = r + t;
            });
        prune(r.root_.get(), 0);
        return r;
    }
    polynomial operator*(coeff s) const {
        polynomial r(ring_);
        if (traits::is_zero(s)) return r;
        mul_scalar(r.root_.get(), root_.get(), 0, s);
        prune(r.root_.get(), 0);
        return r;
    }
    polynomial operator-() const { return *this * coeff(-1); }
    bool operator==(const polynomial& o) const {
        check_ring(o);
        return trees_eq(root_.get(), o.root_.get());
    }
    bool operator!=(const polynomial& o) const { return !(*this == o); }
    bool is_zero() const { return is_zero_node(root_.get()); }
    std::set<degs> support() const {
        std::set<degs> s;
        traverse([&](const degs& d, coeff) { s.insert(d); });
        return s;
    }
    coeff eval(const std::vector<coeff>& pt) const {
        coeff val = coeff(0);
        traverse([&](const degs& d, coeff c) {
            coeff term = c;
            for (size_t i = 0; i < d.size(); ++i) {
                if (d[i] == 0) continue;
                coeff pow = coeff(1);
                for (int k = 0; k < d[i]; ++k) pow = pow * pt[i];
                term = term * pow;
            }
            val = val + term;
            });
        return val;
    }
    bool is_homog() const {
        int td = -1;
        bool hom = true;
        traverse([&](const degs& d, coeff) {
            int deg = poly_ring<coeff>::total_deg(d);
            if (td == -1) td = deg;
            else if (td != deg) hom = false;
            });
        return hom;
    }
    int homog_deg() const {
        int td = -1;
        bool hom = true;
        traverse([&](const degs& d, coeff) {
            int deg = poly_ring<coeff>::total_deg(d);
            if (td == -1) td = deg;
            else if (td != deg) hom = false;
            });
        return hom ? td : -1;
    }
    polynomial homog_part(int d) const {
        polynomial p(ring_);
        traverse([&](const degs& degs, coeff c) {
            if (poly_ring<coeff>::total_deg(degs) == d) p.add_monomial(degs, c);
            });
        return p;
    }
    degs multideg() const { return lead_deg(); }
    coeff lc() const {
        auto deg = lead_deg();
        coeff c = coeff(0);
        get_coeff(root_.get(), 0, deg, c);
        return c;
    }
    polynomial lm() const {
        polynomial m(ring_);
        m.add_monomial(lead_deg(), coeff(1));
        return m;
    }
    polynomial lt() const {
        polynomial t(ring_);
        t.add_monomial(lead_deg(), lc());
        return t;
    }
    void traverse(std::function<void(const degs&, coeff)> f) const {
        degs cur(nvars(), 0);
        trav(root_.get(), 0, cur, f);
    }
    static polynomial s_poly(const polynomial& f, const polynomial& g) {
        f.check_ring(g);
        coeff lcf = f.lc(), lcg = g.lc();
        if (traits::is_zero(lcf) || traits::is_zero(lcg))
            throw std::invalid_argument("s_poly: zero leading coefficient");
        auto df = f.multideg(), dg = g.multideg();
        auto L = poly_ring<coeff>::lcm(df, dg);
        polynomial m1(f.ring_);
        degs d1 = L; for (size_t i = 0; i < d1.size(); ++i) d1[i] -= df[i];
        m1.add_monomial(d1, coeff(1) / lcf);
        polynomial m2(g.ring_);
        degs d2 = L; for (size_t i = 0; i < d2.size(); ++i) d2[i] -= dg[i];
        m2.add_monomial(d2, coeff(1) / lcg);
        return m1 * f - m2 * g;
    }
    static polynomial reduce(const polynomial& p, const std::vector<polynomial>& basis) {
        if (basis.empty()) return p;
        polynomial rem(p.ring_);
        polynomial cur = p;
        while (!cur.is_zero()) {
            auto lt_cur = cur.lt();
            if (lt_cur.is_zero()) break;
            auto dcur = cur.multideg();
            bool divided = false;
            for (const auto& g : basis) {
                if (g.is_zero()) continue;
                auto dg = g.multideg();
                if (poly_ring<coeff>::divides(dcur, dg)) {
                    degs diff = dcur;
                    for (size_t i = 0; i < diff.size(); ++i) diff[i] -= dg[i];
                    coeff factor = lt_cur.lc() / g.lc();
                    polynomial sub(p.ring_);
                    sub.add_monomial(diff, factor);
                    cur = cur - sub * g;
                    divided = true;
                    break;
                }
            }
            if (!divided) {
                rem = rem + lt_cur;
                cur = cur - lt_cur;
            }
        }
        return rem;
    }
    static std::vector<polynomial> buchberger(const std::vector<polynomial>& gens) {
        if (gens.empty()) return {};
        std::vector<polynomial> G = gens;
        G.erase(std::remove_if(G.begin(), G.end(), [](const polynomial& p) { return p.is_zero(); }), G.end());
        if (G.empty()) return {};
        using pair = std::pair<size_t, size_t>;
        std::multimap<int, pair> B;
        auto add_pairs = [&](size_t newIdx) {
            for (size_t i = 0; i < newIdx; ++i) {
                auto lcm = poly_ring<coeff>::lcm(G[i].multideg(), G[newIdx].multideg());
                int deg = poly_ring<coeff>::total_deg(lcm);
                B.emplace(deg, pair(i, newIdx));
            }
            };
        for (size_t i = 0; i < G.size(); ++i)
            for (size_t j = i + 1; j < G.size(); ++j) {
                auto lcm = poly_ring<coeff>::lcm(G[i].multideg(), G[j].multideg());
                int deg = poly_ring<coeff>::total_deg(lcm);
                B.emplace(deg, pair(i, j));
            }
        while (!B.empty()) {
            auto it = B.begin();
            size_t i = it->second.first, j = it->second.second;
            B.erase(it);
            polynomial s = s_poly(G[i], G[j]);
            polynomial r = reduce(s, G);
            if (!r.is_zero()) {
                size_t newIdx = G.size();
                G.push_back(r);
                add_pairs(newIdx);
            }
        }
        return G;
    }
    static std::vector<polynomial> reduce_groebner(std::vector<polynomial> G) {
        for (auto& g : G) {
            coeff lc = g.lc();
            if (!traits::is_zero(lc) && !traits::is_zero(lc - coeff(1)))
                g = g * (coeff(1) / lc);
        }
        bool changed;
        do {
            changed = false;
            std::vector<polynomial> H;
            for (size_t i = 0; i < G.size(); ++i) {
                bool keep = true;
                for (size_t j = 0; j < G.size(); ++j) {
                    if (i == j) continue;
                    if (poly_ring<coeff>::divides(G[i].multideg(), G[j].multideg())) {
                        keep = false;
                        break;
                    }
                }
                if (keep) H.push_back(G[i]);
            }
            if (H.size() != G.size()) {
                G = std::move(H);
                changed = true;
            }
            for (size_t i = 0; i < G.size(); ++i) {
                std::vector<polynomial> others;
                for (size_t j = 0; j < G.size(); ++j)
                    if (j != i) others.push_back(G[j]);
                polynomial r = reduce(G[i], others);
                coeff lc_r = r.lc();
                if (!traits::is_zero(lc_r) && !traits::is_zero(lc_r - coeff(1)))
                    r = r * (coeff(1) / lc_r);
                if (!r.is_zero() && !(r == G[i])) {
                    G[i] = r;
                    changed = true;
                }
            }
            G.erase(std::remove_if(G.begin(), G.end(),
                [](const polynomial& p) { return p.is_zero(); }), G.end());
        } while (changed);
        return G;
    }
    static bool is_groebner_basis(const std::vector<polynomial>& G) {
        for (size_t i = 0; i < G.size(); ++i) {
            for (size_t j = i + 1; j < G.size(); ++j) {
                polynomial s = s_poly(G[i], G[j]);
                polynomial r = reduce(s, G);
                if (!r.is_zero()) return false;
            }
        }
        return true;
    }
};
template<typename coeff>
std::ostream& operator<<(std::ostream& os, const polynomial<coeff>& p) {
    if (p.is_zero()) return os << "0";
    bool first = true;
    p.traverse([&](const std::vector<int>& d, coeff c) {
        if (coeff_traits<coeff>::is_zero(c)) return;
        if (!first) os << (c > coeff(0) ? " + " : " - ");
        else if (c < coeff(0)) os << "-";
        auto ac = coeff_traits<coeff>::abs(c);
        bool isConst = true;
        for (int x : d) if (x) { isConst = false; break; }
        if (isConst || !coeff_traits<coeff>::is_zero(ac - coeff(1))) os << ac;
        for (size_t i = 0; i < d.size(); ++i)
            if (d[i] > 0) { os << p.ring().vars()[i]; if (d[i] != 1) os << "^" << d[i]; }
        first = false;
        });
    return os;
}
template<typename coeff> coeff parse_coeff(const std::string& s);
template<> int parse_coeff<int>(const std::string& s) { return std::stoi(s); }
template<> double parse_coeff<double>(const std::string& s) { return std::stod(s); }
template<> rational parse_coeff<rational>(const std::string& s) {
    size_t slash = s.find('/');
    if (slash == std::string::npos) return rational(std::stoll(s));
    return rational(std::stoll(s.substr(0, slash)), std::stoll(s.substr(slash + 1)));
}
template<typename coeff>
polynomial<coeff> parse_polynomial(const std::string& line, std::shared_ptr<poly_ring<coeff>> ring) {
    polynomial<coeff> p(ring);
    std::string s = line;
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
    if (s.empty() || s == "0") return p;
    size_t pos = 0;
    while (pos < s.length()) {
        coeff sign = 1;
        if (s[pos] == '+') pos++;
        else if (s[pos] == '-') { sign = -1; pos++; }
        size_t start = pos;
        while (pos < s.length() && s[pos] != '+' && s[pos] != '-') pos++;
        std::string term = s.substr(start, pos - start);
        if (term.empty()) continue;
        std::vector<std::string> factors;
        size_t last = 0;
        for (size_t i = 0; i <= term.length(); ++i)
            if (i == term.length() || term[i] == '*') {
                factors.push_back(term.substr(last, i - last));
                last = i + 1;
            }
        coeff coeff = 1;
        typename polynomial<coeff>::degs degs(ring->nvars(), 0);
        for (auto& f : factors) {
            if (f.empty()) continue;
            bool isNum = true;
            for (size_t k = 0; k < f.size(); ++k) {
                char ch = f[k];
                if (k == 0 && ch == '-') continue;
                if (!isdigit(ch) && ch != '/' && ch != '.') { isNum = false; break; }
            }
            if (isNum) coeff = coeff * parse_coeff<coeff>(f);
            else {
                size_t caret = f.find('^');
                std::string vname;
                int exp = 1;
                if (caret == std::string::npos) vname = f;
                else { vname = f.substr(0, caret); exp = std::stoi(f.substr(caret + 1)); }
                int idx = -1;
                for (size_t j = 0; j < ring->vars().size(); ++j)
                    if (ring->vars()[j] == vname) { idx = (int)j; break; }
                if (idx == -1) throw std::runtime_error("Unknown variable: " + vname);
                degs[idx] += exp;
            }
        }
        p.add_monomial(degs, coeff * sign);
    }
    return p;
}
template<typename coeff>
void run() {
    std::cout << "Enter number of variables: "; size_t n; std::cin >> n; std::cin.ignore();
    std::vector<std::string> vars(n);
    std::cout << "Enter variable names separated by space: ";
    for (size_t i = 0; i < n; ++i) std::cin >> vars[i]; std::cin.ignore();
    std::cout << "Ordering (0:lex, 1:grlex, 2:grevlex, 3:invlex, 4:rinvlex): ";
    int o; std::cin >> o; std::cin.ignore();
    ordering_type ord = static_cast<ordering_type>(o);
    auto ring = std::make_shared<poly_ring<coeff>>(vars, ord);
    std::cout << "Enter generators (one per line, empty line to finish):\n";
    std::vector<polynomial<coeff>> gens;
    while (true) {
        std::string line; std::getline(std::cin, line);
        if (line.empty()) break;
        gens.push_back(parse_polynomial<coeff>(line, ring));
    }
    std::cout << "\nGenerators:\n";
    for (auto& g : gens) std::cout << g << "\n";
    if (!gens.empty()) {
        auto& f = gens[0];
        std::cout << "\n    Information about f    \n";
        std::cout << "supp(f) = { ";
        for (auto& d : f.support()) {
            std::cout << "("; for (size_t i = 0; i < d.size(); ++i) std::cout << d[i] << (i + 1 < d.size() ? "," : "");
            std::cout << ") ";
        }
        std::cout << "}\n";
        std::cout << "Homogeneous? " << (f.is_homog() ? "yes" : "no") << ", degree = " << f.homog_deg() << "\n";
        std::cout << "lm(f) = " << f.lm() << "\nlt(f) = " << f.lt() << "\nlc(f) = " << f.lc() << "\n";
        std::cout << "multideg(f) = ("; auto md = f.multideg();
        for (size_t i = 0; i < md.size(); ++i) std::cout << md[i] << (i + 1 < md.size() ? "," : "");
        std::cout << ")\n";
    }
    if (gens.size() >= 2) {
        auto S = polynomial<coeff>::s_poly(gens[0], gens[1]);
        std::cout << "\nS(f1, f2) = " << S << "\n";
        std::cout << "multideg(S) = ("; auto mdS = S.multideg();
        for (size_t i = 0; i < mdS.size(); ++i) std::cout << mdS[i] << (i + 1 < mdS.size() ? "," : "");
        std::cout << ")\n";
    }
    bool isGB = polynomial<coeff>::is_groebner_basis(gens);
    std::cout << "\nOriginal set " << (isGB ? "IS" : "IS NOT") << " a Groebner basis.\n";
    auto G = polynomial<coeff>::buchberger(gens);
    std::cout << "\nGroebner basis (unreduced, " << G.size() << " elements):\n";
    for (auto& g : G) std::cout << g << "\n";
    auto redG = polynomial<coeff>::reduce_groebner(G);
    std::cout << "\nReduced Groebner basis:\n";
    for (auto& g : redG) std::cout << g << "\n";
    std::cout << "\nEvaluate a polynomial at a point? (y/n): ";
    char ans; std::cin >> ans; std::cin.ignore();
    if (ans == 'y' || ans == 'Y') {
        std::cout << "Enter point coordinates separated by space: ";
        std::vector<coeff> pt(n);
        for (size_t i = 0; i < n; ++i) {
            std::string val; std::cin >> val;
            pt[i] = parse_coeff<coeff>(val);
        }
        if (!gens.empty()) {
            coeff val = gens[0].eval(pt);
            std::cout << "f(point) = " << val << "\n";
        }
    }
}
int main() {
    std::cout << "    Polynomial Manipulation and Buchberger's Algorithm    \n";
    std::cout << "Choose coefficient type:\n";
    std::cout << "1: int (WARNING: integer division may break Groebner basis)\n";
    std::cout << "2: double (may be unstable due to rounding)\n";
    std::cout << "3: Rational (recommended)\n";
    std::cout << "Your choice: ";
    int typeChoice; std::cin >> typeChoice; std::cin.ignore();
    if (typeChoice == 1) {
        std::cout << "\n!!! Warning: int division is exact only for factors. Use Rational for reliability.\n";
        std::cout << "Continue? (y/n): ";
        char c; std::cin >> c; std::cin.ignore();
        if (c == 'y' || c == 'Y')
            run<int>();
        else
            std::cout << "Exiting.\n";
    }
    else if (typeChoice == 2) {
        std::cout << "\nNote: floating-point arithmetic may cause inaccuracies.\n";
        run<double>();
    }
    else if (typeChoice == 3) {
        run<rational>();
    }
    else {
        std::cout << "Invalid choice.\n";
    }
    return 0;
}