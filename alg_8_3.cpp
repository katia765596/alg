#include <algorithm>
#include <cmath>
#include <cctype>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

enum class OrderingType { lex, grlex, grevlex, invlex, rinvlex };

class PolynomialRing {
public:
    using DegreeVector = std::vector<int>;

    PolynomialRing(std::vector<std::string> vars, OrderingType ord)
        : vars_(std::move(vars)), ord_(ord) {
    }

    const std::vector<std::string>& vars() const { return vars_; }
    size_t nvars() const { return vars_.size(); }
    OrderingType ord() const { return ord_; }

    int cmp(const DegreeVector& a, const DegreeVector& b) const {
        if (a.size() != nvars() || b.size() != nvars())
            throw std::invalid_argument("wrong size");
        switch (ord_) {
        case OrderingType::lex:   return cmpLex(a, b);
        case OrderingType::grlex: return cmpGrlex(a, b);
        case OrderingType::grevlex: return cmpGrevlex(a, b);
        case OrderingType::invlex: return cmpInvlex(a, b);
        case OrderingType::rinvlex: return cmpRinvlex(a, b);
        }
        return 0;
    }

    static bool divides(const DegreeVector& a, const DegreeVector& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i)
            if (a[i] < b[i]) return false;
        return true;
    }

    static DegreeVector lcm(const DegreeVector& a, const DegreeVector& b) {
        DegreeVector res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = std::max(a[i], b[i]);
        return res;
    }

    static int totalDeg(const DegreeVector& d) {
        return std::accumulate(d.begin(), d.end(), 0);
    }

private:
    std::vector<std::string> vars_;
    OrderingType ord_;

    static int cmpLex(const DegreeVector& a, const DegreeVector& b) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) return a[i] - b[i];
        }
        return 0;
    }

    static int cmpGrlex(const DegreeVector& a, const DegreeVector& b) {
        int sa = totalDeg(a), sb = totalDeg(b);
        if (sa != sb) return sa - sb;
        return cmpLex(a, b);
    }

    static int cmpGrevlex(const DegreeVector& a, const DegreeVector& b) {
        int sa = totalDeg(a), sb = totalDeg(b);
        if (sa != sb) return sa - sb;
        for (int i = (int)a.size() - 1; i >= 0; --i) {
            if (a[i] != b[i]) return b[i] - a[i];
        }
        return 0;
    }

    static int cmpInvlex(const DegreeVector& a, const DegreeVector& b) {
        for (int i = (int)a.size() - 1; i >= 0; --i) {
            if (a[i] != b[i]) return a[i] - b[i];
        }
        return 0;
    }

    static int cmpRinvlex(const DegreeVector& a, const DegreeVector& b) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) return b[i] - a[i];
        }
        return 0;
    }
};

class Polynomial {
public:
    using Coeff = double;
    using DegreeVector = std::vector<int>;
    using RingPtr = std::shared_ptr<PolynomialRing>;

    static constexpr Coeff EPS = 1e-12;

    explicit Polynomial(RingPtr ring) : ring_(ring), root_(std::make_unique<Node>()) {}

    Polynomial(const Polynomial& other) : ring_(other.ring_), root_(cloneNode(other.root_.get())) {}
    Polynomial(Polynomial&&) noexcept = default;
    Polynomial& operator=(const Polynomial& other) {
        if (this != &other) { ring_ = other.ring_; root_ = cloneNode(other.root_.get()); }
        return *this;
    }
    Polynomial& operator=(Polynomial&&) noexcept = default;

    const PolynomialRing& ring() const { return *ring_; }
    size_t nvars() const { return ring_->nvars(); }

    void addMonomial(const DegreeVector& degs, Coeff c) {
        if (degs.size() != nvars()) throw std::invalid_argument("wrong size");
        if (std::abs(c) < EPS) return;
        insert(root_.get(), 0, degs, c);
    }

    Polynomial operator+(const Polynomial& other) const {
        checkSameRing(other);
        Polynomial res(ring_);
        addTrees(res.root_.get(), root_.get(), other.root_.get(), 0);
        return res;
    }

    Polynomial operator-(const Polynomial& other) const {
        checkSameRing(other);
        Polynomial res(ring_);
        subTrees(res.root_.get(), root_.get(), other.root_.get(), 0);
        return res;
    }

    Polynomial operator*(const Polynomial& other) const {
        checkSameRing(other);
        Polynomial res(ring_);
        traverse([&](const DegreeVector& d, Coeff c) {
            Polynomial t = other * c;
            t.shift(d);
            res = res + t;
            });
        return res;
    }

    Polynomial operator*(Coeff s) const {
        Polynomial res(ring_);
        if (std::abs(s) < EPS) return res;
        mulScalar(res.root_.get(), root_.get(), 0, s);
        return res;
    }

    Polynomial operator-() const { return *this * (-1.0); }

    bool operator==(const Polynomial& other) const {
        checkSameRing(other);
        return treesEq(root_.get(), other.root_.get(), 0);
    }
    bool operator!=(const Polynomial& other) const { return !(*this == other); }

    std::set<DegreeVector> support() const {
        std::set<DegreeVector> s;
        traverse([&](const DegreeVector& d, Coeff) { s.insert(d); });
        return s;
    }

    Coeff eval(const std::vector<Coeff>& pt) const {
        if (pt.size() != nvars()) throw std::invalid_argument("wrong size");
        Coeff val = 0.0;
        traverse([&](const DegreeVector& d, Coeff c) {
            Coeff term = c;
            for (size_t i = 0; i < d.size(); ++i)
                if (d[i] != 0) term *= std::pow(pt[i], d[i]);
            val += term;
            });
        return val;
    }

    bool isHomog() const {
        bool first = true;
        int total = 0;
        bool hom = true;
        traverse([&](const DegreeVector& d, Coeff) {
            int deg = PolynomialRing::totalDeg(d);
            if (first) { total = deg; first = false; }
            else if (total != deg) hom = false;
            });
        return hom;
    }

    int homogDeg() const {
        bool first = true;
        int total = 0;
        bool hom = true;
        traverse([&](const DegreeVector& d, Coeff) {
            int deg = PolynomialRing::totalDeg(d);
            if (first) { total = deg; first = false; }
            else if (total != deg) hom = false;
            });
        return hom ? total : -1;
    }

    Polynomial homogPart(int d) const {
        Polynomial p(ring_);
        traverse([&](const DegreeVector& degs, Coeff c) {
            if (PolynomialRing::totalDeg(degs) == d) p.addMonomial(degs, c);
            });
        return p;
    }

    DegreeVector multideg() const { return leadDeg(); }

    Coeff lc() const {
        auto deg = leadDeg();
        Coeff c = 0.0;
        getCoeff(root_.get(), 0, deg, c);
        return c;
    }

    Polynomial lm() const {
        Polynomial m(ring_);
        m.addMonomial(leadDeg(), 1.0);
        return m;
    }

    Polynomial lt() const {
        Polynomial t(ring_);
        t.addMonomial(leadDeg(), lc());
        return t;
    }

    static Polynomial S_poly(const Polynomial& f, const Polynomial& g) {
        f.checkSameRing(g);
        Coeff lc_f = f.lc();
        Coeff lc_g = g.lc();
        if (std::abs(lc_f) < EPS || std::abs(lc_g) < EPS)
            throw std::runtime_error("Leading coefficient is zero in S-polynomial");

        auto df = f.multideg(), dg = g.multideg();
        auto l = PolynomialRing::lcm(df, dg);

        Polynomial m1(f.ring_);
        DegreeVector diff1 = l;
        for (size_t i = 0; i < diff1.size(); ++i) diff1[i] -= df[i];
        m1.addMonomial(diff1, 1.0 / lc_f);

        Polynomial m2(g.ring_);
        DegreeVector diff2 = l;
        for (size_t i = 0; i < diff2.size(); ++i) diff2[i] -= dg[i];
        m2.addMonomial(diff2, 1.0 / lc_g);

        return m1 * f - m2 * g;
    }

    static Polynomial reduce(const Polynomial& p, const std::vector<Polynomial>& basis) {
        if (basis.empty()) return p;
        auto ring = p.ring_;
        Polynomial rem(ring);
        Polynomial cur = p;

        while (!cur.isZero()) {
            auto lt_cur = cur.lt();
            if (lt_cur.isZero()) break;  // защита от зацикливания

            auto dcur = cur.multideg();
            bool divided = false;

            for (const auto& g : basis) {
                if (g.isZero()) continue;
                auto dg = g.multideg();
                if (PolynomialRing::divides(dcur, dg)) {
                    DegreeVector diff = dcur;
                    for (size_t i = 0; i < diff.size(); ++i) diff[i] -= dg[i];
                    Coeff factor = lt_cur.lc() / g.lc();
                    Polynomial sub(ring);
                    sub.addMonomial(diff, factor);
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

    static bool isGroebner(const std::vector<Polynomial>& basis) {
        if (basis.empty()) return true;
        for (size_t i = 0; i < basis.size(); ++i) {
            for (size_t j = i + 1; j < basis.size(); ++j) {
                Polynomial s = S_poly(basis[i], basis[j]);
                if (s.isZero()) continue;
                Polynomial rem = reduce(s, basis);
                if (!rem.isZero()) {
                    return false;
                }
            }
        }
        return true;
    }

    bool isZero() const {
        return isZeroNode(root_.get());
    }

    void traverse(std::function<void(const DegreeVector&, Coeff)> f) const {
        DegreeVector cur(nvars(), 0);
        trav(root_.get(), 0, cur, f);
    }

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
        if (p.isZero()) return os << "0";
        bool first = true;
        p.traverse([&](const DegreeVector& d, Coeff c) {
            if (std::abs(c) < Polynomial::EPS) return;
            if (!first) os << (c > 0 ? " + " : " - ");
            else if (c < 0) os << "-";
            double ac = std::abs(c);
            if (ac != 1.0 || PolynomialRing::totalDeg(d) == 0) os << ac;
            bool any = false;
            for (size_t i = 0; i < d.size(); ++i) {
                if (d[i] > 0) {
                    os << p.ring_->vars()[i];
                    if (d[i] != 1) os << "^" << d[i];
                    any = true;
                }
            }
            if (!any && ac == 1.0) os << "1";
            first = false;
            });
        return os;
    }

private:
    struct Node {
        Coeff coeff = 0.0;
        std::map<int, std::unique_ptr<Node>> children;
    };

    RingPtr ring_;
    std::unique_ptr<Node> root_;

    void checkSameRing(const Polynomial& o) const {
        if (ring_ != o.ring_) throw std::invalid_argument("different rings");
    }

    static std::unique_ptr<Node> cloneNode(const Node* src) {
        if (!src) return nullptr;
        auto n = std::make_unique<Node>();
        n->coeff = src->coeff;
        for (const auto& p : src->children) {
            int deg = p.first;
            const auto& ch = p.second;
            n->children[deg] = cloneNode(ch.get());
        }
        return n;
    }

    void insert(Node* node, size_t depth, const DegreeVector& degs, Coeff c) {
        if (depth == nvars()) { node->coeff += c; return; }
        int d = degs[depth];
        auto& child = node->children[d];
        if (!child) child = std::make_unique<Node>();
        insert(child.get(), depth + 1, degs, c);
    }

    static bool getCoeff(const Node* node, size_t depth, const DegreeVector& degs, Coeff& out) {
        if (!node) return false;
        if (depth == degs.size()) { out = node->coeff; return true; }
        auto it = node->children.find(degs[depth]);
        if (it == node->children.end()) return false;
        return getCoeff(it->second.get(), depth + 1, degs, out);
    }

    static void addTrees(Node* res, const Node* a, const Node* b, size_t depth) {
        if (!a && !b) return;
        if (a) res->coeff += a->coeff;
        if (b) res->coeff += b->coeff;
        std::set<int> degs;
        if (a) for (const auto& p : a->children) degs.insert(p.first);
        if (b) for (const auto& p : b->children) degs.insert(p.first);
        for (int d : degs) {
            const Node* ca = nullptr, * cb = nullptr;
            if (a) { auto it = a->children.find(d); if (it != a->children.end()) ca = it->second.get(); }
            if (b) { auto it = b->children.find(d); if (it != b->children.end()) cb = it->second.get(); }
            auto& rc = res->children[d];
            if (!rc) rc = std::make_unique<Node>();
            addTrees(rc.get(), ca, cb, depth + 1);
        }
    }

    static void subTrees(Node* res, const Node* a, const Node* b, size_t depth) {
        if (!a && !b) return;
        if (a) res->coeff += a->coeff;
        if (b) res->coeff -= b->coeff;
        std::set<int> degs;
        if (a) for (const auto& p : a->children) degs.insert(p.first);
        if (b) for (const auto& p : b->children) degs.insert(p.first);
        for (int d : degs) {
            const Node* ca = nullptr, * cb = nullptr;
            if (a) { auto it = a->children.find(d); if (it != a->children.end()) ca = it->second.get(); }
            if (b) { auto it = b->children.find(d); if (it != b->children.end()) cb = it->second.get(); }
            auto& rc = res->children[d];
            if (!rc) rc = std::make_unique<Node>();
            subTrees(rc.get(), ca, cb, depth + 1);
        }
    }

    static void mulScalar(Node* res, const Node* src, size_t depth, Coeff s) {
        if (!src) return;
        res->coeff = src->coeff * s;
        for (const auto& p : src->children) {
            int deg = p.first;
            const auto& ch = p.second;
            auto& rc = res->children[deg];
            if (!rc) rc = std::make_unique<Node>();
            mulScalar(rc.get(), ch.get(), depth + 1, s);
        }
    }

    static bool treesEq(const Node* a, const Node* b, size_t depth) {
        if (!a && !b) return true;
        if (!a || !b) return false;
        if (std::abs(a->coeff - b->coeff) > EPS) return false;
        if (a->children.size() != b->children.size()) return false;
        for (const auto& p : a->children) {
            int deg = p.first;
            const auto& ca = p.second;
            auto it = b->children.find(deg);
            if (it == b->children.end()) return false;
            if (!treesEq(ca.get(), it->second.get(), depth + 1)) return false;
        }
        return true;
    }

    bool isZeroNode(const Node* node) const {
        if (!node) return true;
        if (std::abs(node->coeff) > EPS) return false;
        for (const auto& p : node->children) {
            if (!isZeroNode(p.second.get())) return false;
        }
        return true;
    }

    void shift(const DegreeVector& s) {
        Polynomial sh(ring_);
        traverse([&](const DegreeVector& d, Coeff c) {
            DegreeVector nd = d;
            for (size_t i = 0; i < nd.size(); ++i) nd[i] += s[i];
            sh.addMonomial(nd, c);
            });
        *this = std::move(sh);
    }

    void trav(const Node* node, size_t depth, DegreeVector& cur,
        const std::function<void(const DegreeVector&, Coeff)>& f) const {
        if (!node) return;
        if (std::abs(node->coeff) > EPS) f(cur, node->coeff);
        if (depth == nvars()) return;
        for (const auto& p : node->children) {
            int deg = p.first;
            const auto& ch = p.second;
            cur[depth] = deg;
            trav(ch.get(), depth + 1, cur, f);
        }
        cur[depth] = 0;
    }

    DegreeVector leadDeg() const {
        DegreeVector best;
        bool first = true;
        traverse([&](const DegreeVector& d, Coeff c) {
            if (std::abs(c) < EPS) return;
            if (first) { best = d; first = false; }
            else if (ring_->cmp(d, best) > 0) best = d;
            });
        if (first) return DegreeVector(nvars(), 0);
        return best;
    }
};

std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

Polynomial parsePolynomial(const std::string& line, Polynomial::RingPtr ring) {
    Polynomial p(ring);
    std::string s = trim(line);
    if (s.empty() || s == "0") return p;

    size_t pos = 0;
    while (pos < s.length()) {
        double coeff = 1.0;
        int sign = 1;
        if (s[pos] == '+') { sign = 1; pos++; }
        else if (s[pos] == '-') { sign = -1; pos++; }

        size_t start = pos;
        while (pos < s.length() && s[pos] != '+' && s[pos] != '-') pos++;
        std::string term = trim(s.substr(start, pos - start));
        if (term.empty()) continue;

        size_t varPos = 0;
        while (varPos < term.length() && !isalpha(term[varPos])) varPos++;
        std::string coeffStr = term.substr(0, varPos);
        if (!coeffStr.empty()) {
            coeff = std::stod(coeffStr);
        }
        else {
            coeff = 1.0;
        }
        coeff *= sign;

        Polynomial::DegreeVector degs(ring->nvars(), 0);
        size_t i = varPos;
        while (i < term.length()) {
            while (i < term.length() && isspace(term[i])) i++;
            size_t nameStart = i;
            while (i < term.length() && isalpha(term[i])) i++;
            std::string varName = term.substr(nameStart, i - nameStart);
            int varIdx = -1;
            for (size_t j = 0; j < ring->vars().size(); j++) {
                if (ring->vars()[j] == varName) { varIdx = j; break; }
            }
            if (varIdx == -1) throw std::runtime_error("Unknown variable: " + varName);

            int exp = 1;
            if (i < term.length() && term[i] == '^') {
                i++;
                size_t numStart = i;
                while (i < term.length() && isdigit(term[i])) i++;
                exp = std::stoi(term.substr(numStart, i - numStart));
            }
            degs[varIdx] += exp;
        }
        p.addMonomial(degs, coeff);
    }
    return p;
}

int main() {
    std::cout << "=== Polynomial Manipulation Application ===\n";
    std::cout << "Enter number of variables: ";
    size_t n;
    std::cin >> n;
    std::cin.ignore();

    std::vector<std::string> varNames(n);
    std::cout << "Enter variable names (separated by spaces): ";
    for (size_t i = 0; i < n; ++i) std::cin >> varNames[i];
    std::cin.ignore();

    std::cout << "Choose monomial ordering:\n";
    std::cout << "0: lex, 1: grlex, 2: grevlex, 3: invlex, 4: rinvlex\n";
    int ordChoice;
    std::cin >> ordChoice;
    std::cin.ignore();
    OrderingType ord;
    switch (ordChoice) {
    case 0: ord = OrderingType::lex; break;
    case 1: ord = OrderingType::grlex; break;
    case 2: ord = OrderingType::grevlex; break;
    case 3: ord = OrderingType::invlex; break;
    case 4: ord = OrderingType::rinvlex; break;
    default: std::cerr << "Invalid choice, using grlex\n"; ord = OrderingType::grlex;
    }

    auto ring = std::make_shared<PolynomialRing>(varNames, ord);

    std::cout << "Enter polynomial f (e.g., 3x^2y - 2z + xyz):\n";
    std::string line;
    std::getline(std::cin, line);
    Polynomial f = parsePolynomial(line, ring);

    std::cout << "Enter polynomial g:\n";
    std::getline(std::cin, line);
    Polynomial g = parsePolynomial(line, ring);

    std::cout << "\n--- Results ---\n";
    std::cout << "f = " << f << "\n";
    std::cout << "g = " << g << "\n";
    std::cout << "f + g = " << f + g << "\n";
    std::cout << "f - g = " << f - g << "\n";
    std::cout << "f * g = " << f * g << "\n\n";

    std::cout << "supp(f) = { ";
    for (const auto& d : f.support()) {
        std::cout << "(";
        for (size_t i = 0; i < d.size(); ++i) {
            std::cout << d[i] << (i + 1 < d.size() ? "," : "");
        }
        std::cout << ") ";
    }
    std::cout << "}\n\n";

    Polynomial h = f;
    std::cout << "f == h : " << (f == h) << "\n";
    std::cout << "f == g : " << (f == g) << "\n\n";

    std::cout << "Enter evaluation point (space-separated values): ";
    std::vector<double> pt(n);
    for (size_t i = 0; i < n; ++i) std::cin >> pt[i];
    std::cin.ignore();
    std::cout << "f(";
    for (size_t i = 0; i < n; ++i) std::cout << pt[i] << (i + 1 < n ? "," : "");
    std::cout << ") = " << f.eval(pt) << "\n\n";

    std::cout << "Is f homogeneous? " << f.isHomog() << "\n";
    std::cout << "Homogeneity degree of f: " << f.homogDeg() << "\n";
    std::cout << "Enter degree for homogeneous part extraction: ";
    int degPart;
    std::cin >> degPart;
    std::cin.ignore();
    std::cout << "Homogeneous part of degree " << degPart << ": " << f.homogPart(degPart) << "\n\n";

    std::cout << "Leading monomial info for f:\n";
    std::cout << "  lm(f) = " << f.lm() << "\n";
    std::cout << "  lt(f) = " << f.lt() << "\n";
    std::cout << "  lc(f) = " << f.lc() << "\n";
    auto md_f = f.multideg();
    std::cout << "  multideg(f) = (";
    for (size_t i = 0; i < md_f.size(); ++i) std::cout << md_f[i] << (i + 1 < md_f.size() ? "," : "");
    std::cout << ")\n\n";

    std::cout << "Leading monomial info for g:\n";
    std::cout << "  lm(g) = " << g.lm() << "\n";
    std::cout << "  lt(g) = " << g.lt() << "\n";
    std::cout << "  lc(g) = " << g.lc() << "\n";
    auto md_g = g.multideg();
    std::cout << "  multideg(g) = (";
    for (size_t i = 0; i < md_g.size(); ++i) std::cout << md_g[i] << (i + 1 < md_g.size() ? "," : "");
    std::cout << ")\n\n";

    Polynomial S = Polynomial::S_poly(f, g);
    std::cout << "S(f,g) = " << S << "\n";
    auto md_S = S.multideg();
    std::cout << "multideg(f) = (";
    for (size_t i = 0; i < md_f.size(); ++i) std::cout << md_f[i] << (i + 1 < md_f.size() ? "," : "");
    std::cout << ")\n";
    std::cout << "multideg(g) = (";
    for (size_t i = 0; i < md_g.size(); ++i) std::cout << md_g[i] << (i + 1 < md_g.size() ? "," : "");
    std::cout << ")\n";
    std::cout << "multideg(S) = (";
    for (size_t i = 0; i < md_S.size(); ++i) std::cout << md_S[i] << (i + 1 < md_S.size() ? "," : "");
    std::cout << ")\n\n";

    std::cout << "Check if {f, g} is a Groebner basis... ";
    std::vector<Polynomial> basis = { f, g };
    bool gb = Polynomial::isGroebner(basis);
    std::cout << (gb ? "Yes" : "No") << "\n";
    if (!gb) {
        basis.push_back(S);
        std::cout << "After adding S-polynomial: " << (Polynomial::isGroebner(basis) ? "Yes" : "No") << "\n";
    }

    return 0;
}