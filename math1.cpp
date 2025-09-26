#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cctype>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <functional>

const double PI = std::acos(-1.0);
const double E = std::exp(1.0);

static inline std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

double toRadians(double deg) { return deg * PI / 180.0; }
double toDegrees(double rad) { return rad * 180.0 / PI; }

double parsePowerExpr(const std::string &raw) {
    std::string s = trim(raw);
    auto pos = s.find('^');
    if (pos == std::string::npos) {
        if (s == "pi") return PI;
        if (s == "e") return E;
        try { return std::stod(s); }
        catch (...) { throw std::runtime_error("Invalid number: '" + s + "'"); }
    } else {
        std::string left = s.substr(0, pos);
        std::string right = s.substr(pos + 1);
        double base = parsePowerExpr(left);
        double exponent = parsePowerExpr(right);
        return std::pow(base, exponent);
    }
}

double readDoubleWithPowers(const std::string &prompt) {
    while (true) {
        std::cout << prompt;
        std::string line;
        if (!std::getline(std::cin, line)) throw std::runtime_error("Input closed.");
        if (trim(line).empty()) continue;
        try { return parsePowerExpr(line); }
        catch (const std::exception &e) { std::cout << "Bad input. Try again. E.g. 2, -3.5, 2^3, 4^0.5, pi, e\n"; }
    }
}

static std::string fmtNum(double x) {
    std::ostringstream ss; ss << std::fixed << std::setprecision(6) << x; return ss.str();
}

void solveQuadratic(double a, double b, double c) {
    const double EPS = 1e-12;
    if (std::abs(a) < EPS) {
        if (std::abs(b) < EPS) {
            if (std::abs(c) < EPS) std::cout << "Infinite solutions.\n";
            else std::cout << "No solution.\n";
        } else {
            double x = -c / b;
            std::cout << "Linear. x = " << std::fixed << std::setprecision(6) << x << "\n";
        }
        return;
    }
    double disc = b*b - 4*a*c;
    if (disc > EPS) {
        double r1 = (-b + std::sqrt(disc)) / (2*a);
        double r2 = (-b - std::sqrt(disc)) / (2*a);
        std::cout << "Real distinct roots: x1 = " << fmtNum(r1) << ", x2 = " << fmtNum(r2) << "\n";
    } else if (std::abs(disc) <= EPS) {
        double r = -b / (2*a);
        std::cout << "Real repeated root: x = " << fmtNum(r) << "\n";
    } else {
        double re = -b / (2*a);
        double im = std::sqrt(-disc) / (2*a);
        std::cout << "Complex roots: x1 = " << fmtNum(re) << " + " << fmtNum(im) << "i, x2 = " << fmtNum(re) << " - " << fmtNum(im) << "i\n";
    }
}

double evalPolynomial(const std::vector<double> &coeffs, double x) {
    double res = 0.0;
    for (double c : coeffs) res = res * x + c;
    return res;
}

std::vector<double> derivativeCoeffs(const std::vector<double> &coeffs) {
    size_t nplus1 = coeffs.size();
    if (nplus1 <= 1) return {0.0};
    size_t n = nplus1 - 1;
    std::vector<double> d; d.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        double power = static_cast<double>(n - i);
        d.push_back(coeffs[i] * power);
    }
    return d;
}

std::vector<double> integralCoeffs(const std::vector<double> &coeffs) {
    std::vector<double> integral;
    integral.reserve(coeffs.size() + 1);
    
    for (size_t i = 0; i < coeffs.size(); ++i) {
        int power = static_cast<int>(coeffs.size() - 1 - i);
        double newCoeff = coeffs[i] / (power + 1);
        integral.push_back(newCoeff);
    }
    integral.push_back(0.0); // constant of integration
    
    return integral;
}

double numericalDerivative(std::function<double(double)> f, double x, double h = 1e-8) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

double trapezoidalRule(std::function<double(double)> f, double a, double b, int n = 1000) {
    double h = (b - a) / n;
    double sum = (f(a) + f(b)) / 2.0;
    
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }
    
    return sum * h;
}

double simpsonsRule(std::function<double(double)> f, double a, double b, int n = 1000) {
    if (n % 2 != 0) n++; // ensure even number
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; i += 2) {
        sum += 4.0 * f(a + i * h);
    }
    
    for (int i = 2; i < n; i += 2) {
        sum += 2.0 * f(a + i * h);
    }
    
    return sum * h / 3.0;
}

double newtonRaphson(std::function<double(double)> f, std::function<double(double)> df, 
                     double x0, double tolerance = 1e-10, int maxIter = 100) {
    double x = x0;
    for (int i = 0; i < maxIter; ++i) {
        double fx = f(x);
        double dfx = df(x);
        
        if (std::abs(dfx) < tolerance) {
            throw std::runtime_error("Derivative too small, cannot continue");
        }
        
        double newX = x - fx / dfx;
        
        if (std::abs(newX - x) < tolerance) {
            return newX;
        }
        
        x = newX;
    }
    
    throw std::runtime_error("Newton-Raphson failed to converge");
}

bool evaluateTrig(const std::string &func, double angle, bool angleIsDegrees, double &outVal) {
    double rad = angleIsDegrees ? toRadians(angle) : angle;
    std::string f = func;
    for (auto &ch : f) ch = std::tolower((unsigned char)ch);
    if (f == "sin") { outVal = std::sin(rad); return true; }
    if (f == "cos") { outVal = std::cos(rad); return true; }
    if (f == "tan") { outVal = std::tan(rad); return true; }
    if (f == "csc") { double s = std::sin(rad); if (std::abs(s) < 1e-15) return false; outVal = 1.0 / s; return true; }
    if (f == "sec") { double c = std::cos(rad); if (std::abs(c) < 1e-15) return false; outVal = 1.0 / c; return true; }
    if (f == "cot") { double t = std::tan(rad); if (std::abs(t) < 1e-15) return false; outVal = 1.0 / t; return true; }
    return false;
}

std::vector<double> inverseTrig(const std::string &func, double value) {
    std::string f = func;
    for (auto &ch : f) ch = std::tolower((unsigned char)ch);
    std::vector<double> sols;
    auto norm = [](double x){ while (x < 0) x += 360.0; while (x >= 360.0) x -= 360.0; return x; };
    if (f == "sin" || f == "arcsin") {
        if (value < -1.0 || value > 1.0) return sols;
        double r = std::asin(value);
        double d1 = toDegrees(r);
        double d2 = 180.0 - d1;
        double n1 = norm(d1), n2 = norm(d2);
        if (std::abs(n1 - n2) < 1e-9) sols.push_back(n1);
        else { sols.push_back(n1); sols.push_back(n2); }
        return sols;
    }
    if (f == "cos" || f == "arccos") {
        if (value < -1.0 || value > 1.0) return sols;
        double r = std::acos(value);
        double d1 = toDegrees(r);
        double d2 = 360.0 - d1;
        double n1 = norm(d1), n2 = norm(d2);
        if (std::abs(n1 - n2) < 1e-9) sols.push_back(n1);
        else { sols.push_back(n1); sols.push_back(n2); }
        return sols;
    }
    if (f == "tan" || f == "arctan") {
        double r = std::atan(value);
        double d1 = toDegrees(r);
        double n1 = norm(d1), n2 = norm(d1 + 180.0);
        sols.push_back(n1);
        if (std::abs(n1 - n2) > 1e-9) sols.push_back(n2);
        return sols;
    }
    return sols;
}

void solveRightAngledTriangle() {
    std::cout << "Right-angled triangle. All sides labelled relative to an acute angle θ.\n";
    std::cout << "Which knowns? 1) two sides  2) one side & angle θ (deg)\n";
    std::cout << "Enter 1 or 2: ";
    int ch;
    if (!(std::cin >> ch)) { std::cin.clear(); std::string d; std::getline(std::cin,d); std::cout<<"Bad\n"; return; }
    std::string rest; std::getline(std::cin, rest);

    if (ch == 1) {
        std::cout << "Type pair: 'opp adj', 'opp hyp' or 'adj hyp': ";
        std::string pair; std::getline(std::cin, pair); pair = trim(pair);
        for (auto &c : pair) c = std::tolower((unsigned char)c);

        if (pair == "opp adj" || pair == "adj opp") {
            double opp = readDoubleWithPowers("opposite: ");
            double adj = readDoubleWithPowers("adjacent: ");
            double hyp = std::sqrt(opp*opp + adj*adj);
            double theta = toDegrees(std::atan2(opp, adj));
            std::cout << "hyp = sqrt(opp^2 + adj^2), θ = arctan(opp/adj)\n";
            std::cout << "Hypotenuse = " << fmtNum(hyp) << ", θ = " << fmtNum(theta) << "°\n";
        } else if (pair == "opp hyp" || pair == "hyp opp") {
            double opp = readDoubleWithPowers("opposite: ");
            double hyp = readDoubleWithPowers("hypotenuse: ");
            if (hyp <= std::abs(opp)) { std::cout << "Invalid (hypotenuse must be longest).\n"; return; }
            double adj = std::sqrt(hyp*hyp - opp*opp);
            double theta = toDegrees(std::asin(opp / hyp));
            std::cout << "adj = sqrt(hyp^2 - opp^2), θ = arcsin(opp/hyp)\n";
            std::cout << "Adjacent = " << fmtNum(adj) << ", θ = " << fmtNum(theta) << "°\n";
        } else if (pair == "adj hyp" || pair == "hyp adj") {
            double adj = readDoubleWithPowers("adjacent: ");
            double hyp = readDoubleWithPowers("hypotenuse: ");
            if (hyp <= std::abs(adj)) { std::cout << "Invalid (hypotenuse must be longest).\n"; return; }
            double opp = std::sqrt(hyp*hyp - adj*adj);
            double theta = toDegrees(std::acos(adj / hyp));
            std::cout << "opp = sqrt(hyp^2 - adj^2), θ = arccos(adj/hyp)\n";
            std::cout << "Opposite = " << fmtNum(opp) << ", θ = " << fmtNum(theta) << "°\n";
        } else {
            std::cout << "Use 'opp adj', 'opp hyp' or 'adj hyp'.\n";
            return;
        }
    } else if (ch == 2) {
        std::cout << "Which side known? 'opposite', 'adjacent', or 'hypotenuse': ";
        std::string s; std::getline(std::cin, s); s = trim(s);
        for (auto &c : s) c = std::tolower((unsigned char)c);
        if (s != "opposite" && s != "adjacent" && s != "hypotenuse") { std::cout << "Bad label\n"; return; }
        double side = readDoubleWithPowers("side length: ");
        double theta = readDoubleWithPowers("θ (degrees): ");
        if (!(theta > 0.0 && theta < 90.0)) { std::cout << "θ must be acute.\n"; return; }
        double rad = toRadians(theta);
        if (s == "opposite") {
            double opp = side;
            double hyp = opp / std::sin(rad);
            double adj = opp / std::tan(rad);
            std::cout << "Hypotenuse = " << fmtNum(hyp) << ", Adjacent = " << fmtNum(adj) << "\n";
        } else if (s == "adjacent") {
            double adj = side;
            double hyp = adj / std::cos(rad);
            double opp = adj * std::tan(rad);
            std::cout << "Hypotenuse = " << fmtNum(hyp) << ", Opposite = " << fmtNum(opp) << "\n";
        } else {
            double hyp = side;
            double opp = hyp * std::sin(rad);
            double adj = hyp * std::cos(rad);
            std::cout << "Opposite = " << fmtNum(opp) << ", Adjacent = " << fmtNum(adj) << "\n";
        }
        std::cout << "Other acute angle = " << fmtNum(90.0 - theta) << "°\n";
    } else {
        std::cout << "Pick 1 or 2.\n";
    }
}

static void handleLinearIneq() {
    std::cout << "Solve ax + b [<, <=, >, >=] 0\n";
    double a = readDoubleWithPowers("a: ");
    double b = readDoubleWithPowers("b: ");
    std::string op;
    std::cout << "Operator (<, <=, >, >=): ";
    std::getline(std::cin, op); op = trim(op);
    if (op.empty()) { std::cout << "No operator.\n"; return; }

    const double EPS = 1e-12;
    if (std::abs(a) < EPS) {
        if (std::abs(b) < EPS) {
            std::cout << "All real numbers.\n";
        } else {
            bool ok=false;
            if (op == "<") ok = (b < 0);
            if (op == "<=") ok = (b <= 0);
            if (op == ">") ok = (b > 0);
            if (op == ">=") ok = (b >= 0);
            std::cout << (ok ? "All real numbers.\n" : "No solution.\n");
        }
        return;
    }

    double x0 = -b / a;
    if (op == "<") {
        if (a > 0) std::cout << "x < " << fmtNum(x0) << "\n";
        else std::cout << "x > " << fmtNum(x0) << "\n";
    } else if (op == "<=") {
        if (a > 0) std::cout << "x ≤ " << fmtNum(x0) << "\n";
        else std::cout << "x ≥ " << fmtNum(x0) << "\n";
    } else if (op == ">") {
        if (a > 0) std::cout << "x > " << fmtNum(x0) << "\n";
        else std::cout << "x < " << fmtNum(x0) << "\n";
    } else if (op == ">=") {
        if (a > 0) std::cout << "x ≥ " << fmtNum(x0) << "\n";
        else std::cout << "x ≤ " << fmtNum(x0) << "\n";
    } else {
        std::cout << "Unknown operator.\n";
    }
}

static std::string joinOr(const std::vector<std::string>& parts) {
    if (parts.empty()) return "";
    std::string s = parts[0];
    for (size_t i=1;i<parts.size();++i) s += " or " + parts[i];
    return s;
}

static void handleQuadraticIneq() {
    std::cout << "Solve ax^2 + bx + c [<, <=, >, >=] 0\n";
    double a = readDoubleWithPowers("a: ");
    double b = readDoubleWithPowers("b: ");
    double c = readDoubleWithPowers("c: ");
    std::string op;
    std::cout << "Operator (<, <=, >, >=): ";
    std::getline(std::cin, op); op = trim(op);
    if (op.empty()) { std::cout<<"No operator.\n"; return; }

    const double EPS = 1e-12;
    if (std::abs(a) < EPS) {
        if (std::abs(b) < EPS) {
            if (std::abs(c) < EPS) std::cout << "All real numbers.\n";
            else std::cout << "No solution.\n";
            return;
        }
        double x0 = -c / b;
        if (op == "<") { if (b>0) std::cout<<"x < "<<fmtNum(x0)<<"\n"; else std::cout<<"x > "<<fmtNum(x0)<<"\n"; }
        else if (op == "<=") { if (b>0) std::cout<<"x ≤ "<<fmtNum(x0)<<"\n"; else std::cout<<"x ≥ "<<fmtNum(x0)<<"\n"; }
        else if (op == ">") { if (b>0) std::cout<<"x > "<<fmtNum(x0)<<"\n"; else std::cout<<"x < "<<fmtNum(x0)<<"\n"; }
        else if (op == ">=") { if (b>0) std::cout<<"x ≥ "<<fmtNum(x0)<<"\n"; else std::cout<<"x ≤ "<<fmtNum(x0)<<"\n"; }
        else std::cout<<"Unknown operator.\n";
        return;
    }

    double disc = b*b - 4*a*c;
    if (disc < -EPS) {
        bool positiveEverywhere = (a > 0);
        if (op == ">" || op == ">=") {
            std::cout << (positiveEverywhere ? "All real numbers.\n" : "No solution.\n");
        } else {
            std::cout << (!positiveEverywhere ? "All real numbers.\n" : "No solution.\n");
        }
        return;
    }

    if (std::abs(disc) <= EPS) {
        double r = -b / (2*a);
        if (op == ">" || op == "<") {
            if ((op == ">" && a > 0) || (op == "<" && a < 0))
                std::cout << "x < " << fmtNum(r) << " or x > " << fmtNum(r) << "\n";
            else std::cout << "No solution.\n";
        } else if (op == ">=") {
            if (a > 0) std::cout << "All real numbers.\n"; else std::cout << "x = " << fmtNum(r) << "\n";
        } else if (op == "<=") {
            if (a < 0) std::cout << "All real numbers.\n"; else std::cout << "x = " << fmtNum(r) << "\n";
        } else std::cout << "Unknown operator.\n";
        return;
    }

    double sqrtD = std::sqrt(disc);
    double r1 = (-b - sqrtD) / (2*a);
    double r2 = (-b + sqrtD) / (2*a);
    if (r1 > r2) std::swap(r1, r2);

    std::vector<std::string> parts;
    if (op == ">") {
        if (a > 0) { parts.push_back("x < " + fmtNum(r1)); parts.push_back("x > " + fmtNum(r2)); }
        else { parts.push_back(fmtNum(r1) + " < x < " + fmtNum(r2)); }
    } else if (op == "<") {
        if (a > 0) { parts.push_back(fmtNum(r1) + " < x < " + fmtNum(r2)); }
        else { parts.push_back("x < " + fmtNum(r1)); parts.push_back("x > " + fmtNum(r2)); }
    } else if (op == ">=") {
        if (a > 0) { parts.push_back("x ≤ " + fmtNum(r1)); parts.push_back("x ≥ " + fmtNum(r2)); }
        else { parts.push_back(fmtNum(r1) + " ≤ x ≤ " + fmtNum(r2)); }
    } else if (op == "<=") {
        if (a > 0) { parts.push_back(fmtNum(r1) + " ≤ x ≤ " + fmtNum(r2)); }
        else { parts.push_back("x ≤ " + fmtNum(r1)); parts.push_back("x ≥ " + fmtNum(r2)); }
    } else {
        std::cout << "Unknown operator.\n"; return;
    }

    std::cout << "Solution: " << joinOr(parts) << "\n";
}

static void handleDomainRange() {
    std::cout << "Pick function: 1) Linear 2) Quadratic 3) Rational (linear/linear)\n";
    std::cout << "Choice: ";
    int ch; if (!(std::cin >> ch)) { std::cin.clear(); std::string d; std::getline(std::cin,d); std::cout<<"Bad\n"; return; }
    std::string rest; std::getline(std::cin, rest);

    if (ch == 1) {
        std::cout << "f(x) = ax + b\n";
        double a = readDoubleWithPowers("a: ");
        double b = readDoubleWithPowers("b: ");
        std::cout << "Domain: all real numbers\nRange: all real numbers\n";
    } else if (ch == 2) {
        std::cout << "f(x) = ax^2 + bx + c\n";
        double a = readDoubleWithPowers("a: ");
        double b = readDoubleWithPowers("b: ");
        double c = readDoubleWithPowers("c: ");
        double xv = -b / (2*a);
        double yv = a*xv*xv + b*xv + c;
        std::cout << "Domain: all real numbers\n";
        if (a > 0) std::cout << "Range: y ≥ " << fmtNum(yv) << " (min at x = " << fmtNum(xv) << ")\n";
        else std::cout << "Range: y ≤ " << fmtNum(yv) << " (max at x = " << fmtNum(xv) << ")\n";
    } else if (ch == 3) {
        std::cout << "f(x) = (a1 x + b1) / (a2 x + b2)\n";
        double a1 = readDoubleWithPowers("a1: ");
        double b1 = readDoubleWithPowers("b1: ");
        double a2 = readDoubleWithPowers("a2: ");
        double b2 = readDoubleWithPowers("b2: ");
        if (std::abs(a2) < 1e-12 && std::abs(b2) < 1e-12) { std::cout << "Denominator zero always -> no domain.\n"; return; }
        if (std::abs(a2) < 1e-12) {
            std::cout << "Domain: all real numbers\n";
            if (std::abs(a1) < 1e-12) std::cout << "Range: single value " << fmtNum(b1 / b2) << "\n";
            else std::cout << "Range: all real numbers\n";
        } else {
            double xr = -b2 / a2;
            std::cout << "Domain: all real numbers except x = " << fmtNum(xr) << "\n";
            double excluded = a1 / a2;
            std::cout << "Range: all real numbers except y = " << fmtNum(excluded) << "\n";
        }
    } else {
        std::cout << "Not supported.\n";
    }
}

void runInequalitiesAndDomainRange() {
    while (true) {
        std::cout << "\nInequalities & domain/range menu:\n";
        std::cout << "1) Linear inequality\n2) Quadratic inequality\n3) Domain & range (simple types)\n4) Back\nPick: ";
        int choice; if (!(std::cin >> choice)) { std::cin.clear(); std::string d; std::getline(std::cin,d); std::cout<<"Bad\n"; continue; }
        std::string rest; std::getline(std::cin, rest);
        if (choice == 1) handleLinearIneq();
        else if (choice == 2) handleQuadraticIneq();
        else if (choice == 3) handleDomainRange();
        else if (choice == 4) break;
        else std::cout << "Pick 1-4.\n";
    }
}

void runCalculusMenu() {
    while (true) {
        std::cout << "\nCalculus menu:\n";
        std::cout << "1) Polynomial integration\n2) Definite integral (numerical)\n3) Find roots (Newton-Raphson)\n4) Critical points (polynomial)\n5) Back\nPick: ";
        int choice; 
        if (!(std::cin >> choice)) { 
            std::cin.clear(); 
            std::string d; 
            std::getline(std::cin, d); 
            std::cout << "Bad\n"; 
            continue; 
        }
        std::string rest; 
        std::getline(std::cin, rest);

        try {
            if (choice == 1) {
                std::cout << "Polynomial integration: ∫P(x)dx\n";
                std::cout << "power (non-negative integer): ";
                int deg;
                if (!(std::cin >> deg) || deg < 0) { 
                    std::cin.clear(); 
                    std::getline(std::cin, rest); 
                    std::cout << "Bad power.\n"; 
                    continue; 
                }
                std::getline(std::cin, rest);
                
                std::vector<double> coeffs; 
                coeffs.reserve(deg + 1);
                std::cout << "Enter coeffs from highest power down to constant.\n";
                for (int i = deg; i >= 0; --i) {
                    std::ostringstream p; 
                    p << "coeff for x^" << i << ": ";
                    double c = readDoubleWithPowers(p.str());
                    coeffs.push_back(c);
                }
                
                auto integral = integralCoeffs(coeffs);
                
                std::cout << "Integral coefficients (+ C):\n";
                for (size_t i = 0; i < integral.size() - 1; ++i) {
                    int power = static_cast<int>(integral.size() - 2 - i);
                    std::cout << "coeff x^" << power << " = " << fmtNum(integral[i]) << "\n";
                }
                std::cout << "constant = C (arbitrary)\n";
                
            } else if (choice == 2) {
                std::cout << "Numerical integration of polynomial\n";
                std::cout << "Highest power: ";
                int deg;
                if (!(std::cin >> deg) || deg < 0) { 
                    std::cin.clear(); 
                    std::getline(std::cin, rest); 
                    std::cout << "Bad power.\n"; 
                    continue; 
                }
                std::getline(std::cin, rest);
                
                std::vector<double> coeffs; 
                coeffs.reserve(deg + 1);
                std::cout << "Enter coeffs from highest power down to constant.\n";
                for (int i = deg; i >= 0; --i) {
                    std::ostringstream p; 
                    p << "coeff for x^" << i << ": ";
                    double c = readDoubleWithPowers(p.str());
                    coeffs.push_back(c);
                }
                
                double a = readDoubleWithPowers("Lower bound a: ");
                double b = readDoubleWithPowers("Upper bound b: ");
                
                auto f = [&coeffs](double x) { return evalPolynomial(coeffs, x); };
                
                double trapResult = trapezoidalRule(f, a, b);
                double simpResult = simpsonsRule(f, a, b);
                
                std::cout << "Trapezoidal rule result: " << fmtNum(trapResult) << "\n";
                std::cout << "Simpson's rule result: " << fmtNum(simpResult) << "\n";
                
            } else if (choice == 3) {
                std::cout << "Find roots using Newton-Raphson method\n";
                std::cout << "Highest power: ";
                int deg;
                if (!(std::cin >> deg) || deg < 1) { 
                    std::cin.clear(); 
                    std::getline(std::cin, rest); 
                    std::cout << "Bad power (must be ≥ 1).\n"; 
                    continue; 
                }
                std::getline(std::cin, rest);
                
                std::vector<double> coeffs; 
                coeffs.reserve(deg + 1);
                std::cout << "Enter coeffs from highest power down to constant.\n";
                for (int i = deg; i >= 0; --i) {
                    std::ostringstream p; 
                    p << "coeff for x^" << i << ": ";
                    double c = readDoubleWithPowers(p.str());
                    coeffs.push_back(c);
                }
                
                auto derivCoeffs = derivativeCoeffs(coeffs);
                
                auto f = [&coeffs](double x) { return evalPolynomial(coeffs, x); };
                auto df = [&derivCoeffs](double x) { return evalPolynomial(derivCoeffs, x); };
                
                double x0 = readDoubleWithPowers("Initial guess x0: ");
                
                try {
                    double root = newtonRaphson(f, df, x0);
                    std::cout << "Root found: x = " << fmtNum(root) << "\n";
                    std::cout << "Verification: f(" << fmtNum(root) << ") = " << fmtNum(f(root)) << "\n";
                } catch (const std::exception& e) {
                    std::cout << "Error: " << e.what() << "\n";
                }
                
            } else if (choice == 4) {
                std::cout << "Find critical points (where derivative = 0)\n";
                std::cout << "Highest power (≥ 2): ";
                int deg;
                if (!(std::cin >> deg) || deg < 2) { 
                    std::cin.clear(); 
                    std::getline(std::cin, rest); 
                    std::cout << "Bad power (must be ≥ 2).\n"; 
                    continue; 
                }
                std::getline(std::cin, rest);
                
                std::vector<double> coeffs; 
                coeffs.reserve(deg + 1);
                std::cout << "Enter coeffs from highest power down to constant.\n";
                for (int i = deg; i >= 0; --i) {
                    std::ostringstream p; 
                    p << "coeff for x^" << i << ": ";
                    double c = readDoubleWithPowers(p.str());
                    coeffs.push_back(c);
                }
                
                auto derivCoeffs = derivativeCoeffs(coeffs);
                auto secondDerivCoeffs = derivativeCoeffs(derivCoeffs);
                
                std::cout << "First derivative coefficients:\n";
                for (size_t i = 0; i < derivCoeffs.size(); ++i) {
                    int power = static_cast<int>(derivCoeffs.size() - 1 - i);
                    std::cout << "coeff x^" << power << " = " << fmtNum(derivCoeffs[i]) << "\n";
                }
                
                if (deg == 2) {
                    double a = derivCoeffs[0];
                    double b = derivCoeffs[1];
                    if (std::abs(a) > 1e-12) {
                        double criticalPoint = -b / a;
                        auto f = [&coeffs](double x) { return evalPolynomial(coeffs, x); };
                        auto f2 = [&secondDerivCoeffs](double x) { return evalPolynomial(secondDerivCoeffs, x); };
                        
                        double fVal = f(criticalPoint);
                        double f2Val = f2(criticalPoint);
                        
                        std::cout << "Critical point: x = " << fmtNum(criticalPoint) << "\n";
                        std::cout << "Function value: f(" << fmtNum(criticalPoint) << ") = " << fmtNum(fVal) << "\n";
                        
                        if (f2Val > 0) std::cout << "This is a local minimum.\n";
                        else if (f2Val < 0) std::cout << "This is a local maximum.\n";
                        else std::cout << "Second derivative test inconclusive.\n";
                    } else {
                        std::cout << "No critical points (constant derivative).\n";
                    }
                } else {
                    std::cout << "For higher power polynomials, use Newton-Raphson on the derivative.\n";
                    std::cout << "Try different initial guesses to find multiple critical points.\n";
                }
                
            } else if (choice == 5) {
                break;
            } else {
                std::cout << "Pick 1-5.\n";
            }
        } catch (const std::exception& e) {
            std::cout << "Error: " << e.what() << "\n";
        }
    }
}

int main() {
    std::cout << "math stuff advanced calc \n";
    std::cout << "You can use '^' in numeric input (e.g. 2^3, 4^0.5), plus 'pi' and 'e'.\nAngles default to degrees when prompted.\n";

    while (true) {
        std::cout << "\nChoose:\n";
        std::cout << "1) Solve quadratic\n2) Evaluate polynomial\n3) Show derivative of polynomial\n4) Compute a^b\n5) Trigonometry\n6) Right-angled triangle solver\n7) Inequalities & domain/range\n8) Calculus\n9) Exit\n";
        std::cout << "Enter 1-9: ";
        int choice = 0;
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::string dummy; std::getline(std::cin, dummy);
            std::cout << "Bad choice. Try again.\n";
            continue;
        }
        std::string rest; std::getline(std::cin, rest);

        try {
            if (choice == 1) {
                double a = readDoubleWithPowers("a (x^2 coeff): ");
                double b = readDoubleWithPowers("b (x coeff): ");
                double c = readDoubleWithPowers("c (const): ");
                solveQuadratic(a, b, c);
            } else if (choice == 2) {
                std::cout << "power (non-negative integer): ";
                int deg;
                if (!(std::cin >> deg) || deg < 0) { std::cin.clear(); std::getline(std::cin, rest); std::cout << "Bad power.\n"; continue; }
                std::getline(std::cin, rest);
                std::vector<double> coeffs; coeffs.reserve(deg+1);
                std::cout << "Enter coeffs from highest power down to constant.\n";
                for (int i = deg; i >= 0; --i) {
                    std::ostringstream p; p << "coeff for x^" << i << ": ";
                    double c = readDoubleWithPowers(p.str());
                    coeffs.push_back(c);
                }
                std::cout << "Enter x (may use ^): ";
                std::string xs; std::getline(std::cin, xs);
                double x = parsePowerExpr(xs);
                double val = evalPolynomial(coeffs, x);
                std::cout << "P(" << x << ") = " << std::fixed << std::setprecision(6) << val << "\n";
            } else if (choice == 3) {
                std::cout << "power (non-negative integer): ";
                int deg;
                if (!(std::cin >> deg) || deg < 0) { std::cin.clear(); std::getline(std::cin, rest); std::cout << "Bad power.\n"; continue; }
                std::getline(std::cin, rest);
                std::vector<double> coeffs; coeffs.reserve(deg+1);
                std::cout << "Enter coeffs from highest power down to constant.\n";
                for (int i = deg; i >= 0; --i) {
                    std::ostringstream p; p << "coeff for x^" << i << ": ";
                    double c = readDoubleWithPowers(p.str());
                    coeffs.push_back(c);
                }
                auto d = derivativeCoeffs(coeffs);
                int ddeg = static_cast<int>(d.size()) - 1;
                std::cout << "Derivative power: " << (ddeg >= 0 ? ddeg : 0) << "\n";
                for (size_t i = 0; i < d.size(); ++i) {
                    int power = static_cast<int>(d.size() - 1 - i);
                    std::cout << "coeff x^" << power << " = " << std::fixed << std::setprecision(6) << d[i] << "\n";
                }
            } else if (choice == 4) {
                double base = readDoubleWithPowers("Base: ");
                double exponent = readDoubleWithPowers("Exponent: ");
                double res = std::pow(base, exponent);
                std::cout << base << "^" << exponent << " = " << std::fixed << std::setprecision(8) << res << "\n";
            } else if (choice == 5) {
                std::cout << "\nTrig options:\n1) Evaluate trig (sin/cos/tan/csc/sec/cot)\n2) Inverse trig (arcsin/arccos/arctan)\nEnter 1 or 2: ";
                int tchoice;
                if (!(std::cin >> tchoice)) { std::cin.clear(); std::getline(std::cin, rest); std::cout << "Bad\n"; continue; }
                std::getline(std::cin, rest);
                if (tchoice == 1) {
                    std::cout << "Function name: ";
                    std::string f; std::getline(std::cin, f); f = trim(f);
                    std::cout << "Angle (may use ^): ";
                    std::string angs; std::getline(std::cin, angs);
                    double ang = parsePowerExpr(angs);
                    std::cout << "Angle in degrees? (y/n) [default y]: ";
                    std::string yn; std::getline(std::cin, yn);
                    bool deg = true;
                    if (!trim(yn).empty()) { char c = std::tolower(yn[0]); deg = (c == 'y'); }
                    double out;
                    bool ok = evaluateTrig(f, ang, deg, out);
                    if (!ok) std::cout << "Undefined or unknown function.\n";
                    else std::cout << f << "(" << ang << (deg ? "°" : " rad") << ") = " << std::fixed << std::setprecision(8) << out << "\n";
                } else if (tchoice == 2) {
                    std::cout << "Inverse (arcsin/arccos/arctan) or sin/cos/tan: ";
                    std::string f; std::getline(std::cin, f); f = trim(f);
                    std::cout << "Value: ";
                    double val = readDoubleWithPowers("");
                    auto sols = inverseTrig(f, val);
                    if (sols.empty()) std::cout << "No real solutions (out of domain or unknown function).\n";
                    else {
                        std::cout << "Solutions in degrees (0 ≤ θ < 360):\n";
                        for (double s : sols) std::cout << std::fixed << std::setprecision(6) << s << "°\n";
                    }
                } else {
                    std::cout << "Bad option.\n";
                }
            } else if (choice == 6) {
                solveRightAngledTriangle();
            } else if (choice == 7) {
                runInequalitiesAndDomainRange();
            } else if (choice == 8) {
                runCalculusMenu();
            } else if (choice == 9) {
                std::cout << "Bye.\n";
                break;
            } else {
                std::cout << "Pick 1-9.\n";
            }
        } catch (const std::exception &e) {
            std::cout << "Error: " << e.what() << "\n";
        }
    }

    return 0;
}