#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>



class BigInteger {
    friend bool operator<(const BigInteger& a, const BigInteger& b);
    friend bool operator==(const BigInteger& a, const BigInteger& b);

private:
    bool is_negative;
    std::vector<short int> digits;
    static const short int base = 100;

    void remove_leading_zeroes() {
        int i = digits.size() - 1;
        while (i > 0 && digits[i] == 0) {
            digits.pop_back();
            --i;
        }
    }

    void fft(std::vector<std::complex<double>>& a, bool invert) {
        int n = a.size();
        for (int i = 1, j = 0; i < n; ++i) {
            int bit = n >> 1;
            for (; j >= bit; bit >>= 1) j -= bit;
            j += bit;
            if (i < j) std::swap(a[i], a[j]);
        }

        for (int len = 2; len <= n; len <<= 1) {
            double ang = 2 * acos(-1) / len * (invert ? -1 : 1);
            std::complex<double> wlen(cos(ang), sin(ang));
            for (int i = 0; i < n; i += len) {
                std::complex<double> w(1);
                for (int j = 0; j < len / 2; ++j) {
                    std::complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert)
            for (int i = 0; i < n; ++i) a[i] /= n;
    }


public:
    BigInteger() : is_negative(false), digits(std::vector<short int>(1, 0)) {};

    BigInteger(const long long a) {
        is_negative = a < 0;
        if (a == 0) {
            digits.push_back(0);
            return;
        }
        long long a_copy;
        a_copy = is_negative ? -a : a;
        while (a_copy / base > 0 || a_copy > 0) {
            digits.push_back(a_copy % base);
            a_copy /= base;
        }
    }

    explicit BigInteger(const std::string& s) {
        is_negative = (s[0] == '-');
        int n = s.length() - 1;
        for (; n >= is_negative + 1; n -= 2) {
            digits.push_back((s[n] - '0') +
                10 * (s[n - 1] - '0'));
        }
        if (n >= static_cast<int>(is_negative)) digits.push_back(s[n] - '0');
    }

    explicit operator bool() const {
        return !(digits.size() == 1 && digits[0] == 0);
    }

    std::string toString() const {
        std::string res = "";
        if (is_negative) res += "-";
        int n = digits.size();
        for (int i = n - 1; i >= 0; --i) {
            std::string p = std::to_string(digits[i]);
            if (p.size() == 1 && i != n - 1) res += '0';
            res += p;
        }
        return res;
    }

    BigInteger abs() const {
        if (is_negative) return -*this;
        return *this;
    }

    BigInteger operator-() const {
        BigInteger t = *this;
        if (!(digits.size() == 1 && digits[0] == 0)) t.is_negative = !t.is_negative;
        return t;
    }

    BigInteger& operator+=(const BigInteger& a) {
        if (is_negative != a.is_negative) {
            if (is_negative) {
                is_negative = false;
                *this -= a;
                is_negative = !is_negative;
                return *this;
            }
            return *this -= -a;
        }
        if (a.digits.size() > digits.size()) digits.resize(a.digits.size());
        size_t n = a.digits.size();
        short int carry = 0;
        size_t i = 0;
        for (; i < n; ++i) {
            digits[i] += (a.digits[i] + carry);
            carry = digits[i] / base;
            digits[i] %= base;
        }
        while (carry > 0) {
            if (digits.size() == i) digits.push_back(0);
            digits[i] += carry;
            carry = digits[i] / base;
            digits[i] %= base;
            ++i;
        }
        return *this;
    }

    BigInteger& operator-=(const BigInteger& a) {
        if (*this == a) return *this = BigInteger();
        if (is_negative != a.is_negative) {
            if (is_negative) {
                is_negative = false;
                *this += a;
                is_negative = true;
                return *this;
            }
            return *this += -a;
        }
        if (abs() < a.abs()) {
            bool f = !is_negative;
            BigInteger t = abs();
            *this = a.abs();
            *this -= t;
            if (f) is_negative = !is_negative;
            return *this;
        }
        size_t n = a.digits.size();
        short int carry = 0;
        size_t i = 0;
        for (; i < n; ++i) {
            digits[i] -= (a.digits[i] + carry);
            carry = digits[i] < 0 ? 1 : 0;
            digits[i] = (digits[i] + base) % base;
        }
        while (carry > 0) {
            if (digits.size() == i) digits.push_back(0);
            digits[i] -= carry;
            carry = digits[i] < 0 ? 1 : 0;
            digits[i] = (digits[i] + base) % base;
            ++i;
        }
        remove_leading_zeroes();
        return *this;
    }

    BigInteger& operator*=(const BigInteger& a) {
        if (*this == 0 || a == 0) {
            *this = 0;
            return *this;
        }
        is_negative ^= a.is_negative;
        if (a.digits.size() == 1) {
            int n = digits.size();
            int carry = 0;
            int i = 0;
            for (; i < n; ++i) {
                digits[i] *= a.digits[0];
                digits[i] += carry;
                carry = digits[i] / base;
                digits[i] %= base;
            }
            while (carry > 0) {
                digits.push_back(0);
                digits[i] += carry;
                carry = digits[i] / base;
                digits[i] %= base;
                ++i;
            }
            return *this;
        }
        if (a.abs() == base) {
            digits.insert(digits.begin(), 1, 0);
            return *this;
        }
        std::vector<std::complex<double>> ft(digits.begin(), digits.end());
        std::vector<std::complex<double>> fa(a.digits.begin(), a.digits.end());
        size_t n = 1;
        while (n < std::max(digits.size(), a.digits.size())) n <<= 1;
        n <<= 1;
        fa.resize(n);
        ft.resize(n);
        fft(fa, false);
        fft(ft, false);
        for(size_t i = 0; i < n; ++i) ft[i] *= fa[i];
        fft(ft, true);
        digits.resize(n);
        int carry = 0;
        for (size_t i = 0; i < n; ++i) {
            int t = static_cast<int>(ft[i].real() + 0.5);
            t += carry;
            carry = t / base;
            digits[i] = t % base;
        }
        remove_leading_zeroes();
        return *this;
    }

    BigInteger& operator/=(const BigInteger& a) {
        if (*this == 0) return *this;
        if (this->abs() < a.abs()) {
            *this = 0;
            return *this;
        }
        is_negative ^= a.is_negative;
        int n = digits.size();
        BigInteger cur = 0;
        std::vector<int> ans;
        for (int i = n - 1; i >= 0; i--) {
            cur *= base;
            cur += digits[i];
            int l = 0;
            int r = base;
            BigInteger t;
            BigInteger d;
            if (cur < a.abs()) {
                ans.push_back(0);
                continue;
            }
            while (r - l > 1) {
                int m = (l + r) / 2;
                t = a.abs();
                t *= m;
                if (cur < t)
                    r = m;
                else {
                    l = m;
                    d = t;
                }
            }
            cur -= d;
            ans.push_back(l);
        }
        int s = ans.size();
        digits.resize(s);
        for (int i = 0; i < s; ++i) digits[i] = ans[s - 1 - i];
        remove_leading_zeroes();
        return *this;
    }

    BigInteger& operator%=(const BigInteger& a) {
        BigInteger t = *this;
        return *this -= ((t /= a) *= a);
    }

    BigInteger operator++(int) {
        BigInteger copy = *this;
        *this += 1;
        return copy;
    }

    BigInteger& operator++() { return *this += 1; }

    BigInteger operator--(int) {
        BigInteger copy = *this;
        *this -= 1;
        return copy;
    }

    BigInteger& operator--() { return *this -= 1; }
};

BigInteger operator+(const BigInteger& a, const BigInteger& b) {
    BigInteger a_copy = a;
    return a_copy += b;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b) {
    BigInteger a_copy = a;
    return a_copy -= b;
}

BigInteger operator*(const BigInteger& a, const BigInteger& b) {
    BigInteger a_copy = a;
    a_copy *= b;
    return a_copy;
}

BigInteger operator/(const BigInteger& a, const BigInteger& b) {
    BigInteger a_copy = a;
    a_copy /= b;
    return a_copy;
}

BigInteger operator%(const BigInteger& a, const BigInteger& b) {
    BigInteger a_copy = a;
    a_copy %= b;
    return a_copy;
}

bool operator<(const BigInteger& a, const BigInteger& b) {
    if (a.is_negative != b.is_negative) return a.is_negative > b.is_negative;
    if (a.digits.size() != b.digits.size()) return a.is_negative == (a.digits.size() > b.digits.size());
    int n = b.digits.size();
    for (int i = 1; i <= n; ++i) 
        if (a.digits[n - i] != b.digits[n - i]) return a.is_negative == (a.digits[n - i] > b.digits[n - i]);
    return false;
}

bool operator>(const BigInteger& a, const BigInteger& b) { return b < a; }

bool operator==(const BigInteger& a, const BigInteger& b) {
    if (a.is_negative != b.is_negative || a.digits.size() != b.digits.size())
        return false;
    for (size_t i = 0; i < a.digits.size(); ++i)
        if (a.digits[i] != b.digits[i]) return false;
    return true;
}

bool operator<=(const BigInteger& a, const BigInteger& b) {
    return !(b < a);
}

bool operator>=(const BigInteger& a, const BigInteger& b) {
    return !(a < b);
}

bool operator!=(const BigInteger& a, const BigInteger& b) { return !(a == b); }

std::ostream& operator<<(std::ostream& out, const BigInteger& a) {
    out << a.toString();
    return out;
}

std::istream& operator>>(std::istream& in, BigInteger& a) {
    std::string s;
    in >> s;
    a = BigInteger(s);
    return in;
}

class Rational {
    friend bool operator<(const Rational& a, const Rational& b);
    friend bool operator==(const Rational& a, const Rational& b);

private:
    BigInteger numerator;
    BigInteger denominator;

    BigInteger gcd(const BigInteger& a, const BigInteger& b) {
        return b ? gcd(b, a % b) : a;
    }

    void reduce() {
        BigInteger GCD = gcd(numerator >= 0 ? numerator : -numerator, denominator);
        numerator /= GCD;
        denominator /= GCD;
    }

public:
    Rational() : numerator(0), denominator(1) {};
    Rational(BigInteger n) : numerator(n), denominator(1) {};
    Rational(long long n) : numerator(n), denominator(1) {};

    Rational& operator+=(const Rational& r) {
        numerator *= r.denominator;
        numerator += (r.numerator * denominator);
        denominator *= r.denominator;
        reduce();
        return *this;
    }

    Rational& operator-=(const Rational& r) {
        numerator *= r.denominator;
        numerator -= (r.numerator * denominator);
        denominator *= r.denominator;
        reduce();
        return *this;
    }

    Rational& operator*=(const Rational& r) {
        numerator *= r.numerator;
        denominator *= r.denominator;
        reduce();
        return *this;
    }

    Rational& operator/=(const Rational& r) {
        if (r.numerator > 0) {
            numerator *= r.denominator;
            denominator *= r.numerator;
        }
        else {
            numerator *= -r.denominator;
            denominator *= -r.numerator;
        }
        reduce();
        return *this;
    }

    Rational operator-() const {
        Rational t = *this;
        t.numerator = -t.numerator;
        return t;
    }

    std::string toString() const {
        if (denominator == 1) return numerator.toString();
        return numerator.toString() + "/" + denominator.toString();
    }

    std::string asDecimal(size_t precision = 0) const {
        BigInteger deg = 1;
        for (size_t i = 0; i < precision; ++i) deg *= 10;
        std::string ans = "";
        if (numerator < 0) ans += '-';
        BigInteger t = numerator >= 0 ? numerator : -numerator;
        BigInteger c = t / denominator;
        ans += c.toString();
        if (precision == 0) return ans;
        t -= (c * denominator);
        t *= deg;
        t /= denominator;
        std::string d = t.toString();
        std::string dd = "";
        for (size_t i = 0; i < precision - d.size(); ++i) dd += '0';
        dd += d;
        ans += "." + dd;
        return ans;
    }

    explicit operator double() const {
        double res;
        std::stringstream(asDecimal(16)) >> res;
        return res;
    }
};

bool operator<(const Rational& a, const Rational& b) {
    return a.numerator * b.denominator < b.numerator* a.denominator;
}

bool operator==(const Rational& a, const Rational& b) {
    return a.numerator == b.numerator && a.denominator == b.denominator;
}

bool operator>(const Rational& a, const Rational& b) { return b < a; }

bool operator<=(const Rational& a, const Rational& b) {
    return (a < b) || (a == b);
}

bool operator>=(const Rational& a, const Rational& b) {
    return (a > b) || (a == b);
}

bool operator!=(const Rational& a, const Rational& b) { return !(a == b); }

Rational operator+(const Rational& a, const Rational& b) {
    Rational a_copy = a;
    return a_copy += b;
}

Rational operator-(const Rational& a, const Rational& b) {
    Rational a_copy = a;
    return a_copy -= b;
}

Rational operator*(const Rational& a, const Rational& b) {
    Rational a_copy = a;
    a_copy *= b;
    return a_copy;
}

Rational operator/(const Rational& a, const Rational& b) {
    Rational a_copy = a;
    a_copy /= b;
    return a_copy;
}