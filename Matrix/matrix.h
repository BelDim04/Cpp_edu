#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>



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
        for (; n >= is_negative + 3; n -= 4) {
            int cur = 0;
            int p = 1;
            for(int cn=n; cn>n-4; --cn){
                cur += (s[cn] - '0') * p;
                p *= 10;
            }
            digits.push_back(cur);
        }
        int cur = 0;
        int p = 1;
        while (n >= static_cast<int>(is_negative)) {
            cur += (s[n] - '0') * p;
            p *= 10;
            --n;
        }
        digits.push_back(cur);
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
        if (a.digits.size() == 1) {
            int n = digits.size();
            int carry = 0;
            int i = n - 1;
            for (; i >= 0; --i) {
                digits[i] += carry * base;
                carry = digits[i] % a.digits[0];
                digits[i] /= a.digits[0];
            }
            remove_leading_zeroes();
            return *this;
        }
        if (a.abs() == base) {
            digits.erase(digits.begin());
            return *this;
        }
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
        if(a == 2){
            *this = digits[0] % 2;
            return *this;
        }
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
        BigInteger ca = a;
        BigInteger cb = b;
        BigInteger ans = 1;
        while(true) {
            if(ca == 0){
                ans *= cb;
                return ans;
            }
            if(cb == 0){
                ans *= ca;
                return ans;
            }
            if ((ca % 2 == 0) && (cb % 2 == 0)) {
                ca /= 2;
                cb /= 2;
                ans *= 2;
                continue;
            }
            if ((ca % 2 == 0)) {
                ca /= 2;
                continue;
            }
            if ((cb % 2 == 0)) {
                cb /= 2;
                continue;
            }
            if(ca >= cb) ca -= cb;
            else cb -= ca;
        }
        return ans;
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

std::ostream& operator<<(std::ostream& out, const Rational& x) {
    out << x.toString();
    return out;
}

std::istream& operator>>(std::istream& in, Rational& x) {
    BigInteger y;
    in >> y;
    x = Rational(y);
    return in;
}

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

template<bool B, size_t K, size_t N>
struct IsPrimeHelper{
    static const bool value = N % K != 0 && IsPrimeHelper< (K*K > N), K + 1, N>::value;
};

template<size_t K, size_t N>
struct IsPrimeHelper<true, K, N>{
    static const bool value = true;
};

template<size_t N>
struct IsPrime{
    static const bool value = IsPrimeHelper<false, 2, N>::value;
};

template<size_t N>
class Residue{
    template<size_t M>
    friend bool operator ==(const Residue<M>& a, const Residue<M>& b);
private:
    size_t val;

    size_t deg(size_t num, size_t d){
        if(d == 1) return num;
        return d%2 == 0? (deg(num, d/2) * deg(num, d/2))%N: (num * deg(num, d-1))%N;
    }
public:
    explicit Residue(int k): val(k >= 0? k%(long long)N : ((long long)N - (-k) % (long long)N) % (long long)N) {}

    Residue(): val(0){}

    Residue& operator =(const Residue& other){
        val = other.val;
        return *this;
    }

    Residue& operator+=(const Residue& other){
        val = (val + other.val) % N;
        return  *this;
    }

    Residue& operator-=(const Residue& other){
        val = (val - other.val + (long long)N) % (long long)N ;
        return *this;
    }

    Residue& operator*=(const Residue& other){
        val = (val * other.val) % N;
        return  *this;
    }

    Residue& operator/=(const Residue& other){
        static_assert(IsPrime<N>::value, "Prime");
        *this *= Residue(deg(other.val, N-2));
        return *this;
    }

    explicit operator int() const{
        return val;
    }

};

template<size_t N>
Residue<N> operator +(const Residue<N>& a, const Residue<N>& b){
    Residue c = a;
    c+=b;
    return c;
}

template<size_t N>
Residue<N> operator -(const Residue<N>& a, const Residue<N>& b){
    Residue c = a;
    c-=b;
    return c;
}

template<size_t N>
Residue<N> operator *(const Residue<N>& a, const Residue<N>& b){
    Residue c = a;
    c*=b;
    return c;
}

template<size_t N>
Residue<N> operator /(const Residue<N>& a, const Residue<N>& b){
    Residue c = a;
    c/=b;
    return c;
}

template<size_t N>
bool operator ==(const Residue<N>& a, const Residue<N>& b){
    return a.val == b.val;
}

template<size_t N>
bool operator !=(const Residue<N>& a, const Residue<N>& b){
    return ! (a == b);
}


template<size_t Mod>
std::ostream& operator<<(std::ostream& out, const Residue<Mod>& x) {
    long long y = static_cast<int>(x);
    out << y;
    return out;
}

template<size_t Mod>
std::istream& operator>>(std::istream& in, Residue<Mod>& x) {
    long long y;
    in >> y;
    x = Residue<Mod>(y);
    return in;
}

template<size_t N, size_t M, typename Field = Rational>
class Matrix{
private:
    Field body[N][M];

    template<size_t K, typename F>
    void check_square(const Matrix<K,K,F>&) const{};
public:

    Matrix() {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] = Field(0);
            }
        }
    }

    explicit Matrix(const std::vector<std::vector<Field>>& v) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] = v[i][j];
            }
        }
    }

    Matrix(const std::vector<std::vector<long long>>& A) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] = Field(A[i][j]);
            }
        }
    }


    Matrix(const std::initializer_list<std::vector<long long>>& A) {
        size_t i = 0;
        for (auto row : A) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] = Field(row[j]);
            }
            ++i;
        }
    }

    Matrix(const std::initializer_list<long long>& A) {
        size_t i = 0;
        size_t j = 0;
        for (auto id : A) {
            body[i][j] = Field(id);
            ++j;
            i += j / M;
            j %= M;
        }
    }

    Matrix(const Matrix& other) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] = other.body[i][j];
            }
        }
    }

    Matrix& operator=(const Matrix& other) {
        if (&other == this) {
            return *this;
        }
        Matrix copyOfOther = other;
        std::swap(body, copyOfOther.body);
        return *this;
    }

    Matrix& operator +=(const Matrix& other){
        for(size_t i=0; i<N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] += other.body[i][j];
            }
        }
        return * this;
    }

    Matrix& operator-=(const Matrix& other) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] -= other.body[i][j];
            }
        }
        return *this;
    }

    Matrix& operator *=(const Field& num){
        for(size_t i=0; i<N; ++i)
            for(size_t j=0;j<M;++j)
                body[i][j] *= num;
        return * this;
    }

    Matrix<N, N, Field>& operator*=(const Matrix<N, N, Field>& other) {
        check_square<N,Field>(*this);
        return *this = *this * other;
    }

    Matrix<M,N,Field> transposed() const{
        Matrix<M,N,Field> tr;
        for(size_t i=0; i<M; ++i)
            for(size_t j=0; j<N; ++j)
                tr[i][j] = body[j][i];
        return tr;
    }

    Field trace() const{
        check_square<N,Field>(*this);
        Field ans(0);
        for(size_t i=0; i<N; ++i) ans += body[i][i];
        return ans;
    }

    std::vector<Field> getRow(size_t k) const{
        return std::vector(body[k], body[k] + M);
    }

    std::vector<Field> getColumn(size_t k) const{
        std::vector<Field> ans;
        ans.reserve(N);
        for(size_t i=0; i<M; ++i) ans.push_back(body[i][k]);
        return ans;
    }

    Field* operator [](size_t n){
        return body[n];
    }

    const Field* operator [](size_t n) const{
        return body[n];
    }

    std::pair<Matrix, size_t> triangleByGauss() const{
        Matrix ans = *this;
        size_t row = 0;
        size_t column = 0;
        size_t cnt = 0;
        while (row < N && column < M) {
            for (size_t i = row + 1; i < N; ++i) {
                if (ans.body[i][column] == Field(0)) {
                    continue;
                }
                if (ans.body[row][column] == Field(0)) {
                    for (size_t k = column; k < M; ++k) {
                        std::swap(ans.body[i][k], ans.body[row][k]);
                    }
                    ++cnt;
                    continue;
                }
                Field cur = ans.body[i][column] / ans.body[row][column];
                for (size_t k = column; k < M; ++k) {
                    ans.body[i][k] -= cur * ans.body[row][k];
                }
            }
            ++row;
            ++column;
        }
        return {ans, cnt};
    }

    size_t rank() const{
        size_t ans = 0;
        auto triRes = this->triangleByGauss();
        Matrix& triangle = triRes.first;
        for(size_t i=0; i<N; ++i) ans += (triangle.body[i][M-1] != Field(0));
        return ans;
    }

    Field det() const{
        auto triRes = this->triangleByGauss();
        Field ans = Field(triRes.second%2 ? -1: 1);
        Matrix& triangle = triRes.first;
        for(size_t i=0; i<N; ++i) ans *= triangle.body[i][i];
        return ans;
    }

    Matrix& invert(){
        check_square<N,Field>(*this);
        Matrix<N, M*2, Field> extended;
        for(size_t i=0; i<N; ++i){
            for(size_t j=0; j<M; ++j){
                extended[i][j] = body[i][j];
            }
            extended[i][i+M] = Field(1);
        }

        extended = extended.triangleByGauss().first;
        for(size_t i=0; i<N; ++i) {
            Field k = Field(1) / extended[i][i];
            for(size_t j=0; j<2*M; ++j) extended[i][j] *= k;
        }
        for(int j=M-1; j>=0; --j){
            for(int i=j-1; i>=0; --i){
                Field mul = extended[i][j] / extended[j][j];
                for(size_t k=j; k<2*M; ++k){
                    extended[i][k] -= mul * extended[j][k];
                }
            }
        }

        for(size_t i=0; i<N; ++i)
            for(size_t j=0; j<M; ++j) body[i][j] = extended[i][j+M];

        return *this;
    }

    Matrix inverted() const{
        Matrix copy = *this;
        return copy.invert();
    }
};

template<size_t N, size_t M, typename Field>
Matrix<N,M,Field> operator -(const Matrix<N,M,Field>& a, const Matrix<N,M,Field>& b){
    Matrix c = a;
    return c -= b;
}

template<size_t N, size_t M, typename Field>
Matrix<N,M,Field> operator +(const Matrix<N,M,Field>& a, const Matrix<N,M,Field>& b){
    Matrix c = a;
    return c += b;
}

template<size_t N, size_t M, typename Field>
Matrix<N,M,Field> operator *(const Field& num, const Matrix<N,M,Field>& b){
    Matrix c = b;
    return c *= num;
}

template<size_t N, size_t M, typename Field>
Matrix<N,M,Field> operator *(const Matrix<N,M,Field>& a, const Field& num){
    Matrix c = a;
    return c *= num;
}

template<size_t N, size_t M, size_t K, typename Field>
Matrix<N,K,Field> operator *(const Matrix<N,M,Field>& a, const Matrix<M,K,Field>& b){
    Matrix<N,K,Field> ans;
    const size_t T = std::max({N,M,K});
    if(T >= 64){
        const size_t exT = T + (T % 2);
        Matrix<exT, exT, Field> A;
        Matrix<exT, exT, Field> B;
        for (size_t i = 0; i < exT; ++i) {
            for (size_t j = 0; j < exT; ++j) {
                A[i][j] = (i < N && j < M) ? a[i][j] : Field(0);
                B[i][j] = (i < M && j < K) ? b[i][j] : Field(0);
            }
        }
        Matrix<exT / 2, exT / 2, Field> A11;
        Matrix<exT / 2, exT / 2, Field> A12;
        Matrix<exT / 2, exT / 2, Field> A21;
        Matrix<exT / 2, exT / 2, Field> A22;
        Matrix<exT / 2, exT / 2, Field> B11;
        Matrix<exT / 2, exT / 2, Field> B12;
        Matrix<exT / 2, exT / 2, Field> B21;
        Matrix<exT / 2, exT / 2, Field> B22;
        for (size_t i = 0; i < exT / 2; ++i) {
            for (size_t j = 0; j < exT / 2; ++j) {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + exT / 2];
                A21[i][j] = A[i + exT / 2][j];
                A22[i][j] = A[i + exT / 2][j + exT / 2];

                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + exT / 2];
                B21[i][j] = B[i + exT / 2][j];
                B22[i][j] = B[i + exT / 2][j + exT / 2];
            }
        }

        Matrix<exT / 2, exT / 2, Field> P1 = (A11 + A22) * (B11 + B22);
        Matrix<exT / 2, exT / 2, Field> P2 = (A21 + A22) * B11;
        Matrix<exT / 2, exT / 2, Field> P3 = A11 * (B12 - B22);
        Matrix<exT / 2, exT / 2, Field> P4 = A22 * (B21 - B11);
        Matrix<exT / 2, exT / 2, Field> P5 = (A11 + A12) * B22;
        Matrix<exT / 2, exT / 2, Field> P6 = (A21 - A11) * (B11 + B12);
        Matrix<exT / 2, exT / 2, Field> P7 = (A12 - A22) * (B21 + B22);

        Matrix<exT / 2, exT / 2, Field> C11 = P1 + P4 - P5 + P7;
        Matrix<exT / 2, exT / 2, Field> C12 = P3 + P5;
        Matrix<exT / 2, exT / 2, Field> C21 = P2 + P4;
        Matrix<exT / 2, exT / 2, Field> C22 = P1 - P2 + P3 + P6;

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < K; ++j) {
                if (i < exT / 2 && j < exT / 2) {
                    ans[i][j] = C11[i][j];
                } else if (i < exT / 2) {
                    ans[i][j] = C12[i][j - exT / 2];
                } else if (j < exT / 2) {
                    ans[i][j] = C21[i - exT / 2][j];
                } else {
                    ans[i][j] = C22[i - exT / 2][j - exT / 2];
                }
            }
        }
        return ans;
    }
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < K; ++j) {
            for (size_t k = 0; k < M; ++k) {
                ans[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return ans;
}

template<size_t N, size_t M, size_t K, size_t T, typename Field>
bool operator==(const Matrix<N,M,Field>& a, const Matrix<K,T,Field>& b){
    if(N != K) return false;
    if(M != T) return  false;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            if (!(a[i][j] == b[i][j])) {
                return false;
            }
        }
    }
    return true;
}

template<size_t N, size_t M, size_t K, size_t T, typename Field>
bool operator!=(const Matrix<N,M,Field>& a, const Matrix<K,T,Field>& b){
    return !(a == b);
}

template<size_t N, size_t M, typename Field>
std::istream& operator>>(std::istream& in, Matrix<N, M, Field>& A) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            in >> A[i][j];
        }
    }
    return in;
}

template<size_t N, size_t M, typename Field>
std::ostream& operator<<(std::ostream& out, const Matrix<N, M, Field>& A) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            out << A[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;
