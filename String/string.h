#include <iostream>
#include <cstring>

class String {
    friend bool operator ==(const String& a, const String& b);

private:
    char* str = new char[1];
    size_t sz = 0;
    size_t capacity = 1;

public:
    String() = default;

    String(size_t n, char c) : str(new char[n]), sz(n), capacity(n) {
        memset(str, c, n);
    }

    String(const String& s) : str(new char[s.sz]), sz(s.sz), capacity(s.sz) {
        memcpy(str, s.str, sz);
    }

    String(const char* s) : str(new char[strlen(s)]), sz(strlen(s)), capacity(strlen(s)) {
        memcpy(str, s, sz);
    }

    String(char c) : String(1, c) {};


    size_t length() const {
        return sz;
    }

    void push_back(char c) {
        if (sz == capacity) {
            char* temp = new char[sz * 2];
            memcpy(temp, str, sz);
            delete[] str;
            str = temp;
            capacity = sz * 2;
        }
        str[sz] = c;
        ++sz;
    }

    void pop_back() {
        --sz;
        if (sz * 4 <= capacity) {
            char* temp = new char[sz * 2];
            memcpy(temp, str, sz);
            delete[] str;
            str = temp;
            capacity = sz * 2;
        }
    }

    char& front() {
        return str[0];
    }

    const char& front() const {
        return str[0];
    }

    char& back() {
        return str[sz - 1];
    }

    const char& back() const {
        return str[sz - 1];
    }

    size_t find(const String& substr) const {
        if (substr.sz > sz) return length();
        for (size_t i = 0; i + substr.sz <= sz; ++i)
            if (memcmp(str + i, substr.str, substr.sz) == 0) return i;
        return length();
    }

    size_t rfind(const String& substr) const {
        if (substr.sz > sz) return length();
        for (long long i = sz - substr.sz; i >= 0; --i)
            if (memcmp(str + i, substr.str, substr.sz) == 0) return i;
        return length();
    }

    String substr(size_t start, size_t count) const {
        String res;
        res.str = new char[count];
        res.sz = count;
        res.capacity = count;
        memcpy(res.str, str + start, count);
        return res;
    }

    bool empty() const {
        return sz == 0;
    }

    void clear() {
        delete[] str;
        str = new char[1];
        sz = 0;
        capacity = 1;
    }

    void swap(String& s) {
        std::swap(str, s.str);
        std::swap(sz, s.sz);
        std::swap(capacity, s.capacity);
    }

    String& operator =(const String& s) {
        String copy = s;
        swap(copy);
        return *this;
    }

    char& operator [](size_t n) {
        return str[n];
    }

    const char& operator [](size_t n) const {
        return str[n];
    }

    String& operator +=(const String& s) {
        if (capacity < sz + s.sz) {
            size_t m = std::max(sz, s.sz);
            char* t = new char[m * 2];
            memcpy(t, str, sz);
            delete[] str;
            str = t;
            capacity = m * 2;
        }
        memcpy(str + sz, s.str, s.sz);
        sz += s.sz;
        return *this;
    }

    ~String() {
        delete[] str;
    }
};

String operator +(const String& a, const String& b) {
    String copy = a;
    return copy += b;
}

bool operator ==(const String& a, const String& b) {
    return (a.length() == b.length()) && (memcmp(a.str, b.str, a.length()) == 0);
}

std::ostream& operator<<(std::ostream& out, const String& s) {
    for (size_t i = 0; i < s.length(); i++)
        out << s[i];
    return out;
}

std::istream& operator>>(std::istream& in, String& s) {
    char c = static_cast<char>(in.get());
    while (!std::isspace(c) && !in.eof()) {
        s.push_back(c);
        c = static_cast<char>(in.get());
    }
    return in;
}
