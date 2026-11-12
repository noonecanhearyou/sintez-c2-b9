#pragma once
#include "orf/orf.hpp"
#include "schema.hpp"
#include <QCryptographicHash>
#include <QFile>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
namespace orf {
namespace fs = std::filesystem;
inline std::string str(const J &j, const std::string &k) {
    auto v = j.is_object() ? j.value(k, J()) : J();
    return v.is_string() ? v.get<std::string>() : "";
}
inline int num(const J &j, const std::string &k) {
    auto v = j.is_object() ? j.value(k, J()) : J();
    return v.is_number() ? v.get<int>() : 0;
}
inline void inc(J &j, const std::string &k, int n = 1) {
    j[k] = num(j, k) + n;
}
inline std::string upper(std::string s) {
    for (auto &c : s)
        c = std::toupper((unsigned char)c);
    return s;
}
inline std::string trim(std::string s) {
    auto a = s.find_first_not_of(" \t\r\n");
    return a == s.npos ? "" : s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
inline bool stop(const std::string &s) {
    return s == "TAA" || s == "TAG" || s == "TGA";
}
inline int invalid(const std::string &s) {
    return std::count_if(s.begin(), s.end(), [](char c) {
        return std::string("ACGTRYSWKMBDHVN").find(c) == std::string::npos;
    });
}
inline int ambiguous(const std::string &s) {
    return std::count_if(s.begin(), s.end(), [](char c) {
        return std::string("RYSWKMBDHVN").find(c) != std::string::npos;
    });
}
inline std::string hash(const std::string &s) {
    return QCryptographicHash::hash(QByteArray::fromStdString(s), QCryptographicHash::Sha256)
        .toHex()
        .toStdString();
}
inline std::string fileHash(const fs::path &p) {
    QFile f(QString::fromStdString(p.string()));
    if (!f.open(QIODevice::ReadOnly))
        throw std::runtime_error("Cannot read " + p.string());
    QCryptographicHash h(QCryptographicHash::Sha256);
    if (!h.addData(&f))
        throw std::runtime_error("Hash failed");
    return h.result().toHex().toStdString();
}
inline std::ofstream output(const fs::path &p) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p);
    f.exceptions(std::ios::badbit | std::ios::failbit);
    return f;
}
inline void save(const fs::path &p, const J &j) {
    auto f = output(p);
    f << j.dump(2) << '\n';
}
inline J load(const fs::path &p) {
    std::ifstream f(p);
    if (!f)
        throw std::runtime_error("Cannot read " + p.string());
    return J::parse(f);
}
inline void jsonl(const fs::path &p, const std::function<void(J)> &cb) {
    std::ifstream f(p);
    if (!f)
        throw std::runtime_error("Cannot read " + p.string());
    std::string s;
    int line = 0;
    while (std::getline(f, s)) {
        ++line;
        if (trim(s).empty())
            continue;
        try {
            auto j = J::parse(s);
            if (!j.is_object())
                throw std::runtime_error("Expected object");
            cb(std::move(j));
        } catch (const std::exception &e) {
            throw std::runtime_error(p.string() + ":" + std::to_string(line) + ": " + e.what());
        }
    }
    if (f.bad())
        throw std::runtime_error("Read failed");
}
inline std::string cell(const J &v) {
    if (v.is_null())
        return "";
    if (v.is_boolean())
        return v.get<bool>() ? "True" : "False";
    if (v.is_number_float()) {
        char buffer[64];
        auto result = std::to_chars(buffer, buffer + sizeof(buffer), v.get<double>());
        std::string text(buffer, result.ptr);
        if (text.find_first_of(".eE") == std::string::npos)
            text += ".0";
        return text;
    }
    return v.is_string() ? v.get<std::string>() : v.dump();
}
inline std::string csvQuote(std::string s) {
    if (s.find_first_of(",\"\r\n") == s.npos)
        return s;
    std::string o = "\"";
    for (char c : s) {
        o += c;
        if (c == '\"')
            o += '\"';
    }
    return o + '\"';
}
struct Csv {
    std::ofstream f;
    J fields;
    Csv(const fs::path &p, J h) : f(output(p)), fields(h) {
        J row;
        for (auto &v : fields)
            row[v.get<std::string>()] = v;
        write(row);
    }
    void write(const J &j) {
        bool first = true;
        for (auto &v : fields) {
            if (!first)
                f << ',';
            first = false;
            f << csvQuote(cell(j.value(v.get<std::string>(), J())));
        }
        f << "\r\n";
    }
};
using Key = std::vector<std::string>;
inline Key key(const J &j, std::initializer_list<const char *> fields) {
    Key k;
    for (auto f : fields)
        k.push_back(str(j, f));
    return k;
}
inline J copyFields(const J &j, std::initializer_list<const char *> fields) {
    J r = J::object();
    for (auto f : fields)
        r[f] = j.value(f, J());
    return r;
}
inline J zeros() {
    return J{{"complete", 0}, {"open_end", 0}, {"ambiguous_boundary", 0}};
}
inline double round6(double x) {
    // Round the represented binary double to six decimal places, as Python round(x, 6).
    // Multiplying by 1e6 first would lose information at half-way boundaries.
    char buffer[64];
    auto end = std::to_chars(buffer, buffer + sizeof(buffer), x, std::chars_format::fixed, 6);
    double result = 0;
    std::from_chars(buffer, end.ptr, result);
    return result;
}
struct Context {
    Settings s;
    Progress log;
    const std::atomic_bool *cancel;
    fs::path root;
    void check() {
        if (cancel && cancel->load())
            throw std::runtime_error("Отменено пользователем");
    }
    void tick(int n) {
        check();
        if (s.progressEvery && n % s.progressEvery == 0 && log)
            log("Обработано записей: " + std::to_string(n));
    }
};
void parserStage(Context &);
void scannerStage(Context &);
void comparisonStage(Context &);
void evidenceStage(Context &);
void profileStage(Context &);
} // namespace orf
