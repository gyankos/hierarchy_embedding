//
// Created by giacomo on 29/12/19.
//

#include "string_utils.h"

std::vector<std::string> string_split_to_stringvector(const std::string &str, const std::string &delim) {
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

std::vector<size_t> string_split_to_sizetvector(const std::string &str, const std::string &delim) {
    std::vector<size_t> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(std::stoull(token));
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

std::string size_vector_to_string(const std::vector<size_t> &vector) {
    StringBuilder<char> sb;
    for (const size_t& x : vector)
        sb.Append(std::to_string(x));
    return sb.Join("_");
}

std::vector<std::vector<size_t>> generateAllPossibleSubpaths(const std::vector<size_t> &x) {
    std::vector<std::vector<size_t>> toReturn;
    toReturn.emplace_back(std::vector<size_t>{});
    for (size_t i = 1, n = x.size(); i<=n; i++) {
        toReturn.emplace_back(std::vector<size_t>(x.begin(), x.begin()+i));
    }
    return toReturn;
}
