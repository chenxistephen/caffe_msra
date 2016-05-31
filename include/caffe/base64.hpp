#pragma once

#include <vector>

namespace IU
{

    /*
    const char* PrintDecodeTable()
    {
    const char* revtable = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    std::vector<char> table;
    table.assign(256, (char)(-1));
    for (size_t i = 0; i < 64; ++i)
    {
    table[(unsigned char)(revtable[i])] = (char)(i);
    }
    table[(unsigned char)('=')] = (char)(0);
    int k = 0;
    for (int i = 0; i < 16; i++)
    {
    for (int j = 0; j<16; j++)
    std::cout << (int)table[k++] << ", ";
    std::cout << std::endl;
    }
    } */
    inline bool Base64Decode(const std::string& str, std::vector<unsigned char>& data)
    {
        const char table[] = {
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -1, -1, 63,
            52, 53, 54, 55, 56, 57, 58, 59, 60, 61, -1, -1, -1, 0, -1, -1,
            -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
            15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, -1, -1, -1, -1, -1,
            -1, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
            41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
        };

        if (str.size() == 0)
        {
            data.resize(0);
            return true;
        }
        if (str.size() % 4 != 0) return false;
        data.resize(3 * (str.size() / 4) - ((str[str.size() - 1] != '=') ? 0 : ((str[str.size() - 2] != '=') ? 1 : 2)));
        size_t i = 0;
        size_t j = 0;
        while (i < str.size())
        {
            unsigned long a = table[str[i++]];
            unsigned long b = table[str[i++]];
            unsigned long c = table[str[i++]];
            unsigned long d = table[str[i++]];
            if (a >= 64 || b >= 64 || c >= 64 || d >= 64) return false;
            unsigned long triple = (a << 18) + (b << 12) + (c << 6) + d;
            if (j < data.size()) data[j++] = (triple >> 16) & 255;
            if (j < data.size()) data[j++] = (triple >> 8) & 255;
            if (j < data.size()) data[j++] = triple & 255;
        }
        return true;
    }

    inline void Base64Encode(const std::vector<unsigned char>& data, std::string& str)
    {
        const char* table = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        str.resize(4 * ((data.size() + 2) / 3));
        size_t i = 0;
        size_t j = 0;
        while (i < data.size())
        {
            unsigned long a = (i < data.size()) ? data[i++] : 0;
            unsigned long b = (i < data.size()) ? data[i++] : 0;
            unsigned long c = (i < data.size()) ? data[i++] : 0;
            unsigned long triple = (a << 16) + (b << 8) + c;
            str[j++] = table[(triple >> 18) & 63];
            str[j++] = table[(triple >> 12) & 63];
            str[j++] = table[(triple >> 6) & 63];
            str[j++] = table[triple & 63];
        }
        j -= ((data.size() % 3 == 0) ? 0 : ((data.size() % 3 == 2) ? 1 : 2));
        while (j < str.size()) str[j++] = '=';
    }

}
