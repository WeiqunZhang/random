
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

struct SFCToken
{
    class Compare
    {
    public:
        bool operator () (const SFCToken& lhs,
                          const SFCToken& rhs) const;
    };
    int index;
    std::array<uint32_t,3> m_morton;
};

bool
SFCToken::Compare::operator () (const SFCToken& lhs,
                                const SFCToken& rhs) const
{
        return (lhs.m_morton[2] <  rhs.m_morton[2]) ||
              ((lhs.m_morton[2] == rhs.m_morton[2]) &&
              ((lhs.m_morton[1] <  rhs.m_morton[1]) ||
              ((lhs.m_morton[1] == rhs.m_morton[1]) &&
               (lhs.m_morton[0] <  rhs.m_morton[0]))));
}

uint32_t make_space (uint32_t x)
{
    // x            : 0000,0000,0000,0000,0000,00a9,8765,4321
    x = (x | (x << 16)) & 0x030000FF;
    // x << 16      : 0000,00a9,8765,4321,0000,0000,0000,0000
    // x | (x << 16): 0000,00a9,8765,4321,0000,00a9,8765,4321
    // 0x030000FF   : 0000,0011,0000,0000,0000,0000,1111,1111
    // x            : 0000,00a9,0000,0000,0000,0000,8765,4321
    x = (x | (x <<  8)) & 0x0300F00F;
    // x << 8       : 0000,0000,0000,0000,8765,4321,0000,0000
    // x | (x << 8) : 0000,00a9,0000,0000,8765,4321,8765,4321
    // 0x0300F00F   : 0000,0011,0000,0000,1111,0000,0000,1111
    // x            : 0000,00a9,0000,0000,8765,0000,0000,4321
    x = (x | (x <<  4)) & 0x030C30C3;
    // x << 4       : 00a9,0000,0000,8765,0000,0000,4321,0000
    // x | (x << 4) : 00a9,00a9,0000,8765,8765,0000,4321,4321
    // 0x030C30C3   : 0000,0011,0000,1100,0011,0000,1100,0011
    // x            : 0000,00a9,0000,8700,0065,0000,4300,0021
    x = (x | (x <<  2)) & 0x09249249;
    // x << 2       : 0000,a900,0087,0000,6500,0043,0000,2100
    // x | (x << 2) : 0000,a9a9,0087,8700,6565,0043,4300,2121
    // 0x09249249   : 0000,1001,0010,0100,1001,0010,0100,1001
    // x            : 0000,a009,0080,0700,6005,0040,0300,2001
    return x;
}

SFCToken makeSFCToken (int i, std::array<int,3> const& iv)
{
    SFCToken token;
    token.index = i;

    constexpr int imin = -static_cast<int>(1 << 29);
    uint32_t x = iv[0] - imin;
    uint32_t y = iv[1] - imin;
    uint32_t z = iv[2] - imin;
    // extract lowest 10 bits and make space for interleaving
    token.m_morton[0] = make_space(x & 0x3FF)
                     | (make_space(y & 0x3FF) << 1)
                     | (make_space(z & 0x3FF) << 2);
    x = x >> 10;
    y = y >> 10;
    z = z >> 10;
    token.m_morton[1] = make_space(x & 0x3FF)
                     | (make_space(y & 0x3FF) << 1)
                     | (make_space(z & 0x3FF) << 2);
    x = x >> 10;
    y = y >> 10;
    z = z >> 10;
    token.m_morton[2] = make_space(x & 0x3FF)
                     | (make_space(y & 0x3FF) << 1)
                     | (make_space(z & 0x3FF) << 2);

    return token;
}


struct OldSFCToken
{
    struct Compare
    {
        bool operator() (const OldSFCToken& lhs, const OldSFCToken& rhs)
        {
            for (int i = 29; i >= 0; --i)
            {
                const int N = (1<<i);

                for (int j = 2; j >= 0; --j)
                {
                    const int il = lhs.m_iv[j]/N;
                    const int ir = rhs.m_iv[j]/N;

                    if (il < ir)
                    {
                        return true;
                    }
                    else if (il > ir)
                    {
                        return false;
                    }
                }
            }
            return false;
        }
    };
    int index;
    std::array<int,3> m_iv;
};

OldSFCToken makeOldSFCToken (int i, std::array<int,3> const& iv)
{
    return OldSFCToken{i, iv};
}

int main (int argc, char* argv[])
{
    constexpr int N = 100000000;
    std::mt19937 gen;
    constexpr int imin = 0; // -(1 << 29);
    constexpr int imax = (1<<29) - 1;
    std::uniform_int_distribution<int> distrib(imin, imax);
    std::vector<SFCToken> tokens;
    std::vector<OldSFCToken> tokens2;
    for (int i = 0; i < N; ++i) {
        auto x = distrib(gen);
        auto y = distrib(gen);
        auto z = distrib(gen);
        tokens.push_back(makeSFCToken(i,std::array<int,3>{x,y,z}));
        tokens2.push_back(makeOldSFCToken(i,std::array<int,3>{x,y,z}));
    }
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::sort(tokens.begin(), tokens.end(), SFCToken::Compare());
        double dt = std::chrono::duration_cast<std::chrono::duration<double> >
            (std::chrono::high_resolution_clock::now() - t0).count();
        std::cout << "New time: " << dt << "\n";
    }
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::sort(tokens2.begin(), tokens2.end(), OldSFCToken::Compare());
        double dt = std::chrono::duration_cast<std::chrono::duration<double> >
            (std::chrono::high_resolution_clock::now() - t0).count();
        std::cout << "Old time: " << dt << "\n";
    }
    for (int i = 0; i < N; ++i) {
        if (tokens[i].index != tokens2[i].index) {
            std::cout << "WRONG: " << i << " " << tokens[i].index << " " << tokens2[i].index
                      << std::endl;
            exit(1);
        }
    }
}
