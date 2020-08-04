
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
    std::array<uint32_t,2> m_morton;
};

bool
SFCToken::Compare::operator () (const SFCToken& lhs,
                                const SFCToken& rhs) const
{
        return (lhs.m_morton[1] <  rhs.m_morton[1]) ||
              ((lhs.m_morton[1] == rhs.m_morton[1]) &&
               (lhs.m_morton[0] <  rhs.m_morton[0]));
}

uint32_t make_space (uint32_t x)
{
    // x           : 0000,0000,0000,0000,gfed,cba9,8765,4321
    x = (x | (x << 8)) & 0x00FF00FF;
    // x << 8      : 0000,0000,gfed,cba9,8765,4321,0000,0000
    // x | (x << 8): 0000,0000,gfed,cba9,????,????,8765,4321
    // 0x00FF00FF  : 0000,0000,1111,1111,0000,0000,1111,1111
    // x           : 0000,0000,gfed,cba9,0000,0000,8765,4321
    x = (x | (x << 4)) & 0x0F0F0F0F;
    // x << 4      : 0000,gfed,cba9,0000,0000,8765,4321,0000
    // x | (x << 4): 0000,gfed,????,cba9,0000,8765,????,4321
    // 0x0F0F0F0F  : 0000,1111,0000,1111,0000,1111,0000,1111
    // x           : 0000,gfed,0000,cba9,0000,8765,0000,4321
    x = (x | (x << 2)) & 0x33333333;
    // x << 2      : 00gf,ed00,00cb,a900,0087,6500,0043,2100
    // x | (x << 2): 00gf,??ed,00cb,??a9,0087,??65,0043,??21
    // 0x33333333  : 0011,0011,0011,0011,0011,0011,0011,0011
    // x           : 00gf,00ed,00cb,00a9,0087,0065,0043,0021
    x = (x | (x << 1)) & 0x55555555;
    // x << 1      : 0gf0,0ed0,0cb0,0a90,0870,0650,0430,0210
    // x | (x << 1): 0g?f,0e?d,0c?b,0a?9,08?7,06?5,04?3,02?1
    // 0x55555555  : 0101,0101,0101,0101,0101,0101,0101,0101
    // x           : 0g0f,0e0d,0c0b,0a09,0807,0605,0403,0201
    return x;
}

SFCToken makeSFCToken (int i, std::array<int,2> const& iv)
{
    SFCToken token;
    token.index = i;

    constexpr uint32_t offset = 1 << 31;
    static_assert(static_cast<uint32_t>(std::numeric_limits<int>::max())+1 == offset,
                  "INT_MAX != (1<<31)-1");
    uint32_t x = (iv[0] >= 0) ? static_cast<uint32_t>(iv[0]) + offset
        : static_cast<uint32_t>(iv[0]-std::numeric_limits<int>::lowest());
    uint32_t y = (iv[1] >= 0) ? static_cast<uint32_t>(iv[1]) + offset
        : static_cast<uint32_t>(iv[1]-std::numeric_limits<int>::lowest());
    // extract lowest 16 bits and make sapce for interleaving
    token.m_morton[0] = make_space(x & 0xFFFF)
        | (make_space(y & 0xFFFF) << 1);
    x = x >> 16;
    y = y >> 16;
    token.m_morton[1] = make_space(x) | (make_space(y) << 1);

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

                for (int j = 1; j >= 0; --j)
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
    std::array<int,2> m_iv;
};

OldSFCToken makeOldSFCToken (int i, std::array<int,2> const& iv)
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
        tokens.push_back(makeSFCToken(i,std::array<int,2>{x,y}));
        tokens2.push_back(makeOldSFCToken(i,std::array<int,2>{x,y}));
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
