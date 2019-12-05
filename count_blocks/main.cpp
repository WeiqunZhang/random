#include <iostream>
#include <memory>
#include <cstdio>

#define NUMBLOCKS  512
#define NUMTHREADS 1024

struct BlockMutex
{
    union state_t
    {
        struct { int blockid; int count; };
        unsigned long long ull;
    };

    __host__ __device__
    static constexpr state_t FreeState () noexcept {
        return state_t{{-1,0}};
    }

    BlockMutex () {
        static_assert(sizeof(unsigned long long) == 2*sizeof(int) and
                      sizeof(unsigned long long) == sizeof(state_t),
                      "BlockMutex: wrong size");
        // The first 4 bytes of unsigned long stores blockIdx.
        // The second 4 bytes, count.
        // The initial values are -1 and 0.
        cudaMalloc(&m_state, sizeof(state_t));
        state_t h_state = FreeState();
        cudaMemcpy(m_state, &h_state, sizeof(state_t), cudaMemcpyHostToDevice);
    }

    ~BlockMutex () {
        cudaFree(m_state);
    }

    void operator= (BlockMutex const&) = delete;

    __device__
    void lock () {
        int blockid = blockIdx.z*blockDim.x*blockDim.y + blockIdx.y*blockDim.x + blockIdx.x;
        state_t old = *m_state;
        state_t assumed;
        do {
            assumed = old;
            state_t val;
            val.blockid = blockid;
            if (assumed.blockid == blockid) {
                // Already locked by another thread in this block. Need to ++count.
                val.count = assumed.count + 1;
            } else {
                // Currently unlocked or locked by another block.  Need to lock.
                val.count = 1;
                assumed = FreeState();
            }
            old.ull = atomicCAS((unsigned long long*)m_state, assumed.ull, val.ull);
        } while (assumed.ull != old.ull);
    }

    __device__
    void unlock () {
        state_t old = *m_state;
        state_t assumed;
        do {
            assumed = old;
            state_t val;
            if (assumed.count == 1) {
                // Need to unlock
                val = FreeState();
            } else {
                // --count, but do NOT unlock
                val = assumed;
                --val.count;
            }
            old.ull = atomicCAS((unsigned long long*)m_state, assumed.ull, val.ull);
        } while (assumed.ull != old.ull);
    }

private:

    state_t* m_state;
};

__global__
void test (BlockMutex* mut, int *numBlocks)
{
    mut->lock();
    if (threadIdx.x == 0) { numBlocks[0] = numBlocks[0] + 1; }
    __threadfence();
    mut->unlock();
}

int main (int argc, char* argv[])
{
    for (int i = 0; i < 10; ++i) {
        int h_counting, *d_counting;
        cudaMalloc(&d_counting, sizeof(int));
        h_counting = 0;
        cudaMemcpy(d_counting, &h_counting, sizeof(int), cudaMemcpyHostToDevice);

        BlockMutex h_mut, *d_mut;
        cudaMalloc(&d_mut, sizeof(BlockMutex));
        cudaMemcpy(d_mut, &h_mut, sizeof(BlockMutex), cudaMemcpyHostToDevice);

        test<<<NUMBLOCKS,NUMBLOCKS>>>(d_mut, d_counting);

        cudaPeekAtLastError();

        cudaDeviceSynchronize();

        cudaMemcpy(&h_counting, d_counting, sizeof(int), cudaMemcpyDeviceToHost);

        std::cout << "Number of blocks is: " << h_counting << std::endl;
        if (h_counting != NUMBLOCKS) std::abort();
        
        cudaFree(d_counting);
    }
}
