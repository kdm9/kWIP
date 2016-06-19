#include <string.h>
#include <assert.h>

#include "kmerhash.h"
#include "xxhash.h"

void
kmer_iter_init(kmer_iter_t *ctx, const char *seq, size_t n, size_t k)
{
    assert(ctx != NULL);
    assert(k > 1 && k <= 32);
    ctx->seq = seq;
    ctx->len = n;
    ctx->i = 0;
    ctx->k = k;
    ctx->mask = (1 << (2*k)) - 1;
}

int
kmer_iter_next(kmer_iter_t *ctx, uint64_t *hash)
{
    assert(ctx != NULL);
    assert(hash != NULL);
    size_t len = ctx->len, k = ctx->k;
    size_t skip = 0;
    do {
        if (ctx->i >= len) {
            /* Prevent an out-of-bounds read */
            return 0;
        }
        if (skip > 0) {
            skip--;
        }
        char nucl = ctx->seq[ctx->i++] & 0x5f;  // Force uppercase
        uint64_t n = 0; // Numeric nucleotide repr
        switch (nucl) {
            case 'A':
                n = 0;
                break;
            case 'C':
                n = 1;
                break;
            case 'G':
                n = 2;
                break;
            case 'T':
                n = 3;
                break;
            default:
                skip = k;
                break;
        }
        ctx->last_hash = ((ctx->last_hash << 2) | n) & ctx->mask;
    } while (skip > 0 || ctx->i < k);
    *hash = ctx->last_hash;
    return 1;
}

int
kmer_iter_next_xxh(kmer_iter_t *ctx, uint64_t *hash, uint64_t seed)
{
    assert(ctx != NULL);
    assert(hash != NULL);
    size_t len = ctx->len, k = ctx->k;
    size_t skip = 0;
    do {
        if (ctx->i >= len) {
            /* Prevent an out-of-bounds read */
            return 0;
        }
        char nucl = ctx->seq[ctx->i++] & 0x5f;  // Force uppercase
        switch (nucl) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                if (skip > 0) {
                    skip--;
                }
                break;
            default:
                skip = k;
                break;
        }
    } while (skip > 0 || ctx->i < k);
    size_t offset = ctx->i - k + 1;
    *hash = XXH64(ctx->seq + offset, k, seed);
    return 1;
}
void
kmer_iter_destroy(kmer_iter_t *ctx)
{
}
