#include "kwip_kmerhash.h"

#include <string.h>
#include <assert.h>
#include <stdio.h>

#include <xxhash.h>



void
kmer_iter_init(kmer_iter_t *ctx, size_t k, uint64_t seed, bool canonicalise)
{
    assert(ctx != NULL);
    assert(k > 0 && k <= 32);
    ctx->k = k;
    ctx->mask = (1 << (2*k)) - 1;
    ctx->canonicalise = canonicalise;
    ctx->seq = NULL;
    ctx->len = 0;
    ctx->i = 0;
    ctx->rcseq =  NULL;
    ctx->seed = seed;
}

void
kmer_iter_set_seq(kmer_iter_t *ctx, char *seq, size_t seqlen)
{
    assert(ctx != NULL);
    assert(seq != NULL);
    ctx->i = 0;
    ctx->seq = seq;
    for (size_t i = 0; i < seqlen; i++ ) {
        ctx->seq[i] &= 0x5f;  // Force uppercase
    }
    if (ctx->canonicalise) {
        if (seqlen > ctx->len) {
            ctx->rcseq = realloc(ctx->rcseq, seqlen + 1);
        }
        assert(ctx->rcseq != NULL);
        for (size_t i = 0; i < seqlen; i++) {
            char nt = seq[seqlen - i - 1];
            switch (nt) {
                case 'A':
                    ctx->rcseq[i] = 'T';
                    break;
                case 'C':
                    ctx->rcseq[i] = 'G';
                    break;
                case 'G':
                    ctx->rcseq[i] = 'C';
                    break;
                case 'T':
                    ctx->rcseq[i] = 'A';
                    break;
                default:
                    ctx->rcseq[i] = 'N';
                    break;
            }
        }
        ctx->rcseq[seqlen] = '\0';
    } else {
        ctx->rcseq = NULL;
    }
    ctx->len = seqlen;
}

int
kmer_iter_next_nthash(kmer_iter_t *ctx, uint64_t *hash)
{
    // Numeric encoding of kmer into uint64 in 2-bit encoding, lexographically
    // ordered
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
        char nucl = ctx->seq[ctx->i++];
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
kmer_iter_next_xxh(kmer_iter_t *ctx, uint64_t *hash)
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
        char nucl = ctx->seq[ctx->i++];  // Force uppercase
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
    size_t offset = ctx->i - k;
    if (ctx->canonicalise) {
        uint64_t fwd, rev;
        fwd = XXH64(ctx->seq + offset, k, ctx->seed);
        rev = XXH64(ctx->rcseq + offset, k, ctx->seed);
        *hash = fwd > rev ? rev : fwd;
    } else {
        *hash = XXH64(ctx->seq + offset, k, ctx->seed);
    }
    return 1;
}

void
kmer_iter_destroy(kmer_iter_t *ctx)
{
    if (ctx->rcseq != NULL) {
        free(ctx->rcseq);
        ctx->rcseq = NULL;
    }
    ctx->seq = NULL;
    ctx->len = 0;
}

uint64_t
kmer_xxh(char *seq, size_t len, uint64_t seed, bool canonicalise)
{
    uint64_t hash;
    kmer_iter_t itr;
    kmer_iter_init(&itr, len, seed, canonicalise);
    kmer_iter_set_seq(&itr, seq, len);

    kmer_iter_next_xxh(&itr, &hash);

    kmer_iter_destroy(&itr);
    return hash;
}
