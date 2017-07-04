#ifndef PARAMS_H
#define PARAMS_H

#define SOFIA_Q 4
#define SOFIA_N 128
#define SOFIA_M SOFIA_N
#define SOFIA_FLEN (SOFIA_M * (((SOFIA_N * (SOFIA_N + 1)) >> 1) + SOFIA_N))
// when restructuring F such that every chunk aligns to 32 bytes;
#define SOFIA_FLEN_GAPPED (SOFIA_M * (((((SOFIA_N * (SOFIA_N + 1)) + 255) & ~255) >> 1) + SOFIA_N))

// in GF4:
#define SOFIA_NBYTES (SOFIA_N >> 2)
#define SOFIA_MBYTES (SOFIA_M >> 2)
#define SOFIA_FBYTES (SOFIA_FLEN >> 2)
#define SOFIA_N64 (SOFIA_NBYTES >> 3)
#define SOFIA_M64 (SOFIA_MBYTES >> 3)
#define SOFIA_F64 (SOFIA_FBYTES >> 3)
// when restructuring F such that every chunk aligns to 32 bytes;
#define SOFIA_FBYTES_GAPPED (SOFIA_FLEN_GAPPED >> 2)
#define SOFIA_F64_GAPPED (SOFIA_FBYTES_GAPPED >> 3)

#define SOFIA_T 3
#define SOFIA_ROUNDS 438

#define SOFIA_HASHBYTES 32  // TODO test if we can vary this without breaking things
#define SOFIA_SEEDBYTES SOFIA_HASHBYTES  // TODO make usage more consistent / semantically correct

#define SOFIA_SECRETKEYBYTES SOFIA_SEEDBYTES
#define SOFIA_PUBLICKEYBYTES (SOFIA_SEEDBYTES + SOFIA_MBYTES)

// Transcript that is used to generate challenges contains;
//  - Public key pk; SOFIA_PUBLICKEYBYTES
//  - The commits; 2*SOFIA_HASHBYTES*SOFIA_ROUNDS
//  - Length-preserving commits over responses: SOFIA_ROUNDS times
//     - G(resp1): SOFIA_T * (SOFIA_NBYTES + SOFIA_MBYTES
//     - G(resp2): 2 * SOFIA_NBYTES
// It does not explicitly contain the message; that is appended when hashing
#define SOFIA_TRANSCRIPTBYTES (\
    SOFIA_PUBLICKEYBYTES +\
    SOFIA_ROUNDS * (2 * SOFIA_HASHBYTES+\
                    SOFIA_T * (SOFIA_NBYTES + SOFIA_MBYTES) +\
                    2 * SOFIA_NBYTES))

// signature consists of:
//   - SOFIA_HASHBYTES hash over transcript from which challenges are derived
//   - SOFIA_ROUNDS * SOFIA_T responses / commitments for challenge 1
//   - SOFIA_ROUNDS * 2       responses / commitments for challenge 2
//   - SOFIA_ROUNDS * SOFIA_HASHBYTES missing commitments
#define SOFIA_BYTES (SOFIA_HASHBYTES + SOFIA_ROUNDS * (SOFIA_T * (SOFIA_MBYTES + SOFIA_NBYTES) + 2 * SOFIA_NBYTES + SOFIA_HASHBYTES))

#endif
