#ifndef MQ_H
#define MQ_H

#include <stdint.h>

void MQ(uint64_t *fx, const uint64_t *x, const uint64_t *F);
void G(uint64_t *gx, const uint64_t *x, const uint64_t *y, const uint64_t *F);

void MQ3(uint64_t *fx, const uint64_t *x, const uint64_t *F);
void G3(uint64_t *gx, const uint64_t *x, const uint64_t *y, const uint64_t *F);

#endif
