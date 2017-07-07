#ifndef F_H
#define F_H

#include <stdint.h>

void expand_F(uint64_t *F, const unsigned char *seed);
void expand_F_gapped(uint64_t *F, const unsigned char *seed);

#endif
