
#ifndef MQ4_H
#define MQ4_H

#include <stdint.h>


#ifdef  __cplusplus
extern  "C" {
#endif

void mq_gf4_n128_m128_ref( uint64_t * fx , const uint64_t * x , const uint64_t * F );
void mq_gf4_n128_m128( uint64_t * fx , const uint64_t * x , const uint64_t * F );
void mq_gf4_n128_m128_x3( uint64_t * fx , const uint64_t * x , const uint64_t * F );
void mq_gf4_n128_m128_x4( uint64_t * fx , const uint64_t * x , const uint64_t * F );
void mq_gf4_n128_m128_vartime( uint64_t * fx , const uint64_t * x , const uint64_t * F );


void G_gf4_n128_m128_ref( uint64_t * fx , const uint64_t * x , const uint64_t * y , const uint64_t * F  );
void G_gf4_n128_m128( uint64_t * fx , const uint64_t * x , const uint64_t * y , const uint64_t * F  );
void G_gf4_n128_m128_x3( uint64_t * fx , const uint64_t * x , const uint64_t * y , const uint64_t * F  );



#ifdef  __cplusplus
}
#endif


#endif
