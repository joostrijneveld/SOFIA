#include <stdint.h>


#include "emmintrin.h"
#include "tmmintrin.h"

#include <immintrin.h>



#ifdef  __cplusplus
extern  "C" {
#endif



static uint8_t extract_bit[32] __attribute__((aligned(32))) = {
0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80, 0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80, 0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80, 0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80 };

static uint32_t extract_01byte[8] __attribute__((aligned(32))) = {0,0,0x01010101,0x01010101, 0,0,0x01010101,0x01010101};

/*
static uint32_t high_lane[8] __attribute__((aligned(32))) = {0,0,0,0, 0xffffffff,0xffffffff,0xffffffff,0xffffffff};
*/

static uint8_t _mask_1[32] __attribute__((aligned(32))) = {1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1};




static void _mq_gf4_n128_m128( uint8_t * z , const uint8_t * w , const uint8_t * pk_mat  )
{
	__m256i zero = _mm256_setzero_si256();
	__m256i r0 = _mm256_setzero_si256();
	__m256i r1 = _mm256_setzero_si256();

	__m256i inp = _mm256_load_si256( (__m256i const *) w );
	unsigned i,j,k;

	/* vertical broadcast bit */
	__m256i x0br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x0br[(i+j)] = br0;
		}
	}

	const __m256i *eqs = (const __m256i *) pk_mat;
	/* linear terms */
	for(i=0;i<128;i++) {
		__m256i br0 = x0br[i];
		__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
		r0 ^= ( br0 & (*eqs) );
		r1 ^= ( br1 & (*eqs) );
		eqs += 1;
	}

	/* quadratic terms */
	for(k=0;k<128;k++) {
		__m256i br0;
		__m256i br1;
		__m256i t[2];
		t[0] = r0;
		t[1] = r1;
		r0 ^= r0;
		r1 ^= r1;
		for(i=0;i<=k;i++) {
			br0 = x0br[i];
			br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			eqs += 1;
		}
		__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
		r0 = (rr&br0);
		r1 = (rr&br1);
		r0 ^= t[0];
		r1 ^= t[1];
	}

	///  r0 = r11, r00 ,  r1 = r10 , r01
        ///  r = (r10,r01,r11) , (r00,r11)

	__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)z , rr );
}





static void _mq_gf4_n128_m128_x3( uint8_t * z , const uint8_t * w , const uint8_t * pk_mat )
{
	__m256i zero = _mm256_setzero_si256();
	__m256i r0 = _mm256_setzero_si256();
	__m256i r1 = _mm256_setzero_si256();
	__m256i r2 = _mm256_setzero_si256();
	__m256i r3 = _mm256_setzero_si256();
	__m256i r4 = _mm256_setzero_si256();
	__m256i r5 = _mm256_setzero_si256();

	__m256i inp = _mm256_load_si256( (__m256i const *) w );
	/// vertical broadcast bit
	__m256i x0br[128];
	unsigned i,j,k;
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x0br[(i+j)] = br0;
		}
	}
	inp = _mm256_load_si256( (__m256i const *) (w+32) );
	__m256i x1br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x1br[(i+j)] = br0;
		}
	}

	inp = _mm256_load_si256( (__m256i const *) (w+64) );
	__m256i x2br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x2br[(i+j)] = br0;
		}
	}


	const __m256i *eqs = (const __m256i *) pk_mat;
	/// linear terms
	for(i=0;i<128;i++) {
		__m256i br0 = x0br[i];
		__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
		__m256i br2 = x1br[i];
		__m256i br3 = _mm256_permute4x64_epi64( br2 , 0x4e );
		__m256i br4 = x2br[i];
		__m256i br5 = _mm256_permute4x64_epi64( br4 , 0x4e );
		r0 ^= ( br0 & (*eqs) );
		r1 ^= ( br1 & (*eqs) );
		r2 ^= ( br2 & (*eqs) );
		r3 ^= ( br3 & (*eqs) );
		r4 ^= ( br4 & (*eqs) );
		r5 ^= ( br5 & (*eqs) );
		eqs += 1;
	}

	/// quadratic terms
	for(k=0;k<128;k++) {
		__m256i br0;
		__m256i br1;
		__m256i br2;
		__m256i br3;
		__m256i br4;
		__m256i br5;
		__m256i t[6];
		t[0] = r0;
		t[1] = r1;
		t[2] = r2;
		t[3] = r3;
		t[4] = r4;
		t[5] = r5;
		r0 ^= r0;
		r1 ^= r1;
		r2 ^= r2;
		r3 ^= r3;
		r4 ^= r4;
		r5 ^= r5;
		for(i=0;i<=k;i++) {
			br0 = x0br[i];
			br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			br2 = x1br[i];
			br3 = _mm256_permute4x64_epi64( br2 , 0x4e );
			br4 = x2br[i];
			br5 = _mm256_permute4x64_epi64( br4 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			r2 ^= ( br2 & (*eqs) );
			r3 ^= ( br3 & (*eqs) );
			r4 ^= ( br4 & (*eqs) );
			r5 ^= ( br5 & (*eqs) );
			eqs += 1;
		}
		__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
		r0 = (rr&br0);
		r1 = (rr&br1);

		rr =  _mm256_permute2x128_si256(r2,r3,0x21) ^ r2 ^ _mm256_permute2x128_si256(r3,r3,0x18);
		r2 = (rr&br2);
		r3 = (rr&br3);

		rr =  _mm256_permute2x128_si256(r4,r5,0x21) ^ r4 ^ _mm256_permute2x128_si256(r5,r5,0x18);
		r4 = (rr&br4);
		r5 = (rr&br5);
		r0 ^= t[0];
		r1 ^= t[1];
		r2 ^= t[2];
		r3 ^= t[3];
		r4 ^= t[4];
		r5 ^= t[5];
	}
	__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)z , rr );

	r0 = r2;
	r1 = r3;
	rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)(z+32) , rr );

	r0 = r4;
	r1 = r5;
	rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)(z+64) , rr );
}




static void _mq_gf4_n128_m128_x4( uint8_t * z , const uint8_t * w , const uint8_t * pk_mat )
{
	__m256i zero = _mm256_setzero_si256();
	__m256i r0 = _mm256_setzero_si256();
	__m256i r1 = _mm256_setzero_si256();
	__m256i r2 = _mm256_setzero_si256();
	__m256i r3 = _mm256_setzero_si256();
	__m256i r4 = _mm256_setzero_si256();
	__m256i r5 = _mm256_setzero_si256();
	__m256i r6 = _mm256_setzero_si256();
	__m256i r7 = _mm256_setzero_si256();

	__m256i inp = _mm256_load_si256( (__m256i const *) w );
	unsigned i,j,k;
	/// vertical broadcast bit
	__m256i x0br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x0br[(i+j)] = br0;
		}
	}
	inp = _mm256_load_si256( (__m256i const *) (w+32) );
	__m256i x1br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x1br[(i+j)] = br0;
		}
	}

	inp = _mm256_load_si256( (__m256i const *) (w+64) );
	__m256i x2br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x2br[(i+j)] = br0;
		}
	}
	inp = _mm256_load_si256( (__m256i const *) (w+96) );
	__m256i x3br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x3br[(i+j)] = br0;
		}
	}


	const __m256i *eqs = (const __m256i *) pk_mat;
	/// linear terms
	for(i=0;i<128;i++) {
		__m256i br0 = x0br[i];
		__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
		__m256i br2 = x1br[i];
		__m256i br3 = _mm256_permute4x64_epi64( br2 , 0x4e );
		__m256i br4 = x2br[i];
		__m256i br5 = _mm256_permute4x64_epi64( br4 , 0x4e );
		__m256i br6 = x3br[i];
		__m256i br7 = _mm256_permute4x64_epi64( br6 , 0x4e );
		r0 ^= ( br0 & (*eqs) );
		r1 ^= ( br1 & (*eqs) );
		r2 ^= ( br2 & (*eqs) );
		r3 ^= ( br3 & (*eqs) );
		r4 ^= ( br4 & (*eqs) );
		r5 ^= ( br5 & (*eqs) );
		r6 ^= ( br6 & (*eqs) );
		r7 ^= ( br7 & (*eqs) );
		eqs += 1;
	}

	/// quadratic terms
	for(k=0;k<128;k++) {
		__m256i br0;
		__m256i br1;
		__m256i br2;
		__m256i br3;
		__m256i br4;
		__m256i br5;
		__m256i br6;
		__m256i br7;
		__m256i t[8];
		t[0] = r0;
		t[1] = r1;
		t[2] = r2;
		t[3] = r3;
		t[4] = r4;
		t[5] = r5;
		t[6] = r6;
		t[7] = r7;
		r0 ^= r0;
		r1 ^= r1;
		r2 ^= r2;
		r3 ^= r3;
		r4 ^= r4;
		r5 ^= r5;
		r6 ^= r6;
		r7 ^= r7;
		for(i=0;i<=k;i++) {
			br0 = x0br[i];
			br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			br2 = x1br[i];
			br3 = _mm256_permute4x64_epi64( br2 , 0x4e );
			br4 = x2br[i];
			br5 = _mm256_permute4x64_epi64( br4 , 0x4e );
			br6 = x3br[i];
			br7 = _mm256_permute4x64_epi64( br6 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			r2 ^= ( br2 & (*eqs) );
			r3 ^= ( br3 & (*eqs) );
			r4 ^= ( br4 & (*eqs) );
			r5 ^= ( br5 & (*eqs) );
			r6 ^= ( br6 & (*eqs) );
			r7 ^= ( br7 & (*eqs) );
			eqs += 1;
		}
		__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
		r0 = (rr&br0);
		r1 = (rr&br1);

		rr =  _mm256_permute2x128_si256(r2,r3,0x21) ^ r2 ^ _mm256_permute2x128_si256(r3,r3,0x18);
		r2 = (rr&br2);
		r3 = (rr&br3);

		rr =  _mm256_permute2x128_si256(r4,r5,0x21) ^ r4 ^ _mm256_permute2x128_si256(r5,r5,0x18);
		r4 = (rr&br4);
		r5 = (rr&br5);

		rr =  _mm256_permute2x128_si256(r6,r7,0x21) ^ r6 ^ _mm256_permute2x128_si256(r7,r7,0x18);
		r6 = (rr&br6);
		r7 = (rr&br7);

		r0 ^= t[0];
		r1 ^= t[1];
		r2 ^= t[2];
		r3 ^= t[3];
		r4 ^= t[4];
		r5 ^= t[5];
		r6 ^= t[6];
		r7 ^= t[7];
	}
	__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)z , rr );

	r0 = r2;
	r1 = r3;
	rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)(z+32) , rr );

	r0 = r4;
	r1 = r5;
	rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)(z+64) , rr );

	r0 = r6;
	r1 = r7;
	rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)(z+96) , rr );
}





static void _G_gf4_n128_m128( uint8_t * z , const uint8_t * w , const uint8_t * _y , const uint8_t * pk_mat )
{
	__m256i zero = _mm256_setzero_si256();
	__m256i r0 = _mm256_setzero_si256();
	__m256i r1 = _mm256_setzero_si256();


	__m256i inp = _mm256_load_si256( (__m256i const *) w );
	unsigned i,j,k;
	/// vertical broadcast bit
	__m256i x0br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x0br[(i+j)] = br0;
		}
	}

	inp = _mm256_load_si256( (__m256i const *) _y );
	/// vertical broadcast bit
	__m256i y0br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			y0br[(i+j)] = br0;
		}
	}

	__m256i x = _mm256_load_si256( (__m256i const *) w );
	__m256i y = _mm256_load_si256( (__m256i const *) _y );
	__m256i x_usd = _mm256_permute4x64_epi64( x , 0x4e );
	__m256i y_usd = _mm256_permute4x64_epi64( y , 0x4e );

	const __m256i *eqs = (const __m256i *) pk_mat;
	/// linear terms
	eqs += 128;

	__m256i mask_1 = _mm256_load_si256((__m256i*)_mask_1);
	/// quadratic terms
	eqs += 1;
	for(k=1;k<64;k++) {
		__m256i a0 = (x & y0br[k]) ^ (y & x0br[k]);
		__m256i a1 = (x_usd & y0br[k] ) ^ (y_usd & x0br[k]);
		__m256i quad_terms = _mm256_permute2x128_si256(a0,a1,0x21) ^ a0 ^ _mm256_permute2x128_si256(a1,a1,0x18);
		for(i=0;i<k;i++) {
			__m256i br0 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms,zero));
			quad_terms = _mm256_srli_epi64( quad_terms , 1 );
			__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			eqs += 1;
		}
		eqs += 1;
	}
	for(k=64;k<128;k++) {
		__m256i a0 = (x & y0br[k]) ^ (y & x0br[k]);
		__m256i a1 = (x_usd & y0br[k] ) ^ (y_usd & x0br[k]);
		__m256i quad_terms = _mm256_permute2x128_si256(a0,a1,0x21) ^ a0 ^ _mm256_permute2x128_si256(a1,a1,0x18);
		__m256i quad_terms_h64 = _mm256_srli_si256( quad_terms , 8 );
		for(i=0;i<64;i++) {
			__m256i br0 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms,zero));
			quad_terms = _mm256_srli_epi64( quad_terms , 1 );
			__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			eqs += 1;
		}
		quad_terms = quad_terms_h64;
		for(i=64;i<k;i++) {
			__m256i br0 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms,zero));
			quad_terms = _mm256_srli_epi64( quad_terms , 1 );
			__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			eqs += 1;
		}
		eqs += 1;
	}
	__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)z , rr );
}










static void _G_gf4_n128_m128_x3( uint8_t * z , const uint8_t * w , const uint8_t * _y , const uint8_t * pk_mat )
{
	__m256i zero = _mm256_setzero_si256();
	__m256i r0 = _mm256_setzero_si256();
	__m256i r1 = _mm256_setzero_si256();
	__m256i r2 = _mm256_setzero_si256();
	__m256i r3 = _mm256_setzero_si256();
	__m256i r4 = _mm256_setzero_si256();
	__m256i r5 = _mm256_setzero_si256();


	__m256i x0 = _mm256_load_si256( (__m256i const *) w );
	__m256i x1 = _mm256_load_si256( (__m256i const *) (w+32) );
	__m256i x2 = _mm256_load_si256( (__m256i const *) (w+64) );

	unsigned i,j,k;
	/// vertical broadcast bit
	__m256i x0br[128*3];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(x0,*(__m256i*)extract_01byte);
		x0 = _mm256_srli_si256( x0 , 2 );
		__m256i tt3 = _mm256_cmpeq_epi8( tt & (*(__m256i*)extract_bit) , zero );
		__m256i inp_run0 = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );

		tt = _mm256_shuffle_epi8(x1,*(__m256i*)extract_01byte);
		x1 = _mm256_srli_si256( x1 , 2 );
		tt3 = _mm256_cmpeq_epi8( tt & (*(__m256i*)extract_bit) , zero );
		__m256i inp_run1 = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );

		tt = _mm256_shuffle_epi8(x2,*(__m256i*)extract_01byte);
		x2 = _mm256_srli_si256( x2 , 2 );
		tt3 = _mm256_cmpeq_epi8( tt & (*(__m256i*)extract_bit) , zero );
		__m256i inp_run2 = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );

		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run0 , zero );
			inp_run0 = _mm256_srli_si256( inp_run0 , 1 );
			x0br[(i+j)*3] = br0;
			br0 = _mm256_shuffle_epi8( inp_run1 , zero );
			inp_run1 = _mm256_srli_si256( inp_run1 , 1 );
			x0br[(i+j)*3+1] = br0;
			br0 = _mm256_shuffle_epi8( inp_run2 , zero );
			inp_run2 = _mm256_srli_si256( inp_run2 , 1 );
			x0br[(i+j)*3+2] = br0;
		}
	}

	__m256i y0 = _mm256_load_si256( (__m256i const *) _y );
	__m256i y1 = _mm256_load_si256( (__m256i const *) (_y+32) );
	__m256i y2 = _mm256_load_si256( (__m256i const *) (_y+64) );

	/// vertical broadcast bit
	__m256i y0br[128*3];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(y0,*(__m256i*)extract_01byte);
		y0 = _mm256_srli_si256( y0 , 2 );
		__m256i tt3 = _mm256_cmpeq_epi8( tt & (*(__m256i*)extract_bit) , zero );
		__m256i inp_run0 = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );

		tt = _mm256_shuffle_epi8(y1,*(__m256i*)extract_01byte);
		y1 = _mm256_srli_si256( y1 , 2 );
		tt3 = _mm256_cmpeq_epi8( tt & (*(__m256i*)extract_bit) , zero );
		__m256i inp_run1 = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );

		tt = _mm256_shuffle_epi8(y2,*(__m256i*)extract_01byte);
		y2 = _mm256_srli_si256( y2 , 2 );
		tt3 = _mm256_cmpeq_epi8( tt & (*(__m256i*)extract_bit) , zero );
		__m256i inp_run2 = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );

		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run0 , zero );
			inp_run0 = _mm256_srli_si256( inp_run0 , 1 );
			y0br[(i+j)*3] = br0;
			br0 = _mm256_shuffle_epi8( inp_run1 , zero );
			inp_run1 = _mm256_srli_si256( inp_run1 , 1 );
			y0br[(i+j)*3+1] = br0;
			br0 = _mm256_shuffle_epi8( inp_run2 , zero );
			inp_run2 = _mm256_srli_si256( inp_run2 , 1 );
			y0br[(i+j)*3+2] = br0;
		}
	}

	x0 = _mm256_load_si256( (__m256i const *) w );
	y0 = _mm256_load_si256( (__m256i const *) _y );
	__m256i x_usd0 = _mm256_permute4x64_epi64( x0 , 0x4e );
	__m256i y_usd0 = _mm256_permute4x64_epi64( y0 , 0x4e );

	x1 = _mm256_load_si256( (__m256i const *) (w+32) );
	y1 = _mm256_load_si256( (__m256i const *) (_y+32) );
	__m256i x_usd1 = _mm256_permute4x64_epi64( x1 , 0x4e );
	__m256i y_usd1 = _mm256_permute4x64_epi64( y1 , 0x4e );

	x2 = _mm256_load_si256( (__m256i const *) (w+64) );
	y2 = _mm256_load_si256( (__m256i const *) (_y+64) );
	__m256i x_usd2 = _mm256_permute4x64_epi64( x2 , 0x4e );
	__m256i y_usd2 = _mm256_permute4x64_epi64( y2 , 0x4e );


	const __m256i *eqs = (const __m256i *) pk_mat;
	/// linear terms
	eqs += 128;

	__m256i mask_1 = _mm256_load_si256((__m256i*)_mask_1);
	/// quadratic terms
	eqs += 1;
	for(k=1;k<64;k++) {
		__m256i a0 = (x0 & y0br[k*3]) ^ (y0 & x0br[k*3]);
		__m256i a1 = (x_usd0 & y0br[k*3] ) ^ (y_usd0 & x0br[k*3]);
		__m256i quad_terms0 = _mm256_permute2x128_si256(a0,a1,0x21) ^ a0 ^ _mm256_permute2x128_si256(a1,a1,0x18);

		__m256i a2 = (x1 & y0br[k*3+1]) ^ (y1 & x0br[k*3+1]);
		__m256i a3 = (x_usd1 & y0br[k*3+1] ) ^ (y_usd1 & x0br[k*3+1]);
		__m256i quad_terms1 = _mm256_permute2x128_si256(a2,a3,0x21) ^ a2 ^ _mm256_permute2x128_si256(a3,a3,0x18);

		__m256i a4 = (x2 & y0br[k*3+2]) ^ (y2 & x0br[k*3+2]);
		__m256i a5 = (x_usd2 & y0br[k*3+2] ) ^ (y_usd2 & x0br[k*3+2]);
		__m256i quad_terms2 = _mm256_permute2x128_si256(a4,a5,0x21) ^ a4 ^ _mm256_permute2x128_si256(a5,a5,0x18);

		for(i=0;i<k;i++) {
			__m256i br0 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms0,zero));
			quad_terms0 = _mm256_srli_epi64( quad_terms0 , 1 );
			__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );

			__m256i br2 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms1,zero));
			quad_terms1 = _mm256_srli_epi64( quad_terms1 , 1 );
			__m256i br3 = _mm256_permute4x64_epi64( br2 , 0x4e );

			__m256i br4 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms2,zero));
			quad_terms2 = _mm256_srli_epi64( quad_terms2 , 1 );
			__m256i br5 = _mm256_permute4x64_epi64( br4 , 0x4e );

			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );

			r2 ^= ( br2 & (*eqs) );
			r3 ^= ( br3 & (*eqs) );

			r4 ^= ( br4 & (*eqs) );
			r5 ^= ( br5 & (*eqs) );
			eqs += 1;
		}
		eqs += 1;
	}

	for(k=64;k<128;k++) {
		__m256i a0 = (x0 & y0br[k*3]) ^ (y0 & x0br[k*3]);
		__m256i a1 = (x_usd0 & y0br[k*3] ) ^ (y_usd0 & x0br[k*3]);
		__m256i quad_terms0 = _mm256_permute2x128_si256(a0,a1,0x21) ^ a0 ^ _mm256_permute2x128_si256(a1,a1,0x18);

		__m256i a2 = (x1 & y0br[k*3+1]) ^ (y1 & x0br[k*3+1]);
		__m256i a3 = (x_usd1 & y0br[k*3+1] ) ^ (y_usd1 & x0br[k*3+1]);
		__m256i quad_terms1 = _mm256_permute2x128_si256(a2,a3,0x21) ^ a2 ^ _mm256_permute2x128_si256(a3,a3,0x18);

		__m256i a4 = (x2 & y0br[k*3+2]) ^ (y2 & x0br[k*3+2]);
		__m256i a5 = (x_usd2 & y0br[k*3+2] ) ^ (y_usd2 & x0br[k*3+2]);
		__m256i quad_terms2 = _mm256_permute2x128_si256(a4,a5,0x21) ^ a4 ^ _mm256_permute2x128_si256(a5,a5,0x18);

		__m256i quad_terms0_h64 = _mm256_srli_si256( quad_terms0 , 8 );
		__m256i quad_terms1_h64 = _mm256_srli_si256( quad_terms1 , 8 );
		__m256i quad_terms2_h64 = _mm256_srli_si256( quad_terms2 , 8 );

		for(i=0;i<64;i++) {
			__m256i br0 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms0,zero));
			quad_terms0 = _mm256_srli_epi64( quad_terms0 , 1 );
			__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			__m256i br2 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms1,zero));
			quad_terms1 = _mm256_srli_epi64( quad_terms1 , 1 );
			__m256i br3 = _mm256_permute4x64_epi64( br2 , 0x4e );
			__m256i br4 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms2,zero));
			quad_terms2 = _mm256_srli_epi64( quad_terms2 , 1 );
			__m256i br5 = _mm256_permute4x64_epi64( br4 , 0x4e );

			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			r2 ^= ( br2 & (*eqs) );
			r3 ^= ( br3 & (*eqs) );
			r4 ^= ( br4 & (*eqs) );
			r5 ^= ( br5 & (*eqs) );
			eqs += 1;
		}
		quad_terms0 = quad_terms0_h64;
		quad_terms1 = quad_terms1_h64;
		quad_terms2 = quad_terms2_h64;
		for(i=64;i<k;i++) {
			__m256i br0 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms0,zero));
			quad_terms0 = _mm256_srli_epi64( quad_terms0 , 1 );
			__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			__m256i br2 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms1,zero));
			quad_terms1 = _mm256_srli_epi64( quad_terms1 , 1 );
			__m256i br3 = _mm256_permute4x64_epi64( br2 , 0x4e );
			__m256i br4 = _mm256_cmpeq_epi8(mask_1,mask_1&_mm256_shuffle_epi8(quad_terms2,zero));
			quad_terms2 = _mm256_srli_epi64( quad_terms2 , 1 );
			__m256i br5 = _mm256_permute4x64_epi64( br4 , 0x4e );

			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			r2 ^= ( br2 & (*eqs) );
			r3 ^= ( br3 & (*eqs) );
			r4 ^= ( br4 & (*eqs) );
			r5 ^= ( br5 & (*eqs) );
			eqs += 1;
		}
		eqs += 1;
	}

	__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)z , rr );

	r0 = r2;
	r1 = r3;
	rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)(z+32) , rr );

	r0 = r4;
	r1 = r5;
	rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)(z+64) , rr );
}











static void _mq_gf4_n128_m128_vartime( uint8_t * z , const uint8_t * w , const uint8_t * pk_mat  )
{
	__m256i zero = _mm256_setzero_si256();
	__m256i r0 = _mm256_setzero_si256();
	__m256i r1 = _mm256_setzero_si256();

	__m256i inp = _mm256_load_si256( (__m256i const *) w );
	__m128i inp_and = _mm256_extracti128_si256( inp , 0 ) | _mm256_extracti128_si256( inp , 1 );
	uint64_t v0 = (uint64_t)_mm_extract_epi64( inp_and , 0 );
	uint64_t v1 = (uint64_t)_mm_extract_epi64( inp_and , 1 );
	unsigned i,j,k;
	/// vertical broadcast bit
	__m256i x0br[128];
	for(i=0;i<128;i+=16) {
		__m256i tt = _mm256_shuffle_epi8(inp,*(__m256i*)extract_01byte);
		inp = _mm256_srli_si256( inp , 2 );
		__m256i tt2 = tt & (*(__m256i*)extract_bit);
		__m256i tt3 = _mm256_cmpeq_epi8( tt2 , zero );
		__m256i inp_run = _mm256_andnot_si256( tt3 , _mm256_cmpeq_epi8(tt3,tt3) );
		for(j=0;j<16;j++) {
			__m256i br0 = _mm256_shuffle_epi8( inp_run , zero );
			inp_run = _mm256_srli_si256( inp_run , 1 );
			x0br[(i+j)] = br0;
		}
	}

	const __m256i *eqs = (const __m256i *) pk_mat;
	/// linear terms
	for(i=0;i<128;i++) {
		__m256i br0 = x0br[i];
		__m256i br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
		r0 ^= ( br0 & (*eqs) );
		r1 ^= ( br1 & (*eqs) );
		eqs += 1;
	}

	/// quadratic terms
	for(k=0;k<64;k++) {
		uint32_t jj = v0&1;
		v0 >>= 1;
		if( 0==jj ) {
			eqs += (k+1);
			continue;
		}
		__m256i br0;
		__m256i br1;
		__m256i t[2];
		t[0] = r0;
		t[1] = r1;
		r0 ^= r0;
		r1 ^= r1;
		for(i=0;i<=k;i++) {
			br0 = x0br[i];
			br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			eqs += 1;
		}
		__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
		r0 = (rr&br0);
		r1 = (rr&br1);
		r0 ^= t[0];
		r1 ^= t[1];
	}
	for(k=64;k<128;k++) {
		uint32_t jj = v1&1;
		v1 >>= 1;
		if( 0==jj ) {
			eqs += (k+1);
			continue;
		}
		__m256i br0;
		__m256i br1;
		__m256i t[2];
		t[0] = r0;
		t[1] = r1;
		r0 ^= r0;
		r1 ^= r1;
		for(i=0;i<=k;i++) {
			br0 = x0br[i];
			br1 = _mm256_permute4x64_epi64( br0 , 0x4e );
			r0 ^= ( br0 & (*eqs) );
			r1 ^= ( br1 & (*eqs) );
			eqs += 1;
		}
		__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
		r0 = (rr&br0);
		r1 = (rr&br1);
		r0 ^= t[0];
		r1 ^= t[1];
	}
	__m256i rr =  _mm256_permute2x128_si256(r0,r1,0x21) ^ r0 ^ _mm256_permute2x128_si256(r1,r1,0x18);
	_mm256_store_si256((__m256i*)z , rr );
}




void mq_gf4_n128_m128( uint64_t * fx , const uint64_t * x , const uint64_t * F )
{
	_mq_gf4_n128_m128((uint8_t*)fx,(const uint8_t *)x, (const uint8_t *)F );
}

void mq_gf4_n128_m128_x3( uint64_t * fx , const uint64_t * x , const uint64_t * F )
{
	_mq_gf4_n128_m128_x3((uint8_t*)fx,(const uint8_t *)x, (const uint8_t *)F );
}

void mq_gf4_n128_m128_x4( uint64_t * fx , const uint64_t * x , const uint64_t * F )
{
	_mq_gf4_n128_m128_x4((uint8_t*)fx,(const uint8_t *)x, (const uint8_t *)F );
}

void G_gf4_n128_m128( uint64_t * fx , const uint64_t * x , const uint64_t * y , const uint64_t * F  )
{
	_G_gf4_n128_m128((uint8_t*)fx,(const uint8_t *)x, (const uint8_t *)y, (const uint8_t *)F );
}

void G_gf4_n128_m128_x3( uint64_t * fx , const uint64_t * x , const uint64_t * y , const uint64_t * F  )
{
	_G_gf4_n128_m128_x3((uint8_t*)fx,(const uint8_t *)x, (const uint8_t *)y, (const uint8_t *)F );
}

void mq_gf4_n128_m128_vartime( uint64_t * fx , const uint64_t * x , const uint64_t * F )
{
	_mq_gf4_n128_m128_vartime((uint8_t*)fx,(const uint8_t *)x, (const uint8_t *)F );
}



#ifdef  __cplusplus
}
#endif
