#include <stdint.h>



#ifdef  __cplusplus
extern  "C" {
#endif




void mq_gf4_n128_m128_ref( uint64_t * z , const uint64_t * w , const uint64_t * pk_mat  )
{
	uint64_t r0_0 = 0;
	uint64_t r0_1 = 0;
	uint64_t r1_0 = 0;
	uint64_t r1_1 = 0;
	uint64_t r2_0 = 0;
	uint64_t r2_1 = 0;

	uint64_t x0_0 = w[0];
	uint64_t x1_0 = w[1];
	uint64_t x0_1 = w[2];
	uint64_t x1_1 = w[3];

	uint64_t x_0[128];
	uint64_t x_1[128];
	unsigned i,k;
	/* broadcast bit */
	for(i=0;i<64;i++) {
		x_0[i] = 0-(x0_0&1);
		x_1[i] = 0-(x0_1&1);
		x0_0 >>= 1;
		x0_1 >>= 1;
		x_0[i+64] = 0-(x1_0&1);
		x_1[i+64] = 0-(x1_1&1);
		x1_0 >>= 1;
		x1_1 >>= 1;
	}

	/* linear terms */
	for(i=0;i<128;i++) {
		uint64_t br0 = x_0[i];
		uint64_t br1 = x_1[i];

		uint64_t eq0_0 = pk_mat[0];
		uint64_t eq0_1 = pk_mat[1];
		uint64_t eq1_0 = pk_mat[2];
		uint64_t eq1_1 = pk_mat[3];

		r0_0 ^= (br0 & eq0_0 );
		r0_1 ^= (br0 & eq0_1 );
		r1_0 ^= (br0 & eq1_0 );
		r1_1 ^= (br0 & eq1_1 );

		r1_0 ^= (br1 & eq0_0 );
		r1_1 ^= (br1 & eq0_1 );
		r2_0 ^= (br1 & eq1_0 );
		r2_1 ^= (br1 & eq1_1 );
		pk_mat += 4;
	}

	/* quadratic terms */
	for(k=0;k<128;k++) {
		uint64_t br0;
		uint64_t br1;
		uint64_t t0_0 = 0;
		uint64_t t0_1 = 0;
		uint64_t t1_0 = 0;
		uint64_t t1_1 = 0;
		uint64_t t2_0 = 0;
		uint64_t t2_1 = 0;
		for(i=0;i<=k;i++) {
			br0 = x_0[i];
			br1 = x_1[i];

			uint64_t eq0_0 = pk_mat[0];
			uint64_t eq0_1 = pk_mat[1];
			uint64_t eq1_0 = pk_mat[2];
			uint64_t eq1_1 = pk_mat[3];

			t0_0 ^= (br0 & eq0_0 );
			t0_1 ^= (br0 & eq0_1 );
			t1_0 ^= (br0 & eq1_0 );
			t1_1 ^= (br0 & eq1_1 );

			t1_0 ^= (br1 & eq0_0 );
			t1_1 ^= (br1 & eq0_1 );
			t2_0 ^= (br1 & eq1_0 );
			t2_1 ^= (br1 & eq1_1 );
			pk_mat += 4;
		}
		t0_0 ^= t2_0;
		t0_1 ^= t2_1;
		t1_0 ^= t2_0;
		t1_1 ^= t2_1;
		r0_0 ^= (t0_0 & br0);
		r0_1 ^= (t0_1 & br0);
		r1_0 ^= (t1_0 & br0);
		r1_1 ^= (t1_1 & br0);
		r1_0 ^= (t0_0 & br1);
		r1_1 ^= (t0_1 & br1);
		r2_0 ^= (t1_0 & br1);
		r2_1 ^= (t1_1 & br1);
	}

	r0_0 ^= r2_0;
	r0_1 ^= r2_1;
	r1_0 ^= r2_0;
	r1_1 ^= r2_1;

	z[0] = r0_0;
	z[1] = r0_1;
	z[2] = r1_0;
	z[3] = r1_1;

}


void G_gf4_n128_m128_ref( uint64_t * z , const uint64_t * w , const uint64_t * _y , const uint64_t * pk_mat )
{

	uint64_t r0_0 = 0;
	uint64_t r0_1 = 0;
	uint64_t r1_0 = 0;
	uint64_t r1_1 = 0;
	uint64_t r2_0 = 0;
	uint64_t r2_1 = 0;

	uint64_t x0_0 = w[0];
	uint64_t x1_0 = w[1];
	uint64_t x0_1 = w[2];
	uint64_t x1_1 = w[3];
	uint64_t y0_0 = _y[0];
	uint64_t y1_0 = _y[1];
	uint64_t y0_1 = _y[2];
	uint64_t y1_1 = _y[3];

	uint64_t x_0[128];
	uint64_t x_1[128];
	uint64_t y_0[128];
	uint64_t y_1[128];
	unsigned i,k;
	/* broadcast bit */
	for(i=0;i<64;i++) {
		x_0[i] = 0-(x0_0&1);
		x_1[i] = 0-(x0_1&1);
		x0_0 >>= 1;
		x0_1 >>= 1;
		x_0[i+64] = 0-(x1_0&1);
		x_1[i+64] = 0-(x1_1&1);
		x1_0 >>= 1;
		x1_1 >>= 1;

		y_0[i] = 0-(y0_0&1);
		y_1[i] = 0-(y0_1&1);
		y0_0 >>= 1;
		y0_1 >>= 1;
		y_0[i+64] = 0-(y1_0&1);
		y_1[i+64] = 0-(y1_1&1);
		y1_0 >>= 1;
		y1_1 >>= 1;

	}

	/* linear terms */
	pk_mat += 4*128;

	/* quadratic terms */
	pk_mat += 4;
	for(k=1;k<128;k++) {
		for(i=0;i<k;i++) {
			uint64_t br0 = (x_0[i]&y_0[k])^(x_1[i]&y_1[k])^ (x_0[k]&y_0[i])^(x_1[k]&y_1[i]);
			uint64_t br1 = (x_1[i]&y_1[k])^(x_0[i]&y_1[k])^(x_1[i]&y_0[k])^ (x_1[k]&y_1[i])^(x_0[k]&y_1[i])^(x_1[k]&y_0[i]);

			uint64_t eq0_0 = pk_mat[0];
			uint64_t eq0_1 = pk_mat[1];
			uint64_t eq1_0 = pk_mat[2];
			uint64_t eq1_1 = pk_mat[3];

			r0_0 ^= (br0 & eq0_0 );
			r0_1 ^= (br0 & eq0_1 );
			r1_0 ^= (br0 & eq1_0 );
			r1_1 ^= (br0 & eq1_1 );

			r1_0 ^= (br1 & eq0_0 );
			r1_1 ^= (br1 & eq0_1 );
			r2_0 ^= (br1 & eq1_0 );
			r2_1 ^= (br1 & eq1_1 );
			pk_mat += 4;
		}
		pk_mat += 4;
	}
	r0_0 ^= r2_0;
	r0_1 ^= r2_1;
	r1_0 ^= r2_0;
	r1_1 ^= r2_1;

	z[0] = r0_0;
	z[1] = r0_1;
	z[2] = r1_0;
	z[3] = r1_1;

}




#ifdef  __cplusplus
}
#endif
