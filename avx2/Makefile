CC = /usr/bin/cc
GCC = /usr/bin/gcc
CFLAGS = -O3 -Wall -Wextra -Wpedantic -march=native -std=c99 -DSOFIA_GAP_F

SOURCES =      F.c randombytes.c sofia.c sampling.c fips202.c fips202x4.c prg.c
HEADERS = mq.h F.h randombytes.h sofia.h sampling.h fips202.h fips202x4.h crypto_stream_aes256ctr.h prg.h params.h api.h

OBJS = mq.s g.s mq3.s g3.s crypto_stream_aes256ctr.o

TESTS = test/test_mq \
		test/test_mq_consistency \
		test/test_3way_mq \
		test/test_F_gap \
		test/test_fips202x4 \
		test/test_sofia \
		test/test_sample_challenges \
		test/speed \

all: $(TESTS)

%.s: %.s.py
	python3 $^ > $@

test/%: test/%.c $(SOURCES) $(OBJS) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(OBJS) $<

test/speed: test/speed.c $(SOURCES) $(OBJS) $(HEADERS) mq_gf4_vertical.h
	$(CC) $(CFLAGS) -DSOFIA_MEASURE_ALTMQ4 -DSOFIA_MEASURE_4X_PERM_COMMIT -DSOFIA_EXPOSE_STATICS -o $@ $(SOURCES) mq_gf4_vertical_avx2.c $(OBJS) $<

# prevent strict aliasing issues/warnings
crypto_stream_aes256ctr.o: crypto_stream_aes256ctr.c
	$(GCC) $(CFLAGS) -fno-strict-aliasing -c $^

.PHONY: clean test

.DELETE_ON_ERROR:

test: $(TESTS)
	for test in $(TESTS); do\
		echo "=== $$test ==="; \
		./$$test; \
	done

clean:
	-find . -name '*.pyc' -delete
	-find . -name '__pycache__' -delete
	-$(RM) $(OBJS)
	-$(RM) $(TESTS)
	-$(RM) test/keypair
	-$(RM) test/sign
	-$(RM) test/open
