CC = /usr/bin/cc
CFLAGS = -O3 -Wall -Wextra -Wpedantic -march=native -std=c99

SOURCES = mq.c F.c randombytes.c sofia.c sampling.c fips202.c
HEADERS = mq.h F.h randombytes.h sofia.h sampling.h fips202.h params.h api.h

TESTS = test/test_mq \
		test/test_mq_consistency \
		test/test_3way_mq \
		test/test_sofia \
		test/test_sample_challenges \
		test/speed \

all: $(TESTS)

test/%: test/%.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $<

.PHONY: clean test

.DELETE_ON_ERROR:

test: $(TESTS)
	for test in $(TESTS); do\
		echo "=== $$test ==="; \
		./$$test; \
	done

clean:
	-$(RM) -r $(TESTS)
	-$(RM) test/keypair
	-$(RM) test/sign
	-$(RM) test/open
