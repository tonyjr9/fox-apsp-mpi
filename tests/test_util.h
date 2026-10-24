#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include <stdio.h>

#include "utils.h"

static u32 tests_run;
static u32 tests_failed;

#define CHECK(cond)                                                 \
  do {                                                              \
      tests_run++;                                                  \
      if (!(cond)) {                                                \
          tests_failed++;                                           \
          fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__,   \
                  #cond);                                           \
      }                                                             \
  } while (0)

static i32 test_report(const char *name) {
  if (tests_failed == 0) {
    printf("PASS %s (%u checks)\n", name, tests_run);
    return 0;
  }

  fprintf(stderr, "FAIL %s (%u/%u checks failed)\n",
          name, tests_failed, tests_run);
  return 1;
}

#endif
