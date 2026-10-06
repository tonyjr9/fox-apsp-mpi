#define CHECK(cond)                                                 \
  do {                                                              \
      tests_run++;                                                  \
      if (!(cond)) {                                                \
          tests_failed++;                                           \
          fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__,   \
                  #cond);                                           \
      }                                                             \
  } while (0)
