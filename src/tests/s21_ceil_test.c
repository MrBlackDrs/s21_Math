#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(ceil_test_1) {
  ck_assert_ldouble_eq_tol(s21_ceil(100.51), ceil(100.51), 1e-6);
}
END_TEST

START_TEST(ceil_test_2) {
  ck_assert_ldouble_eq_tol(s21_ceil(5), ceil(5), 1e-6);
}
END_TEST

START_TEST(ceil_test_3) {
  ck_assert_ldouble_eq_tol(s21_ceil(10000.2), ceil(10000.2), 1e-6);
}
END_TEST

START_TEST(ceil_test_4) {
  ck_assert_ldouble_eq_tol(s21_ceil(-3.3), ceil(-3.3), 1e-6);
}
END_TEST

START_TEST(ceil_test_5) {
  ck_assert_ldouble_eq_tol(s21_ceil(-4732323456789076.5),
                           ceil(-4732323456789076.5), 1e-6);
}
END_TEST

START_TEST(ceil_test_6) {
  ck_assert_ldouble_eq_tol(s21_ceil(4.7323234567890765),
                           ceil(4.7323234567890765), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[9] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e308,
                 0.0,  -1e308,  1e-18,    -1e-18};
  int i = _i;
  double res = ceil(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_ceil(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_ceil(a[i]));
  else
    ck_assert_ldouble_eq(s21_ceil(a[i]), res);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Ceil Test");
  TCase *tc1_1 = tcase_create("Ceil Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, ceil_test_1);
  tcase_add_test(tc1_1, ceil_test_2);
  tcase_add_test(tc1_1, ceil_test_3);
  tcase_add_test(tc1_1, ceil_test_4);
  tcase_add_test(tc1_1, ceil_test_5);
  tcase_add_test(tc1_1, ceil_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 9);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}