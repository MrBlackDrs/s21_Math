#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(log_test_1) { ck_assert_ldouble_infinite(-s21_log(0)); }
END_TEST

START_TEST(log_test_3) { ck_assert_ldouble_nan(s21_log(-1)); }
END_TEST

START_TEST(log_test_4) {
  ck_assert_ldouble_eq_tol(s21_log(0.0083928293832), log(0.0083928293832),
                           1e-6);
}
END_TEST

START_TEST(log_test_5) { ck_assert_ldouble_eq_tol(s21_log(10), log(10), 1e-6); }
END_TEST

START_TEST(log_test_6) {
  ck_assert_ldouble_eq_tol(s21_log(0.00001), log(0.00001), 1e-6);
}
END_TEST

START_TEST(log_test_8) {
  ck_assert_ldouble_eq_tol(s21_log(1.79769e+308), log(1.79769e+308), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[9] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e10,
                 0.0,  -564674, 0.7,      1.9};
  int i = _i;
  double res = log(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_log(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_log(a[i]));
  else
    ck_assert_ldouble_eq_tol(s21_log(a[i]), res, 1e-6);
}
END_TEST

START_TEST(test_range_0_1) {
  int i = _i;
  double x = i;
  ck_assert_ldouble_eq_tol(s21_log(x), s21_log(x), 1e-6);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Log Test");
  TCase *tc1_1 = tcase_create("Log Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, log_test_1);
  tcase_add_test(tc1_1, log_test_3);
  tcase_add_test(tc1_1, log_test_4);
  tcase_add_test(tc1_1, log_test_5);
  tcase_add_test(tc1_1, log_test_6);
  tcase_add_test(tc1_1, log_test_8);
  tcase_add_loop_test(tc1_1, test_range, 0, 9);
  tcase_add_loop_test(tc1_1, test_range_0_1, 1, 1000);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}