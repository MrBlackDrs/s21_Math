#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(pow_test_1) {
  ck_assert_ldouble_eq_tol(s21_pow(5.2, 3.3), pow(5.2, 3.3), TEST_EPS);
}
END_TEST

START_TEST(pow_test_2) { ck_assert_ldouble_nan(s21_pow(-5.2, 3.3)); }
END_TEST

START_TEST(pow_test_3) {
  ck_assert_ldouble_eq_tol(s21_pow(5.9, -3.8), pow(5.9, -3.8), TEST_EPS);
}
END_TEST

START_TEST(pow_test_4) {
  ck_assert_ldouble_eq_tol(s21_pow(0, 0), pow(0, 0), TEST_EPS);
}
END_TEST

START_TEST(pow_test_5) {
  ck_assert_ldouble_eq_tol(s21_pow(0, 5.7), pow(0, 5.7), TEST_EPS);
}
END_TEST

START_TEST(pow_test_6) {
  ck_assert_ldouble_eq_tol(s21_pow(2989839, 0.0004), pow(2989839, 0.0004),
                           TEST_EPS);
}
END_TEST

START_TEST(pow_test_7) {
  ck_assert_ldouble_eq_tol(s21_pow(2989839, 0.0004), pow(2989839, 0.0004),
                           TEST_EPS);
}
END_TEST

START_TEST(pow_test_8) {
  ck_assert_ldouble_eq_tol(s21_pow(2989839.12, 0), pow(2989839.12, 0),
                           TEST_EPS);
}
END_TEST

START_TEST(pow_test_9) {
  ck_assert_ldouble_eq_tol(s21_tolerance(s21_pow(0.083928293832, -7.9)),
                           s21_tolerance(pow(0.083928293832, -7.9)), 1);
}
END_TEST

START_TEST(pow_test_10) {
  ck_assert_ldouble_eq(s21_pow(1.79769e+308, 3), pow(1.79769e+308, 3));
}
END_TEST

START_TEST(pow_test_11) {
  ck_assert_ldouble_eq_tol(s21_pow(-8, 3), pow(-8, 3), TEST_EPS);
}
END_TEST

START_TEST(test_range) {
  double a[7] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e308, 0.0, -1e308};
  int i = _i;
  for (int j = 0; j < 7; j++) {
    double res = pow(a[i], a[j]);
    if (s21_isnan(res))
      ck_assert_ldouble_nan(s21_pow(a[i], a[j]));
    else if (s21_isinf(res))
      ck_assert_ldouble_infinite(s21_pow(a[i], a[j]));
    else
      ck_assert_ldouble_eq_tol(s21_pow(a[i], a[j]), res, 1e-6);
  }
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Pow Test");
  TCase *tc1_1 = tcase_create("Pow Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, pow_test_1);
  tcase_add_test(tc1_1, pow_test_2);
  tcase_add_test(tc1_1, pow_test_3);
  tcase_add_test(tc1_1, pow_test_4);
  tcase_add_test(tc1_1, pow_test_5);
  tcase_add_test(tc1_1, pow_test_6);
  tcase_add_test(tc1_1, pow_test_7);
  tcase_add_test(tc1_1, pow_test_8);
  tcase_add_test(tc1_1, pow_test_9);
  tcase_add_test(tc1_1, pow_test_10);
  tcase_add_test(tc1_1, pow_test_11);
  tcase_add_loop_test(tc1_1, test_range, 0, 7);
  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}