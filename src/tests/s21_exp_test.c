#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(exp_test_1) {
  ck_assert_ldouble_eq_tol(s21_exp(0.1), exp(0.1), 1e-6);
}
END_TEST

START_TEST(exp_test_2) {
  ck_assert_ldouble_eq_tol(s21_exp(1.0), exp(1.0), 1e-6);
}
END_TEST

START_TEST(exp_test_3) {
  ck_assert_ldouble_eq_tol(s21_exp(0.0004), exp(0.0004), 1e-6);
}
END_TEST

START_TEST(exp_test_4) {
  ck_assert_ldouble_eq_tol(s21_exp(-0.0004), exp(-0.0004), 1e-6);
}
END_TEST

START_TEST(exp_test_5) {
  ck_assert_ldouble_eq_tol(s21_exp(-24.4353425), exp(-24.4353425), 1e-6);
}
END_TEST

START_TEST(exp_test_6) {
  ck_assert_ldouble_eq_tol(s21_exp(0.0), exp(0.0), 1e-6);
}
END_TEST
START_TEST(test_range) {
  double a[7] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e10, 0.0, -564674};
  int i = _i;
  double res = exp(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_exp(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_exp(a[i]));
  else
    ck_assert_ldouble_eq_tol(s21_exp(a[i]), res, 1e-6);
}
END_TEST

START_TEST(test_range_0_1) {
  int i = _i;
  double x = i;
  ck_assert_ldouble_eq_tol(s21_exp(x), s21_exp(x), 1e-6);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Exp Test");
  TCase *tc1_1 = tcase_create("Exp Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, exp_test_1);
  tcase_add_test(tc1_1, exp_test_2);
  tcase_add_test(tc1_1, exp_test_3);
  tcase_add_test(tc1_1, exp_test_4);
  tcase_add_test(tc1_1, exp_test_5);
  tcase_add_test(tc1_1, exp_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 7);
  tcase_add_loop_test(tc1_1, test_range_0_1, -300, 300);
  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}