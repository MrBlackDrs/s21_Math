#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(sin_test_1) {
  ck_assert_ldouble_eq_tol(s21_sin(0.98439849385), sin(0.98439849385), 1e-6);
}
END_TEST

START_TEST(sin_test_2) { ck_assert_ldouble_eq_tol(s21_sin(PI), sin(PI), 1e-6); }
END_TEST

START_TEST(sin_test_3) {
  ck_assert_ldouble_eq_tol(s21_sin(-PI), sin(-PI), 1e-6);
}
END_TEST

START_TEST(sin_test_4) {
  ck_assert_ldouble_eq_tol(s21_sin(0.0), sin(0.0), 1e-6);
}
END_TEST

START_TEST(sin_test_5) {
  ck_assert_ldouble_eq_tol(s21_sin(PI / 2), sin(PI / 2), 1e-6);
}
END_TEST

START_TEST(sin_test_6) {
  ck_assert_ldouble_eq_tol(s21_sin(-PI / 2), sin(-PI / 2), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[7] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e10, 0.0, -564674};
  int i = _i;
  double res = sin(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_sin(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_sin(a[i]));
  else
    ck_assert_ldouble_eq_tol(s21_sin(a[i]), res, 1e-6);
}
END_TEST

START_TEST(test_range_0_2_PI) {
  int i = _i;
  double dx = 2 * PI / 360;
  ck_assert_ldouble_eq_tol(s21_sin(i * dx), sin(i * dx), 1e-6);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Sin Test");
  TCase *tc1_1 = tcase_create("Sin Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, sin_test_1);
  tcase_add_test(tc1_1, sin_test_2);
  tcase_add_test(tc1_1, sin_test_3);
  tcase_add_test(tc1_1, sin_test_4);
  tcase_add_test(tc1_1, sin_test_5);
  tcase_add_test(tc1_1, sin_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 7);
  tcase_add_loop_test(tc1_1, test_range_0_2_PI, 0, 360);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}