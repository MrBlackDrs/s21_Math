#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(cos_test_1) {
  ck_assert_ldouble_eq_tol(s21_cos(0.98439849385), cos(0.98439849385), 1e-6);
}
END_TEST

START_TEST(cos_test_2) { ck_assert_ldouble_eq_tol(s21_cos(PI), cos(PI), 1e-6); }
END_TEST

START_TEST(cos_test_3) {
  ck_assert_ldouble_eq_tol(s21_cos(-PI), cos(-PI), 1e-6);
}
END_TEST

START_TEST(cos_test_4) {
  ck_assert_ldouble_eq_tol(s21_cos(0.0), cos(0.0), 1e-6);
}
END_TEST

START_TEST(cos_test_5) {
  ck_assert_ldouble_eq_tol(s21_cos(PI / 2), cos(PI / 2), 1e-6);
}
END_TEST

START_TEST(cos_test_6) {
  ck_assert_ldouble_eq_tol(s21_cos(-PI / 2), cos(-PI / 2), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[7] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e10, 0.0, -564674};
  int i = _i;
  double res = cos(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_cos(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_cos(a[i]));
  else
    ck_assert_ldouble_eq_tol(s21_cos(a[i]), res, 1e-6);
}
END_TEST

START_TEST(test_range_0_2_PI) {
  int i = _i;
  double dx = 2 * PI / 360;
  ck_assert_ldouble_eq_tol(s21_cos(i * dx), cos(i * dx), 1e-6);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Cos Test");
  TCase *tc1_1 = tcase_create("Cos Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, cos_test_1);
  tcase_add_test(tc1_1, cos_test_2);
  tcase_add_test(tc1_1, cos_test_3);
  tcase_add_test(tc1_1, cos_test_4);
  tcase_add_test(tc1_1, cos_test_5);
  tcase_add_test(tc1_1, cos_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 7);
  tcase_add_loop_test(tc1_1, test_range_0_2_PI, 0, 360);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}