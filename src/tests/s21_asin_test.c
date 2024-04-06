#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(asin_test_1) {
  ck_assert_ldouble_eq_tol(s21_asin(1.0), asin(1.0), 1e-6);
}
END_TEST

START_TEST(asin_test_2) {
  ck_assert_ldouble_eq_tol(s21_asin(-1.0), asin(-1.0), 1e-6);
}
END_TEST

START_TEST(asin_test_3) {
  ck_assert_ldouble_eq_tol(s21_asin(0), asin(0), 1e-6);
}
END_TEST

START_TEST(asin_test_4) {
  ck_assert_ldouble_eq_tol(s21_asin(0.23234), asin(0.23234), 1e-6);
}
END_TEST

START_TEST(asin_test_5) {
  ck_assert_ldouble_eq_tol(s21_asin(-0.23234), asin(-0.23234), 1e-6);
}
END_TEST

START_TEST(asin_test_6) {
  ck_assert_ldouble_eq_tol(s21_asin(0.0002), asin(0.0002), 1e-6);
}
END_TEST

START_TEST(asin_test_7) {
  ck_assert_ldouble_eq_tol(s21_asin(-0.0002), asin(-0.0002), 1e-6);
}
END_TEST

START_TEST(asin_test_8) {
  ck_assert_ldouble_eq_tol(s21_asin(2.2250738585072014e-308),
                           asin(2.2250738585072014e-308), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[7] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e10, 0.0, -564674};
  int i = _i;
  double res = asin(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_asin(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_asin(a[i]));
  else
    ck_assert_ldouble_eq_tol(s21_asin(a[i]), res, 1e-6);
}
END_TEST

START_TEST(test_range_0_1) {
  int i = _i;
  double x = i / 300;
  ck_assert_ldouble_eq_tol(s21_asin(x), asin(x), 1e-6);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Asin Test");
  TCase *tc1_1 = tcase_create("Asin Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, asin_test_1);
  tcase_add_test(tc1_1, asin_test_2);
  tcase_add_test(tc1_1, asin_test_3);
  tcase_add_test(tc1_1, asin_test_4);
  tcase_add_test(tc1_1, asin_test_5);
  tcase_add_test(tc1_1, asin_test_6);
  tcase_add_test(tc1_1, asin_test_7);
  tcase_add_test(tc1_1, asin_test_8);
  tcase_add_loop_test(tc1_1, test_range, 0, 7);
  tcase_add_loop_test(tc1_1, test_range_0_1, 0, 300);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}