#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(sqrt_test_1) {
  ck_assert_ldouble_eq_tol(s21_sqrt(5.2), sqrt(5.2), TEST_EPS);
}
END_TEST

START_TEST(sqrt_test_2) { ck_assert_double_nan(s21_sqrt(-5.2)); }
END_TEST

START_TEST(sqrt_test_3) {
  ck_assert_ldouble_eq_tol(s21_sqrt(0), sqrt(0), TEST_EPS);
}
END_TEST

START_TEST(sqrt_test_4) {
  ck_assert_ldouble_eq_tol(s21_sqrt(0.01), sqrt(0.01), TEST_EPS);
}
END_TEST

START_TEST(sqrt_test_5) {
  ck_assert_ldouble_eq_tol(s21_sqrt(0.123456), sqrt(0.123456), TEST_EPS);
}
END_TEST

START_TEST(sqrt_test_6) {
  ck_assert_ldouble_eq_tol(s21_sqrt(987565.9874641213656),
                           sqrt(987565.9874641213656), TEST_EPS);
}
END_TEST

START_TEST(sqrt_test_7) {
  ck_assert_ldouble_eq(s21_sqrt(S21_INF), sqrt(S21_INF));
}
END_TEST

START_TEST(sqrt_test_8) {
  ck_assert_ldouble_eq_tol(s21_sqrt(4.163435e-10), sqrt(4.163435e-10),
                           TEST_EPS);
}
END_TEST

START_TEST(sqrt_test_9) {
  ck_assert_ldouble_eq_tol(s21_tolerance(s21_sqrt(1.79769e+25)),
                           s21_tolerance(sqrt(1.79769e+25)), 1);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Sqrt Test");
  TCase *tc1_1 = tcase_create("Sqrt Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, sqrt_test_1);
  tcase_add_test(tc1_1, sqrt_test_2);
  tcase_add_test(tc1_1, sqrt_test_3);
  tcase_add_test(tc1_1, sqrt_test_4);
  tcase_add_test(tc1_1, sqrt_test_5);
  tcase_add_test(tc1_1, sqrt_test_6);
  tcase_add_test(tc1_1, sqrt_test_7);
  tcase_add_test(tc1_1, sqrt_test_8);
  tcase_add_test(tc1_1, sqrt_test_9);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}