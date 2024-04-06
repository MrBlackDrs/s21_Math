#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(tan_test_1) {
  ck_assert_ldouble_eq_tol(s21_tan(0.98439849385), tan(0.98439849385), 1e-6);
}
END_TEST

START_TEST(tan_test_2) {
  ck_assert_ldouble_eq_tol(s21_tan(PI / 3), tan(PI / 3), 1e-6);
}
END_TEST

START_TEST(tan_test_3) {
  ck_assert_ldouble_eq_tol(s21_tan(-PI / 3), tan(-PI / 3), 1e-6);
}
END_TEST

START_TEST(tan_test_4) {
  ck_assert_ldouble_eq_tol(s21_tan(0.0), tan(0.0), 1e-6);
}
END_TEST

START_TEST(tan_test_5) {
  ck_assert_ldouble_eq_tol(s21_tan(0.0005), tan(0.0005), 1e-6);
}
END_TEST

START_TEST(tan_test_6) {
  ck_assert_ldouble_eq_tol(s21_tan(-0.0005), tan(-0.0005), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[7] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e10, 0.0, -564674};
  int i = _i;
  double res = tan(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_tan(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_tan(a[i]));
  else
    ck_assert_ldouble_eq_tol(s21_tan(a[i]), res, 1e-6);
}
END_TEST

START_TEST(test_range_0_2_PI) {
  int i = _i;
  double dx = 2 * PI / 123;
  ck_assert_ldouble_eq_tol(s21_tan(i * dx), tan(i * dx), 1e-6);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Tan Test");
  TCase *tc1_1 = tcase_create("Tan Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, tan_test_1);
  tcase_add_test(tc1_1, tan_test_2);
  tcase_add_test(tc1_1, tan_test_3);
  tcase_add_test(tc1_1, tan_test_4);
  tcase_add_test(tc1_1, tan_test_5);
  tcase_add_test(tc1_1, tan_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 7);
  tcase_add_loop_test(tc1_1, test_range_0_2_PI, 0, 123);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}