#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(atan_test_1) {
  ck_assert_ldouble_eq_tol(s21_atan(0.1), atan(0.1), 1e-6);
}
END_TEST

START_TEST(atan_test_2) {
  ck_assert_ldouble_eq_tol(s21_atan(1.0), atan(1.0), 1e-6);
}
END_TEST

START_TEST(atan_test_3) {
  ck_assert_ldouble_eq_tol(s21_atan(1000.0), atan(1000.0), 1e-6);
}
END_TEST

START_TEST(atan_test_4) {
  ck_assert_ldouble_eq_tol(s21_atan(-1000.0), atan(-1000.0), 1e-6);
}
END_TEST

START_TEST(atan_test_5) {
  ck_assert_ldouble_eq_tol(s21_atan(0.0), atan(0.0), 1e-6);
}
END_TEST

START_TEST(atan_test_6) {
  ck_assert_ldouble_eq_tol(s21_atan(-0.0005), atan(-0.0005), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[11] = {-0.0,    S21_INF,  -S21_INF, S21_NAN,  1e10, 0.0,
                  -564674, 11. / 16, 19. / 16, 39. / 16, 1.5};
  int i = _i;
  double res = atan(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_atan(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_atan(a[i]));
  else
    ck_assert_ldouble_eq_tol(s21_atan(a[i]), res, 1e-6);
}
END_TEST

START_TEST(test_range_0_1) {
  int i = _i;
  double x = i;
  ck_assert_ldouble_eq_tol(s21_atan(x), atan(x), 1e-6);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Atan Test");
  TCase *tc1_1 = tcase_create("Atan Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, atan_test_1);
  tcase_add_test(tc1_1, atan_test_2);
  tcase_add_test(tc1_1, atan_test_3);
  tcase_add_test(tc1_1, atan_test_4);
  tcase_add_test(tc1_1, atan_test_5);
  tcase_add_test(tc1_1, atan_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 11);
  tcase_add_loop_test(tc1_1, test_range_0_1, -300, 300);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}