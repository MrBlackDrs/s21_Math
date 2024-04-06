#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(fmod_test_1) {
  ck_assert_ldouble_eq_tol(s21_fmod(100.51, 9.0401), fmod(100.51, 9.0401),
                           1e-6);
}
END_TEST

START_TEST(fmod_test_2) { ck_assert_double_nan(s21_fmod(5, 0.0)); }
END_TEST

START_TEST(fmod_test_3) {
  ck_assert_ldouble_eq_tol(s21_fmod(-10, 2), fmod(-10, 2), 1e-6);
}
END_TEST

START_TEST(fmod_test_4) {
  ck_assert_ldouble_eq_tol(s21_fmod(-25, -3.3), fmod(-25, -3.3), 1e-6);
}
END_TEST

START_TEST(fmod_test_5) {
  ck_assert_ldouble_eq_tol(s21_fmod(50, -10), fmod(50, -10), 1e-6);
}
END_TEST

START_TEST(fmod_test_6) {
  ck_assert_ldouble_eq_tol(s21_fmod(0, 5), fmod(0, 5), 1e-6);
}
END_TEST

START_TEST(test_range) {
  double a[7] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e308, 0.0, -1e308};
  int i = _i;
  for (int j = 0; j < 7; j++) {
    double res = fmod(a[i], a[j]);
    if (s21_isnan(res))
      ck_assert_ldouble_nan(s21_fmod(a[i], a[j]));
    else if (s21_isinf(res))
      ck_assert_ldouble_infinite(s21_fmod(a[i], a[j]));
    else
      ck_assert_ldouble_eq_tol(s21_fmod(a[i], a[j]), res, 1e-6);
  }
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Fmod Test");
  TCase *tc1_1 = tcase_create("Fmod Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, fmod_test_1);
  tcase_add_test(tc1_1, fmod_test_2);
  tcase_add_test(tc1_1, fmod_test_3);
  tcase_add_test(tc1_1, fmod_test_4);
  tcase_add_test(tc1_1, fmod_test_5);
  tcase_add_test(tc1_1, fmod_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 7);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}