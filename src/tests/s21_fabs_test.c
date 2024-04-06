#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(fabs_test_1) {
  ck_assert_ldouble_eq(s21_fabs(-100000.001421412), fabs(-100000.001421412));
}
END_TEST

START_TEST(fabs_test_2) { ck_assert_ldouble_eq(s21_fabs(1.0), fabs(1.0)); }
END_TEST

START_TEST(fabs_test_3) { ck_assert_ldouble_eq(s21_fabs(0.0), fabs(0.0)); }
END_TEST

START_TEST(fabs_test_4) {
  ck_assert_ldouble_eq(s21_fabs(-2147483647.0), fabs(-2147483647.0));
}
END_TEST

START_TEST(fabs_test_5) { ck_assert_ldouble_eq(s21_fabs(-0.5), fabs(-0.5)); }
END_TEST

START_TEST(fabs_test_6) { ck_assert_ldouble_eq(s21_fabs(-0.0), fabs(-0.0)); }
END_TEST

START_TEST(test_range) {
  double a[9] = {-0.0, S21_INF, -S21_INF, S21_NAN, 1e308,
                 0.0,  -1e308,  1e-18,    -1e-18};
  int i = _i;
  double res = fabs(a[i]);
  if (s21_isnan(res))
    ck_assert_ldouble_nan(s21_fabs(a[i]));
  else if (s21_isinf(res))
    ck_assert_ldouble_infinite(s21_fabs(a[i]));
  else
    ck_assert_ldouble_eq(s21_fabs(a[i]), res);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Fabs Test");
  TCase *tc1_1 = tcase_create("Fabs Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, fabs_test_1);
  tcase_add_test(tc1_1, fabs_test_2);
  tcase_add_test(tc1_1, fabs_test_3);
  tcase_add_test(tc1_1, fabs_test_4);
  tcase_add_test(tc1_1, fabs_test_5);
  tcase_add_test(tc1_1, fabs_test_6);
  tcase_add_loop_test(tc1_1, test_range, 0, 9);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}