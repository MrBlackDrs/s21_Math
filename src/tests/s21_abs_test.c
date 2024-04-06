#include <check.h>
#include <math.h>

#include "../s21_math.h"

START_TEST(abs_test_1) { ck_assert_int_eq(s21_abs(-10000), abs(-10000)); }
END_TEST

START_TEST(abs_test_2) { ck_assert_int_eq(s21_abs(10000), abs(10000)); }
END_TEST

START_TEST(abs_test_3) { ck_assert_int_eq(s21_abs(0), abs(0)); }
END_TEST

START_TEST(abs_test_4) {
  ck_assert_int_eq(s21_abs(-2147483647), abs(-2147483647));
}
END_TEST

START_TEST(abs_test_5) { ck_assert_int_eq(s21_abs(-05), abs(-05)); }
END_TEST

START_TEST(abs_test_6) { ck_assert_int_eq(s21_abs(-0), abs(-0)); }
END_TEST

int main(void) {
  Suite *s1 = suite_create("Abs Test");
  TCase *tc1_1 = tcase_create("Abs Test");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, abs_test_1);
  tcase_add_test(tc1_1, abs_test_2);
  tcase_add_test(tc1_1, abs_test_3);
  tcase_add_test(tc1_1, abs_test_4);
  tcase_add_test(tc1_1, abs_test_5);
  tcase_add_test(tc1_1, abs_test_6);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}