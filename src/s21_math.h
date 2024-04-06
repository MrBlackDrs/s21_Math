#ifndef MATH
#define MATH
#endif

#include <stdio.h>
#include <stdlib.h>
#define PI 3.141592653589793

#define TEST_EPS 1E-6
#define EPS 1E-16
#define S21_INF __builtin_inf()
#define S21_NAN __builtin_nan("")

int s21_abs(int x);
// вычисляет абсолютное значение целого числа

long double s21_acos(double x);
// вычисляет арккосинус

long double s21_asin(double x);
// вычисляет арксинус

long double s21_atan(double x);
// вычисляет арктангенс

long double s21_ceil(double x);
// возвращает ближайшее целое число, не меньшее заданного значения

long double s21_cos(double x);
// вычисляет косинус

long double s21_exp(double x);
// возвращает значение e, возведенное в заданную степень

long double s21_fabs(double x);
// вычисляет абсолютное значение числа с плавающей точкой

long double s21_floor(double x);
// возвращает ближайшее целое число, не превышающее заданное значение

long double s21_fmod(double x, double y);
// остаток операции деления с плавающей точкой

long double s21_log(double x);
// вычисляет натуральный логарифм

long double s21_pow(double base, double exp);
// возводит число в заданную степень

long double s21_sin(double x);
// вычисляет синус

long double s21_sqrt(double x);
// вычисляет квадратный корень

long double s21_tan(double x);
// вычисляет тангенс

int s21_isnan(long double x);
// проверяет число на nan

int s21_isinf(long double x);
// проверяет число на inf

long double s21_tolerance(long double x);