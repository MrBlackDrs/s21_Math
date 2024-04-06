#include "s21_math.h"

int s21_abs(int x) {
  int res = x;
  if (x < 0) res = -x;
  return res;
}

long double s21_acos(double x) { return PI / 2 - s21_asin(x); }
// вычисляет арккосинус

long double s21_asin(double x) {
  long double res = 0;
  if (x < -1 || x > 1 || s21_isnan(x))
    res = S21_NAN;
  else if (x == 1)
    res = PI / 2;
  else if (x == -1)
    res = -PI / 2;
  else
    res = s21_atan(x / s21_sqrt(1 - x * x));
  return res;
}

// вычисляет арксинус

long double s21_atan(double x) {
  long double cur = 0.;
  int negative = 1;
  if (!s21_isnan(x)) {
    if (x < 0) {
      x = -x;
      negative = -1;
    }
    if (x <= 11. / 16) {
      long double k = 0., den = 0., sign = 1., y = (long double)x, prev = 1.;
      while (s21_fabs(cur - prev) > EPS) {
        prev = cur;
        den = 2 * k + 1;
        cur += sign * y / den;
        y *= x * x;
        k++;
        sign = -sign;
      }
    } else if (x < 19. / 16) {
      cur = PI / 4 + s21_atan((x - 1) / (x + 1));
    } else if (x < 39. / 16) {
      cur = 0.98279372324733 +
            s21_atan((x - 1.5) / (1 + 1.5 * x));  // atan(1.5) + ...
    } else {
      cur = PI / 2 + s21_atan(-1 / x);
    }
  } else {
    cur = S21_NAN;
  }
  return negative * cur;
}
// вычисляет арктангенс в пределах (-1;1), нужно расширить

long double s21_ceil(double x) {
  long long int res = s21_floor(x) + (x != (long long int)x);
  return (long double)((x > 1e17 || x < -1e17 || s21_isnan(x) || x == -0.)
                           ? x
                           : res);
}
// возвращает ближайшее целое число, не меньшее заданного значения по модулю

long double s21_cos(double x) { return s21_sin(PI / 2 + x); }

long double s21_exp(double x) {
  long double prev = 0., cur = 1., y = (long double)x;
  long double k = 1., den = 1.;
  int flag = 0;
  if (x < 0) {
    x = -x;
    y = -y;
    flag = 1;
  }
  while (s21_fabs(cur - prev) > EPS) {
    prev = cur;
    cur += y / den;
    y *= x;
    den *= (k + 1);
    k++;
  }
  return (flag) ? 1 / cur : cur;
}
// возвращает значение e, возведенное в заданную степень

long double s21_fabs(double x) {
  long double res = x;
  if (x == 0)
    res = 0;
  else if (x < 0)
    res = -x;
  return res;
}
// вычисляет абсолютное значение числа с плавающей точкой

long double s21_floor(double x) {
  long long int res = (long long int)x;
  if (x < 0 && x != (double)res) res--;
  return (long double)((x > 1e17 || x < -1e17 || s21_isnan(x) || x == -0.)
                           ? x
                           : res);
}
// возвращает ближайшее целое число, не превышающее заданное значение

long double s21_fmod(double x, double y) {
  long double x_mod = s21_fabs(x);
  long double y_mod = s21_fabs(y);
  long double dif = x_mod;
  if (y_mod == 0 || s21_isinf(x) || s21_isnan(y_mod) || s21_isnan(x_mod)) {
    dif = S21_NAN;
  } else if (s21_isinf(y_mod))
    dif = x;
  else {
    dif = x_mod / y_mod - s21_floor(x_mod / y_mod);
    dif *= y_mod;
    if (x < 0) dif = -dif;
  }
  return dif;
}

// остаток операции деления с плавающей точкой
long double s21_log(double x) {  // разложение для ln(1+z)
  long double ret;
  if (x == 0) {
    ret = -S21_INF;
  } else if (x < 0 || s21_isnan(x)) {  // проверка на NaN и отрицательное число
    ret = S21_NAN;
  } else if (s21_isinf(x)) {
    ret = S21_INF;
  } else {
    long double prev = 1., cur = 0., exp = s21_exp(1);
    int k = 0;
    int n = 0;
    if (x < 0.7) {
      while (x < 0.7) {  // ln(a/e^a) = ln(a)-b
        x *= exp;
        k++;
      }
    } else {
      while (x > 1.9) {  // полином аппроксимирует от 0 до 2 поэтому применяем
                         // ln(e^a * b) = a+ln(b)
        x /= exp;
        n++;
      }
    }

    long double z = x - 1;  // ln(x) = ln(z+1) = Σ(...) x = z+1 => z = x-1
    for (int i = 1; s21_fabs(cur - prev) > EPS; i++) {
      prev = cur;
      cur += ((i % 2 == 0) ? -1 : 1) * z / i;
      z *= x - 1;
    }
    ret = n + cur - k;
  }
  return ret;
}
// вычисляет натуральный логарифм

long double s21_pow(double base, double exp) {  // x^y = e^(y*ln(x))
  long double res;
  if (exp == 0)
    res = 1;
  else if (s21_isnan(base) || s21_isnan(exp))
    res = S21_NAN;
  else if (base == 0) {
    if (exp > 0)
      res = base;
    else
      res = S21_INF;
  } else if (s21_isinf(base) && s21_isinf(exp))
    res = (exp < 0) ? 0 : S21_INF;
  else {
    if (base > 0)
      res = s21_exp(exp * s21_log(base));
    else if (s21_ceil(exp) == exp)
      if (((long long int)exp) % 2 == 1)
        res = -s21_exp(exp * s21_log(-base));
      else
        res = s21_exp(exp * s21_log(-base));
    else
      res = S21_NAN;
  }

  return res;
}

long double s21_sin(double x) {
  long double res;
  int res_sign = (x > 0) ? 1 : -1;  // sin(-x) = -sin(x)
  x = res_sign * x;
  if (s21_isinf(x) || s21_isnan(x)) {
    res = S21_NAN;
  } else {
    // Сокращение на период
    long long int count_PI = (long long int)s21_floor(x / PI);
    x -= count_PI * PI;
    if (count_PI % 2 == 1) res_sign *= -1;
    long double prev = 1., cur = 0., y = (long double)x;
    long double k = 1., den = 1., sign = 1.;
    while (s21_fabs(cur - prev) > EPS) {
      prev = cur;
      cur += sign * y / den;
      y *= x * x;
      den *= (k + 1) * (k + 2);
      k += 2;
      sign = -sign;
    }
    res = cur;
  }

  return res_sign * res;
}

long double s21_sqrt(double x) { return s21_pow(x, 0.5); }
// вычисляет квадратный корень

long double s21_tan(double x) {
  long double res = s21_sin(x) / s21_cos(x);
  return res;
}
// вычисляет тангенс

int s21_isnan(long double x) { return (x != x) ? 1 : 0; }
// проверяет число на nan

int s21_isinf(long double x) {
  return ((x == S21_INF || x == -S21_INF) && !s21_isnan(x)) ? 1 : 0;
}
// проверяет число на inf

long double s21_tolerance(long double x) {
  while (x >= 1e16 && x < S21_INF) x /= 10;
  return x;
}