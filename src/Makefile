CC=gcc
CFLAGS=-Wall -Werror -Wextra -std=c11
GCOV=-fprofile-arcs -ftest-coverage
# CHFLAGS=-lcheck
CHFLAGS=-lcheck -lsubunit -lm # for ubuntu
TEST_FILES_C=tests/*_test.c

all: s21_math.a

test: 
	make cover_lib
	$(CC) $(CFLAGS) $(GCOV) -c $(TEST_FILES_C)
	$(CC) $(CFLAGS) $(GCOV) s21_abs_test.o -o s21_abs_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_fabs_test.o -o s21_fabs_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_sin_test.o -o s21_sin_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_cos_test.o -o s21_cos_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_tan_test.o -o s21_tan_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_exp_test.o -o s21_exp_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_fmod_test.o -o s21_fmod_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_floor_test.o -o s21_floor_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_ceil_test.o -o s21_ceil_test s21_math.a $(CHFLAGS)

	$(CC) $(CFLAGS) $(GCOV) s21_asin_test.o -o s21_asin_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_acos_test.o -o s21_acos_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_atan_test.o -o s21_atan_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_log_test.o -o s21_log_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_pow_test.o -o s21_pow_test s21_math.a $(CHFLAGS)
	$(CC) $(CFLAGS) $(GCOV) s21_sqrt_test.o -o s21_sqrt_test s21_math.a $(CHFLAGS)
	make run_tests
	make clean
	make gcov_report

s21_math.a:
	$(CC) $(CFLAGS) -c s21_math.c -o s21_math.o
	ar rc s21_math.a s21_math.o
	ranlib s21_math.a
	make clean

cover_lib:	
	$(CC) $(CFLAGS) $(GCOV) -c s21_math.c -o s21_math.o
	ar rc s21_math.a s21_math.o
	ranlib s21_math.a
	make clean

gcov_report: clean_folder
	gcov *.c
	mkdir gcov_report
	lcov --capture --directory . --output-file gcov_report/report.info
	genhtml gcov_report/report.info --output-directory gcov_report
	open gcov_report/index.html

run_tests:
# valgrind ./[ваш бинарник] [аргументы] --leak-check=full
	./s21_abs_test
	./s21_fabs_test
	./s21_sin_test
	./s21_cos_test
	./s21_tan_test
	./s21_exp_test
	./s21_fmod_test
	./s21_floor_test
	./s21_ceil_test
	./s21_asin_test
	./s21_acos_test
	./s21_atan_test
	./s21_log_test
	./s21_pow_test
	./s21_sqrt_test

clean:
	rm *.o

clean_report: 
	rm *.gcda *.gcno *_test s21_math.a *.gcov

clean_folder:
	rm -rf gcov_report

rebuild: clean all