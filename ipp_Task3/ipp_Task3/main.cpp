#include <stdio.h>
#include <ctime>
#include <stdio.h>
#include <windows.h>
#include <mmsystem.h>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>

#include <iostream>

// количество строк в исходной квадратной матрице
const int MATRIX_SIZE = 1500;

/// Функция IntiMatrix() заполняет переданную в качестве 
/// параметра квадратную матрицу случайными значениями
/// matrix - исходная матрица СЛАУ
void IntiMatrix(double** matrix)
{
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		for (int j = 0; j <= MATRIX_SIZE; ++j)
		{
			matrix[i][j] = rand() % 250 + 1;
		}
	}
}

/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная квадратная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной квадратной матрице
/// result - массив ответов СЛАУ
int SerialGaussMethod(double **matrix, const int rows, double* result)
{
	int k;
	double koef;

	DWORD starttime, elapsedtime;
	starttime = timeGetTime();

	// прямой ход метода Гаусса
	for (k = 0; k < rows; ++k)
	{
		//
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];
			for (int j = k; j <= rows; ++j)
				matrix[i][j] += koef * matrix[k][j];
		}
	}

	elapsedtime = timeGetTime() - starttime;
	std::cout << "Время работы прямого хода: " << elapsedtime << " миллисекундах.\n";

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];
	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		//
		for (int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}
		result[k] /= matrix[k][k];
	}
	return elapsedtime;
}

/// Функция SerialGaussParallelMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная квадратная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной квадратной матрице
/// result - массив ответов СЛАУ
int SerialGaussParallelMethod(double **matrix, const int rows, double* result)
{
	DWORD starttime, elapsedtime;
	starttime = timeGetTime();

	// прямой ход метода Гаусса
	for (int k = 0; k < rows; ++k)
	{
		//
		cilk_for(int i = k + 1; i < rows; ++i)
		{
			double koef = -matrix[i][k] / matrix[k][k];
			for (int j = k; j <= rows; ++j)
				matrix[i][j] += koef * matrix[k][j];
		}
	}

	elapsedtime = timeGetTime() - starttime;
	std::cout << "Время работы прямого хода: " << elapsedtime << " миллисекундах.\n";

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];
	for (int k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		//
		cilk::reducer_opadd<double> summ(0.0);
		cilk_for (int j = k + 1; j < rows; ++j)
		{
			summ += matrix[k][j] * result[j];
		}
		result[k] -= summ.get_value();
		result[k] /= matrix[k][k];
	}
	return elapsedtime;
}

int main()
{
	int consistent = 0;
	int parallel = 0;

	setlocale(LC_ALL, "Russian");
	srand((unsigned)time(0));

	int i;
	// кол-во строк в матрице, приводимой в качестве примера
	/*const int test_matrix_lines = 4;

	double **test_matrix = new double*[test_matrix_lines];

	// (test_matrix_lines+1) - количество столбцов в тестовой матрице,
	// последний столбец матрицы отведен под правые части уравнений, входящих в СЛАУ
	for (i = 0; i < (test_matrix_lines + 1); ++i)
	{
		test_matrix[i] = new double[(test_matrix_lines + 1)];
	}

	// массив решений СЛАУ
	double *result = new double[test_matrix_lines];

	// инициализация тестовой матрицы
	test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;

	//consistent = SerialGaussMethod(test_matrix, test_matrix_lines, result);
	parallel = SerialGaussParallelMethod(test_matrix, test_matrix_lines, result);
	std::cout << "Ускорение: " << (double)consistent / parallel << ".\n";

	for (i = 0; i < test_matrix_lines; ++i)
	{
		delete[]test_matrix[i];
	}

	printf("Solution:\n");

	for (i = 0; i < test_matrix_lines; ++i)
	{
		printf("x(%d) = %lf\n", i, result[i]);
	}

	delete[] result;*/

	double **matrix = new double*[MATRIX_SIZE];
	for (i = 0; i < (MATRIX_SIZE + 1); ++i)
	{
		matrix[i] = new double[(MATRIX_SIZE + 1)];
	}
	
	double *result = new double[MATRIX_SIZE];

	IntiMatrix(matrix);
	
	consistent = SerialGaussMethod(matrix, MATRIX_SIZE, result);
	parallel = SerialGaussParallelMethod(matrix, MATRIX_SIZE, result);

	printf("\nSolution:\n");
	for (i = 0; i < MATRIX_SIZE; ++i)
		printf("x(%d) = %lf\n", i, result[i]);

	for (i = 0; i < MATRIX_SIZE; ++i)
		delete[]matrix[i];

	delete[] result;

	if ((consistent > 0) && (parallel > 0))
	std::cout << "Ускорение: " << (double)consistent / parallel << ".\n";
	system("pause");

	return 0;
}