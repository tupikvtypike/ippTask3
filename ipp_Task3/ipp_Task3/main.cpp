#include <stdio.h>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <chrono>

#include <iostream>

using namespace std::chrono;

// ���������� ����� � �������� ���������� �������
const int MATRIX_SIZE = 1500;

/// ������� IntiMatrix() ��������� ���������� � �������� 
/// ��������� ���������� ������� ���������� ����������
/// matrix - �������� ������� ����
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

/// ������� SerialGaussMethod() ������ ���� ������� ������ 
/// matrix - �������� ���������� ������� �������������� ���������, �������� � ����,
/// ��������� ������� ������� - �������� ������ ������ ���������
/// rows - ���������� ����� � �������� ���������� �������
/// result - ������ ������� ����
void SerialGaussMethod(double **matrix, const int rows, double* result)
{
	int k;
	double koef;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	// ������ ��� ������ ������
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

	end = std::chrono::system_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "����� ������ ������� ����: " << elapsed_seconds << " �����������.\n";

	// �������� ��� ������ ������
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
}

/// ������� SerialGaussParallelMethod() ������ ���� ������� ������ 
/// matrix - �������� ���������� ������� �������������� ���������, �������� � ����,
/// ��������� ������� ������� - �������� ������ ������ ���������
/// rows - ���������� ����� � �������� ���������� �������
/// result - ������ ������� ����
void SerialGaussParallelMethod(double **matrix, const int rows, double* result)
{
	int k;
	double koef;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	// ������ ��� ������ ������
	for (k = 0; k < rows; ++k)
	{
		//
		cilk_for(int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];
			for (int j = k; j <= rows; ++j)
				matrix[i][j] += koef * matrix[k][j];
		}
	}

	end = std::chrono::system_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "����� ������ ������� ����: " << elapsed_seconds << " �����������.\n";

	// �������� ��� ������ ������
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];
	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		//
		cilk_for (int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}
		result[k] /= matrix[k][k];
	}
}

int main()
{
	setlocale(LC_ALL, "Russian");
	srand((unsigned)time(0));

	int i;
	// ���-�� ����� � �������, ���������� � �������� �������
	//const int test_matrix_lines = 4;

	//double **test_matrix = new double*[test_matrix_lines];

	// (test_matrix_lines+1) - ���������� �������� � �������� �������,
	// ��������� ������� ������� ������� ��� ������ ����� ���������, �������� � ����
	/*for (i = 0; i < (test_matrix_lines + 1); ++i)
	{
		test_matrix[i] = new double[(test_matrix_lines + 1)];
	}*/

	// ������ ������� ����
	//double *result = new double[test_matrix_lines];

	// ������������� �������� �������
	//test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	//test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	//test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	//test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;

	double **matrix = new double*[MATRIX_SIZE];
	for (i = 0; i < (MATRIX_SIZE + 1); ++i)
	{
		matrix[i] = new double[(MATRIX_SIZE + 1)];
	}
	
	double *result = new double[MATRIX_SIZE];

	IntiMatrix(matrix);
	
	SerialGaussMethod(matrix, MATRIX_SIZE, result);

	printf("\nSolution:\n");
	for (i = 0; i < MATRIX_SIZE; ++i)
		printf("x(%d) = %lf\n", i, result[i]);

	for (i = 0; i < MATRIX_SIZE; ++i)
		delete[]matrix[i];

	delete[] result;

	system("pause");

	return 0;
}