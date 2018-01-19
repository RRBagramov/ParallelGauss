// ParallelGauss.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <vector>

class barrier {
	unsigned int const count;
	std::atomic<unsigned int> spaces;
	std::atomic<unsigned int> generation;
public:
	explicit barrier(unsigned nthreads) : count(nthreads), spaces(nthreads),
		generation(0) { }
	void wait() {
		unsigned const my_generation = generation;
		if (!--spaces) {
			spaces = count;
			++generation;
		}
		else {
			while (generation == my_generation)
				std::this_thread::yield();
		}
	}
};

void threadCalculateRow(float **matA, float *matB, float *x, int n, int start, int end, barrier& fowardStrokeBarrier, barrier& reverseMotionBarrier) {
	float coef = 0;
	for (int i = 0; i < n - 1; i++)
	{
		for (int k = i + 1; k <= end - 1; k++)
		{
			if ((k >= start) && (k < end)) {
				coef = matA[k][i] / matA[i][i];
				matB[k] -= coef * matB[i];
				for (int j = i; j <= n - 1; j++)
				{
					matA[k][j] = matA[k][j] - (coef * matA[i][j]);
				}
			}
		}
		fowardStrokeBarrier.wait();
	}

	reverseMotionBarrier.wait();

	float sumA = 0.0;
	for (int i = n - 2; i >= 0; i--)
	{
		if ((i >= start) && (i < end)) {
			for (int j = i + 1; j <= n - 1; j++)
			{
				sumA += matA[i][j] * x[j];
			}
			x[i] = (matB[i] - sumA) / matA[i][i];
			sumA = 0.0;
		}
		fowardStrokeBarrier.wait();
	}
}

void parallelGauss(float **matA, float *matB, float *x, int n, int MAX_THREAD) {
	std::vector<std::thread> threads;
	barrier fowardStrokeBarrier(MAX_THREAD), reverseMotionBarrier(MAX_THREAD);

	int temp = n / (MAX_THREAD - 1);

	for (int i = 0; i < (MAX_THREAD - 1); ++i) {
		int start = temp * i, end = temp * (i + 1);
		threads.emplace_back(threadCalculateRow, std::ref(matA), std::ref(matB), std::ref(x), n, start, end, std::ref(fowardStrokeBarrier), std::ref(reverseMotionBarrier));
	}

	int start, end;

	start = temp * (MAX_THREAD - 1);
	end = n;
	float coef = 0;
	for (int i = 0; i < n - 1; i++)
	{
		for (int k = i + 1; k <= end - 1; k++)
		{
			if ((k >= start) && (k < end)) {
				coef = matA[k][i] / matA[i][i];
				matB[k] -= coef * matB[i];
				for (int j = i; j <= n - 1; j++)
				{
					matA[k][j] = matA[k][j] - (coef * matA[i][j]);
				}
			}
		}
		fowardStrokeBarrier.wait();
	}

	x[n - 1] = matB[n - 1] / matA[n - 1][n - 1];

	reverseMotionBarrier.wait();

	start = temp * (MAX_THREAD - 1);
	end = n;
	float sumA = 0;
	for (int i = n - 2; i >= 0; i--)
	{
		if ((i >= start) && (i < end)) {
			for (int j = i + 1; j <= n - 1; j++)
			{
				sumA += matA[i][j] * x[j];
			}
			x[i] = (matB[i] - sumA) / matA[i][i];
			sumA = 0;
		}
		fowardStrokeBarrier.wait();
	}

	for (auto& thread : threads) {
		if (thread.joinable()){
			thread.join();
		}
	}
}
void matrixAFill(float **matA, int n)
{
	for (int i = 0; i < n; i++)
	{
		matA[i] = new float[n];
		for (int j = 0; j < n; j++)
		{
			matA[i][j] = rand() % 100 + 1;
		}
	}
}

void matrixBFill(float *matB, int n)
{
	for (int i = 0; i < n; i++)
	{
		matB[i] = rand() % 100 + 1;
	}
}

void showMatrixAB(float **matA, float *matB, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << matA[i][j] << "x" << j << " ";
		}
		std::cout << "= " << matB[i] << std::endl;
	}
	std::cout << std::endl;
}
void showX(float *x, int n)
{
	for (int i = 0; i < n; i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;
}

void standardGauss(float **matA, float *matB, float *x, int n)
{
	float sumA = 0;
	float coef = 0;
	for (int i = 0; i < n - 1; i++)
	{
		for (int k = i + 1; k <= n - 1; k++)
		{
			coef = matA[k][i] / matA[i][i];
			matB[k] -= coef * matB[i];
			for (int j = i; j <= n - 1; j++)
			{
				matA[k][j] = matA[k][j] - (coef * matA[i][j]);
			}
		}
	}

	x[n - 1] = matB[n - 1] / matA[n - 1][n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		for (int j = i + 1; j <= n - 1; j++)
		{
			sumA += matA[i][j] * x[j];
		}
		x[i] = (matB[i] - sumA) / matA[i][i];
		sumA = 0;
	}
}

void copyMatrix(float **matA, float **matAP, float *matB, float *matBP, int n)
{
	for (int i = 0; i < n; i++)
	{
		matAP[i] = new float[n];
		for (int j = 0; j < n; j++)
		{
			matAP[i][j] = matA[i][j];
		}
		matBP[i] = matB[i];
	}
}
int _tmain(int argc, _TCHAR* argv[])
{
	int n;
	std::cout << "Enter matrix size:" << std::endl;
	std::cin >> n;

	float **matA = new float*[n];
	float *matB = new float[n];
	float *x = new float[n];

	float **matAP = new float*[n];
	float *matBP = new float[n];
	float *xP = new float[n];

	//--------FILL-----------
	matrixAFill(matA, n);
	matrixBFill(matB, n);

	copyMatrix(matA, matAP, matB, matBP, n);

	//showMatrixAB(matA, matB, n);
	//--------FILL END-----------

	//--------STANDARD GAUSS-----------
	std::cout << "STANDARD GAUSS" << std::endl;
	clock_t t;

	t = clock();

	standardGauss(matA, matB, x, n);

	t = clock() - t;

	std::cout << "TIME: "<< t <<std::endl;
	//showMatrixAB(matA, matB, n);
	//showX(x, n);
	//--------STANDARD GAUSS END-----------

	//--------GAUSS PARALLEL--------

	std::cout << "GAUSS PARALLEL" << std::endl;
	int threadCount = 10;

	t = clock();

	parallelGauss(matAP, matBP, xP, n, threadCount);

	t = clock() - t;

	std::cout << "TIME: " << t << std::endl;
	//showMatrixAB(matAP, matBP, n);
	//showX(xP, n);
	//--------GAUSS PARALLEL END-----------

	return 0;
}

