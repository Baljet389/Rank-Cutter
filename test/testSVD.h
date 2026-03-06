#ifndef TEST_SVD_H
#define TEST_SVD_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "utility.h"
#include "svd.h"
#include "qr.h"


#define ASSERT(condition, msg) \
    do \
    { \
        if (!(condition)) \
        { \
            std::clog << "Fataler Fehler: " << msg << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (0)

namespace test {
    void testSVD(uint32_t r, uint32_t c, uint32_t precision, bool print, bool time = false);
    void testQR(uint32_t r, uint32_t c, uint32_t precision, bool print, bool time = false);
    template<typename A>
    void assertSVD(SVD<A>& svd, Matrix<A>& origiMat);
    template<typename A>
    void assertQR(QR<A>& qr, Matrix<A>& origiMat);
    template<typename A>
    void assertOrthorgonality(const Matrix<A>& mat);
    template<typename T, typename Derived>
    void printMatrixMatlab(MatrixInterface<T, Derived>& mat, uint32_t precision = 5);
}

template<typename A>
void test::assertSVD(SVD<A>& svd, Matrix<A>& origiMat) {
    Matrix<A>&  mat = svd.S;
    std::string whichMatrix =
      (mat.rows == mat.cols) ? "square" : ((mat.rows > mat.cols) ? "thin" : "wide");

    Matrix<A> finalRecon(mat.rows, mat.cols);
    reconstructSVD(svd, finalRecon, std::min(mat.rows, mat.cols));

    A maxDiff = A(0);
    for (uint32_t i = 0; i < mat.rows; i++)
    {
        for (uint32_t j = 0; j < mat.cols; j++)
        {
            A diff = std::abs(finalRecon.get(i, j) - origiMat.get(i, j));
            ASSERT(diff < 2.2e-10, "SVD Test failed; diff: " << diff << "\n");
            maxDiff = std::max(maxDiff, diff);
        }
    }
    std::cout << "SVD: Test successful for " << whichMatrix << " matrix! \n";
    std::cout << "Max difference: " << maxDiff << "\n";
}
template<typename A>
void test::assertQR(QR<A>& qr, Matrix<A>& origiMat) {
    Matrix<A>&  Q           = qr.Q;
    Matrix<A>&  R           = qr.R;
    std::string whichMatrix = (R.rows == R.cols) ? "square" : ((R.rows > R.cols) ? "thin" : "wide");
    Matrix<A>   recon(R.rows, R.cols);
    matrixMultMatrix(Q, R, recon);
    A maxDiff = A(0);
    for (uint32_t i = 0; i < R.rows; i++)
    {
        for (uint32_t j = 0; j < R.cols; j++)
        {
            A diff = std::abs(recon.get(i, j) - origiMat.get(i, j));
            ASSERT(diff < 2.2e-10, "QR Test failed");
            maxDiff = std::max(maxDiff, diff);
        }
    }
    std::cout << "QR: Test successful for " << whichMatrix << " matrix! \n";
    std::cout << "Max difference: " << maxDiff << "\n";
}
template<typename A>
void test::assertOrthorgonality(const Matrix<A>& mat) {
    for (uint32_t i = 0; i < mat.cols; i++)
    {
        for (uint32_t j = 0; j < mat.cols; j++)
        {
            A dot = 0;
            for (uint32_t k = 0; k < mat.rows; k++)
            {
                dot += mat.get(k, i) * mat.get(k, j);
            }
            if (i == j)
                ASSERT(std::abs(dot - A(1)) < 2.2e-10, "Orthogonal Test failed");
            else
                ASSERT(std::abs(dot) < 2.2e-10, "Orthogonal Test failed");
        }
    }
    std::cout << "Orthogonality test successful! \n";
}

template<typename T, typename Derived>
void test::printMatrixMatlab(MatrixInterface<T, Derived>& mat, uint32_t precision) {
    uint32_t rows = mat.getRows();
    uint32_t cols = mat.getCols();

    if (rows == 0 || cols == 0)
        return;

    std::vector<std::vector<std::string>> formattedData(rows, std::vector<std::string>(cols));
    std::vector<size_t>                   colWidths(cols, 0);

    for (uint32_t j = 0; j < cols; j++)
    {
        for (uint32_t i = 0; i < rows; i++)
        {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(precision) << mat.get(i, j);
            formattedData[i][j] = ss.str();

            colWidths[j] = std::max(colWidths[j], formattedData[i][j].length());
        }
    }

    for (uint32_t i = 0; i < rows; i++)
    {
        for (uint32_t j = 0; j < cols; j++)
        {
            std::cout << std::setw(colWidths[j]) << formattedData[i][j];

            if (j < cols - 1)
                std::cout << ", ";
        }
        std::cout << ";\n";
    }
}


#endif