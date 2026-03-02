#ifndef UTILITY_H
#define UTILITY_H

#include <random>
#include <chrono>
#include <cassert>

#include "matrix.h"

struct Timer;
template<typename A>
void fillMatrixRandomValues(Matrix<A>& mat, A disStart, A disStop);
template<typename A>
void reconstructSVD(SVD<A>& svd, Matrix<A>& res, uint32_t rank);
template<typename A>
inline A calculateDot(const A* vec1, const A* vec2, uint32_t size);
template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void transposeMultMatrix(MatrixInterface<A, DerivedA>& matA,
                         MatrixInterface<A, DerivedB>& matB,
                         MatrixInterface<A, DerivedR>& result);

template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void matrixMultMatrix(MatrixInterface<A, DerivedA>& matA,
                      const MatrixInterface<A, DerivedB>& matB,
                      MatrixInterface<A, DerivedR>& result);

// calculates A = A - B * C^T
template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void matrixMinusMatrixMultTranspose(MatrixInterface<A, DerivedA>& matA,
                                    MatrixInterface<A, DerivedB>& matB,
                                    const MatrixInterface<A, DerivedR>& matC);
// calculates A = A - B * C
template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void matrixMinusMatrixMultMatrix(MatrixInterface<A, DerivedA>& matA,
                                 MatrixInterface<A, DerivedB>& matB,
                                 const MatrixInterface<A, DerivedR>& matC);

struct Timer {
    std::chrono::steady_clock::time_point start;

    void                      startTimer() { start = std::chrono::steady_clock::now(); }
    std::chrono::milliseconds stopTimer() const {
        auto end = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    }
};
template<typename A>
void fillMatrixRandomValues(Matrix<A>& mat, A disStart, A disStop) {
    std::random_device          rd;
    std::mt19937                generator(rd());
    std::normal_distribution<A> distribution(disStart, disStop);
    for (uint32_t i = 0; i < mat.rows; i++)
    {
        for (uint32_t j = 0; j < mat.cols; j++)
        {
            mat(i, j) = distribution(generator);
        }
    }
}
template<typename A>
void reconstructSVD(SVD<A>& svd, Matrix<A>& res, uint32_t rank) {
    res.setValue(A(0));

    Matrix<A>& U = svd.U;
    Matrix<A>& S = svd.S;
    Matrix<A>& V = svd.V;

    uint32_t rows = S.getRows();
    uint32_t cols = S.getCols();

    for (uint32_t j = 0; j < cols; j++)
    {
        A* resPtr = res.getColumnPointer(j);
        for (uint32_t k = 0; k < rank; k++)
        {
            A* UColPtr = U.getColumnPointer(k);
            A  scalar  = S(k, k) * V(j, k);
            for (uint32_t i = 0; i < rows; i++)
            {
                resPtr[i] += UColPtr[i] * scalar;
            }
        }
    }
}

template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void transposeMultMatrix(MatrixInterface<A, DerivedA>& matA,
                         MatrixInterface<A, DerivedB>& matB,
                         MatrixInterface<A, DerivedR>& result) {

    static_assert(DerivedA::isContiguous());
    static_assert(DerivedB::isContiguous());
    assert(matA.getRows() == matB.getRows());

    const uint32_t     colsA = matA.getCols();
    const uint32_t     colsB = matB.getCols();
    const uint32_t     rowsB = matB.getRows();
    constexpr uint32_t BS    = 64;

    result.setValue(0);
#pragma omp parallel for schedule(static)
    for (uint32_t jj = 0; jj < colsA; jj += BS)
        for (uint32_t kk = 0; kk < rowsB; kk += BS)
            for (uint32_t ii = 0; ii < colsB; ii += BS)
                for (uint32_t r = jj; r < std::min(jj + BS, colsA); r++)
                {
                    A* colAPtr = matA.getColumnPointer(r);
                    for (uint32_t k = ii; k < std::min(ii + BS, colsB); k++)
                    {
                        A* colBPtr = matB.getColumnPointer(k);
                        A  dot     = 0;
                        for (uint32_t c = kk; c < std::min(kk + BS, rowsB); c++)
                        {
                            dot += colAPtr[c] * colBPtr[c];
                        }
                        result(r, k) += dot;
                    }
                }
}
template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void matrixMultMatrix(MatrixInterface<A, DerivedA>& matA,
                      const MatrixInterface<A, DerivedB>& matB,
                      MatrixInterface<A, DerivedR>& result) {

    static_assert(DerivedA::isContiguous());
    static_assert(DerivedB::isContiguous());
    assert(matA.getCols() == matB.getRows());

    const uint32_t     rowsA = matA.getRows();
    const uint32_t     colsA = matA.getCols();
    const uint32_t     colsB = matB.getCols();
    constexpr uint32_t BS    = 64;

    result.setValue(0);
#pragma omp parallel for schedule(static)
    for (uint32_t kk = 0; kk < rowsA; kk += BS)
        for (uint32_t ii = 0; ii < colsA; ii += BS)
            for (uint32_t jj = 0; jj < colsB; jj += BS)
                for (uint32_t r = jj; r < std::min(jj + BS, colsB); r++)
                {
                    A* resColPtr = result.getColumnPointer(r);
                    for (uint32_t k = ii; k < std::min(ii + BS, colsA); k++)
                    {
                        A  val        = matB(k, r);
                        A* matAColPtr = matA.getColumnPointer(k);
                        for (uint32_t c = kk; c < std::min(kk + BS, rowsA); c++)
                        {
                            resColPtr[c] += matAColPtr[c] * val;
                        }
                    }
                }
}

template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void matrixMinusMatrixMultTranspose(MatrixInterface<A, DerivedA>& matA,
                                    MatrixInterface<A, DerivedB>& matB,
                                    const MatrixInterface<A, DerivedR>& matC) {

    static_assert(DerivedA::isContiguous());
    static_assert(DerivedB::isContiguous());
    assert(matA.getRows() == matB.getRows());
    assert(matA.getCols() == matC.getRows());
    assert(matB.getCols() == matC.getCols());

    const uint32_t     rowsA = matA.getRows();
    const uint32_t     colsA = matA.getCols();
    const uint32_t     colsB = matB.getCols();
    constexpr uint32_t BS    = 64;

#pragma omp parallel for schedule(static)
    for (uint32_t jj = 0; jj < rowsA; jj += BS)
        for (uint32_t ii = 0; ii < colsB; ii += BS)
            for (uint32_t kk = 0; kk < colsA; kk += BS)
                for (uint32_t r = kk; r < std::min(kk + BS, colsA); r++)
                {
                    A* aColPtr = matA.getColumnPointer(r);
                    for (uint32_t k = ii; k < std::min(ii + BS, colsB); k++)
                    {
                        A* bColPtr = matB.getColumnPointer(k);
                        A  val     = matC(r, k);
                        for (uint32_t c = jj; c < std::min(jj + BS, rowsA); c++)
                        {
                            aColPtr[c] -= bColPtr[c] * val;
                        }
                    }
                }
}
template<typename A, typename DerivedA, typename DerivedB, typename DerivedR>
void matrixMinusMatrixMultMatrix(MatrixInterface<A, DerivedA>& matA,
                                 MatrixInterface<A, DerivedB>& matB,
                                 const MatrixInterface<A, DerivedR>& matC) {

    static_assert(DerivedA::isContiguous());
    static_assert(DerivedB::isContiguous());

    assert(matA.getRows() == matB.getRows());
    assert(matB.getCols() == matC.getRows());
    assert(matA.getCols() == matC.getCols());

    const uint32_t     rowsA  = matA.getRows();
    const uint32_t     colsA  = matA.getCols();
    const uint32_t     innerK = matB.getCols();
    constexpr uint32_t BS     = 64;

#pragma omp parallel for schedule(static)
    for (uint32_t kk = 0; kk < rowsA; kk += BS)
        for (uint32_t ii = 0; ii < innerK; ii += BS)
            for (uint32_t jj = 0; jj < colsA; jj += BS)
                for (uint32_t j = jj; j < std::min(jj + BS, colsA); j++)
                {
                    A* aColPtr = matA.getColumnPointer(j);

                    for (uint32_t k = ii; k < std::min(ii + BS, innerK); k++)
                    {
                        const A* bColPtr = matB.getColumnPointer(k);
                        A        val     = matC(k, j);

                        for (uint32_t i = kk; i < std::min(kk + BS, rowsA); i++)
                        {
                            aColPtr[i] -= bColPtr[i] * val;
                        }
                    }
                }
}
template<typename A>
inline A calculateDot(const A* vec1, const A* vec2, uint32_t size) {
    assert(vec1);
    assert(vec2);

    A res = A(0);
    for (uint32_t i = 0; i < size; i++)
    {
        res += vec1[i] * vec2[i];
    }
    return res;
}
#endif
