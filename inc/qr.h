#ifndef QR_H
#define QR_H


#include <vector>
#include <cmath>
#include <iostream>

#include "blockUpdate.h"
#include "matrix.h"
#include "householder.h"

template<typename A>
QR<A> calcQRBlocked(const Matrix<A>& mat, bool updateQ);
template<typename A>
void generatePanel(HouseholderWS<A>& ws, uint32_t currBlockSize, uint32_t col);


template<typename A>
QR<A> calcQRBlocked(const Matrix<A>& mat, bool updateQ) {
    Matrix<A> QMatrix(mat.rows, mat.rows);
    QMatrix.setIdentity();
    QR<A> qr(QMatrix, mat);

    Matrix<A>& Q            = qr.Q;
    Matrix<A>& R            = qr.R;
    uint32_t   minDimension = std::min(R.rows, R.cols);
    uint32_t   blockSize    = 64;

    Matrix<A> V(R.rows, blockSize);
    Matrix<A> T(blockSize, blockSize);

    Matrix<A> W(blockSize, R.cols);
    Matrix<A> Y(blockSize, R.cols);

    Matrix<A> W2(R.rows, blockSize);
    Matrix<A> Y2(R.rows, blockSize);

    Matrix<A> VTMultVStorage(blockSize, blockSize);

    HouseholderWS<A> ws(R, V, T, VTMultVStorage);
    ws.reserve(R.rows);
    for (uint32_t i = 0; i < minDimension; i += blockSize)
    {
        uint32_t currBlockSize = std::min(blockSize, minDimension - i);
        V.setValue(A(0));
        T.setValue(A(0));
        generatePanel(ws, currBlockSize, i);


        const uint32_t nextCol       = i + currBlockSize;
        const uint32_t trailingCols  = R.cols - nextCol;
        const uint32_t currentHeight = R.rows - i;

        SubMatrixView<A> activeV(V, i, 0, currentHeight, currBlockSize);
        SubMatrixView<A> activeT(T, 0, 0, currBlockSize, currBlockSize);
        SubMatrixView<A> trailingR(R, i, nextCol, currentHeight, trailingCols);
        SubMatrixView<A> trailingQ(Q, 0, i, Q.rows, currentHeight);

        BlockUpdateWorkspace<A> blockws(W, Y, W2, Y2, activeV, activeT, trailingR, trailingQ);

        // Update R = R - V * T^T * V^T * R

        if (trailingCols > 0)
        {
            applyBlockReflectorToMiddle(blockws);
        }

        // Update Q = Q - Q * V * T * V^T

        if (updateQ)
        {
            applyBlockReflectorToLeft(blockws);
        }
    }
    return qr;
}

template<typename A>
void generatePanel(HouseholderWS<A>& ws, uint32_t currBlockSize, uint32_t col) {


    for (uint32_t i = col; i < col + currBlockSize; i++)
    {
        applyLeftHouseholder(ws, i, col, col + currBlockSize);
    }
    updateTMatrix(ws, col, currBlockSize);
}


#endif