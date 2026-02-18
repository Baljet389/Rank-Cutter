#ifndef BIDIAGONAL_H
#define BIDIAGONAL_H

#include "householder.h"
#include "blockUpdate.h"
#include "matrix.h"

template<typename A>
void calculateBidiagonalFormBlocked(SVD<A>& svd);
template<typename A>
void generateBidiagPanel(HouseholderWS<A>& wsLeft,
                         HouseholderWS<A>& wsRight,
                         uint32_t          currBlockSize,
                         uint32_t          col);
template<typename A>
A calculateBidiagonalFrobreniusNorm(const SubMatrixView<A>& subMat);

template<typename A>
void calculateBidiagonalFormBlocked(SVD<A>& svd) {
    Matrix<A>& U       = svd.U;
    Matrix<A>& S       = svd.S;
    Matrix<A>& V       = svd.V;
    bool       updateU = (U.data.size() > 0);
    bool       updateV = (V.data.size() > 0);

    const uint32_t     rows         = S.rows;
    const uint32_t     cols         = S.cols;
    const uint32_t     minDimension = std::min(rows, cols);
    constexpr uint32_t blockSize    = 3;

    Matrix<A> VL(rows, blockSize);
    Matrix<A> VR(cols, blockSize);

    Matrix<A> TL(blockSize, blockSize);
    Matrix<A> TR(blockSize, blockSize);
    Matrix<A> VTMultVStorage(blockSize, blockSize);

    Matrix<A>        placeholder(0, 0);
    SubMatrixView<A> subPlaceholder(placeholder, 0, 0, 0, 0);

    Matrix<A> W(rows, blockSize);
    Matrix<A> Y(rows, blockSize);

    HouseholderWS<A> wsLeft(S, VL, TL, VTMultVStorage);
    HouseholderWS<A> wsRight(S, VR, TR, VTMultVStorage);
    wsLeft.reserve(rows);
    wsRight.reserve(cols);

    for (uint32_t i = 0; i < minDimension; i += blockSize)
    {
        const uint32_t currBlockSize = std::min(blockSize, minDimension - i);

        VL.setValue(A(0));
        VR.setValue(A(0));

        TL.setValue(A(0));
        TR.setValue(A(0));

        generateBidiagPanel(wsLeft, wsRight, currBlockSize, i);

        const uint32_t currentHeight = rows - i;

        if (updateU)
        {
            SubMatrixView<A> activeV(VL, i, 0, currentHeight, currBlockSize);
            SubMatrixView<A> activeT(TL, 0, 0, currBlockSize, currBlockSize);
            SubMatrixView<A> trailingU(U, 0, i, rows, currentHeight);

            BlockUpdateWorkspace<A> blockwsLeft(placeholder, placeholder, W, Y, activeV, activeT,
                                                subPlaceholder, trailingU);
            applyBlockReflectorToLeft(blockwsLeft);
        }

        if (!updateV || (i + 1 >= cols))
            continue;

        const uint32_t currentWidth = cols - (i + 1);

        SubMatrixView<A> activeVR(VR, i + 1, 0, currentWidth, currBlockSize);
        SubMatrixView<A> activeTR(TR, 0, 0, currBlockSize, currBlockSize);
        SubMatrixView<A> trailingV(V, 0, i + 1, cols, currentWidth);

        BlockUpdateWorkspace<A> blockwsRight(placeholder, placeholder, W, Y, activeVR, activeTR,
                                             subPlaceholder, trailingV);
        applyBlockReflectorToRight(blockwsRight);
    }
}
template<typename A>
void generateBidiagPanel(HouseholderWS<A>& wsLeft,
                         HouseholderWS<A>& wsRight,
                         uint32_t          currBlockSize,
                         uint32_t          col) {
    for (uint32_t i = col; i < col + currBlockSize; i++)
    {
        applyLeftHouseholder(wsLeft, i, col, wsLeft.R.cols);

        if (i != wsRight.R.cols - 1)
            applyRightHouseholder(wsRight, i, col);
    }
    updateTMatrix(wsLeft, col, currBlockSize);
    updateTMatrix(wsRight, col, currBlockSize);
}
template<typename A>
A calculateBidiagonalFrobreniusNorm(const SubMatrixView<A>& subMat) {
    A        sum          = 0;
    uint32_t minDimension = std::min(subMat.rows, subMat.cols);
    for (uint32_t i = 0; i < minDimension; i++)
    {
        A diag = subMat.get(i, i);
        sum += diag * diag;
        if (i < minDimension - 1)
        {
            A superDiag = subMat.get(i, i + 1);
            sum += superDiag * superDiag;
        }
    }
    return std::sqrt(sum);
}

#endif