#ifndef BLOCK_UPDATE_H
#define BLOCK_UPDATE_H

#include "matrix.h"
#include "utility.h"

template<typename A>
struct BlockUpdateWorkspace;
template<typename A>
void applyBlockReflectorToMiddle(BlockUpdateWorkspace<A>& blockws);
template<typename A>
void applyBlockReflectorToLeft(BlockUpdateWorkspace<A>& blockws);
template<typename A>
void applyBlockReflectorToRight(BlockUpdateWorkspace<A>& ws);


template<typename A>
struct BlockUpdateWorkspace {
    Matrix<A>& W;
    Matrix<A>& Y;

    Matrix<A>& W2;
    Matrix<A>& Y2;

    SubMatrixView<A>& activeV;
    SubMatrixView<A>& activeT;
    SubMatrixView<A>& trailingMiddle;
    SubMatrixView<A>& trailingSide;
    BlockUpdateWorkspace(Matrix<A>&        w,
                         Matrix<A>&        y,
                         Matrix<A>&        w2,
                         Matrix<A>&        y2,
                         SubMatrixView<A>& v,
                         SubMatrixView<A>& t,
                         SubMatrixView<A>& m,
                         SubMatrixView<A>& s) :
        W(w),
        Y(y),
        W2(w2),
        Y2(y2),
        activeV(v),
        activeT(t),
        trailingMiddle(m),
        trailingSide(s) {}
};
template<typename A>
void applyBlockReflectorToMiddle(BlockUpdateWorkspace<A>& blockws) {
    SubMatrixView<A>& trailingMiddle = blockws.trailingMiddle;
    SubMatrixView<A>& activeV        = blockws.activeV;
    SubMatrixView<A>& activeT        = blockws.activeT;

    const uint32_t currBlockSize = activeT.rows;
    // const uint32_t currentHeight = activeV.rows;
    const uint32_t trailingCols = trailingMiddle.cols;

    SubMatrixView<A> W(blockws.W, 0, 0, currBlockSize, trailingCols);
    SubMatrixView<A> Y(blockws.Y, 0, 0, currBlockSize, trailingCols);

    // Step A: W = V^T * R_trailing
    // Size: blockSize x trailingCols = (currBlockSize x currentHeight) * (currentHeight x trailingCols)

    transposeMultMatrix(activeV, trailingMiddle, W);

    // Step B: Y = T^T * W
    // Size: blockSize x trailingCols = (currBlockSize x currBlockSize) * (currBlockSize x trailingCols)

    transposeMultMatrix(activeT, W, Y);

    // Step C: R_trailing = R_trailing - V * Y
    // Size: (currentHeight x trailingCols) = (currentHeight x currBlockSize) * (currBlockSize x trailingCols)
    matrixMinusMatrixMultMatrix(trailingMiddle, activeV, Y);
}
template<typename A>
void applyBlockReflectorToLeft(BlockUpdateWorkspace<A>& blockws) {
    SubMatrixView<A>& trailingLeft = blockws.trailingSide;
    SubMatrixView<A>& activeV      = blockws.activeV;
    SubMatrixView<A>& activeT      = blockws.activeT;

    const uint32_t currBlockSize = activeT.rows;
    const uint32_t rows          = trailingLeft.rows;
    // const uint32_t cols          = trailingLeft.cols;

    SubMatrixView<A> W2(blockws.W2, 0, 0, rows, currBlockSize);
    SubMatrixView<A> Y2(blockws.Y2, 0, 0, rows, currBlockSize);

    // Step A: W = trailingLeft * V
    // Size: (rows, cols) * (cols, currBlockSize) = (rows, currBlockSize)

    matrixMultMatrix(trailingLeft, activeV, W2);

    // Step B: Y = W * T
    // Size: (rows, currBlockSize) * (currBlockSize, currBlockSize) = (rows, currBlockSize)
    matrixMultMatrix(W2, activeT, Y2);

    // Step C: trailingLeft = trailingLeft - W_updated * V^T
    // Size: (rows, currBlockSize) * (currBlockSize, cols) = (rows, cols)

    matrixMinusMatrixMultTranspose(trailingLeft, Y2, activeV);
}
template<typename A>
void applyBlockReflectorToRight(BlockUpdateWorkspace<A>& ws) {


    SubMatrixView<A>& trailingRight = ws.trailingSide;
    SubMatrixView<A>& activeVR      = ws.activeV;
    SubMatrixView<A>& activeTR      = ws.activeT;

    const uint32_t currBlockSize = activeTR.rows;
    const uint32_t rows          = trailingRight.rows;

    SubMatrixView<A> W2(ws.W2, 0, 0, rows, currBlockSize);
    SubMatrixView<A> Y2(ws.Y2, 0, 0, rows, currBlockSize);

    // Step A: W = trailingRight * V
    // Size: (rows, cols) * (cols, currBlockSize) = (rows, currBlockSize)
    matrixMultMatrix(trailingRight, activeVR, W2);

    // Step B: Y = W * T
    // Size: (rows, currBlockSize) * (currBlockSize, currBlockSize) = (rows, currBlockSize)
    matrixMultMatrix(W2, activeTR, Y2);

    // Step C: trailingRight = trailingRight - W_updated * V^T
    // Size: (rows, currBlockSize) * (currBlockSize, cols) = (rows, cols)
    matrixMinusMatrixMultTranspose(trailingRight, Y2, activeVR);
}

#endif
