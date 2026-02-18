#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include "matrix.h"

template<typename A>
struct HouseholderWS;
template<typename A>
A calcNormedHouseholder(std::vector<A>& v);
template<typename A>
void applyRightHouseholder(HouseholderWS<A>& ws, uint32_t i, uint32_t col);
template<typename A>
void applyLeftHouseholder(HouseholderWS<A>& ws, uint32_t i, uint32_t col, uint32_t endColumn);
template<typename A>
void updateTMatrix(HouseholderWS<A>& ws, uint32_t col, uint32_t currBlockSize);

template<typename A>
struct HouseholderWS {
    Matrix<A>&     R;
    Matrix<A>&     V;
    Matrix<A>&     T;
    std::vector<A> vL;
    Matrix<A>&     VTMultVStorage;
    HouseholderWS(Matrix<A>& r, Matrix<A>& v, Matrix<A>& t, Matrix<A>& vTMultVStorage) :
        R(r),
        V(v),
        T(t),
        VTMultVStorage(vTMultVStorage) {}

    void reserve(uint32_t n) { vL.reserve(n); }
};
template<typename A>
void applyRightHouseholder(HouseholderWS<A>& ws, uint32_t i, uint32_t col) {
    std::vector<A>& vR = ws.vL;
    Matrix<A>&      VR = ws.V;
    Matrix<A>&      R  = ws.R;

    const uint32_t startColumn = i + 1;
    const uint32_t n           = R.cols - startColumn;
    const uint32_t m           = R.rows - i;
    vR.resize(n);

    for (uint32_t j = 0; j < n; j++)
        vR[j] = R(i, i + 1 + j);

    A tauRight = calcNormedHouseholder(vR);

    for (uint32_t j = 0; j < n; j++)
        VR(j + i + 1, i - col) = vR[j];


    std::vector<A> dots(m, A(0));
    for (uint32_t k = 0; k < n; ++k)
    {
        const A vk  = vR[k];
        A*      col = R.getColumnPointer(startColumn + k) + i;

        for (uint32_t j = 0; j < m; ++j)
            dots[j] += col[j] * vk;
    }
    for (uint32_t k = 0; k < n; ++k)
    {
        const A scale = tauRight * vR[k];
        A*      col   = R.getColumnPointer(startColumn + k) + i;

        for (uint32_t j = 0; j < m; ++j)
            col[j] -= scale * dots[j];
    }

    ws.T(i - col, i - col) = tauRight;
}
template<typename A>
void applyLeftHouseholder(HouseholderWS<A>& ws, uint32_t i, uint32_t col, uint32_t endColumn) {

    std::vector<A>& vL = ws.vL;
    Matrix<A>&      R  = ws.R;
    Matrix<A>&      V  = ws.V;
    Matrix<A>&      T  = ws.T;

    const uint32_t m = R.rows - i;
    vL.resize(m);

    for (uint32_t j = 0; j < m; j++)
        vL[j] = R(i + j, i);

    A tauLeft = calcNormedHouseholder(vL);

    A* VColumnPtr = V.getColumnPointer(i - col) + i;
    for (uint32_t j = 0; j < m; j++)
        VColumnPtr[j] = vL[j];

    for (uint32_t j = i; j < endColumn; j++)
    {
        A* RColumnPtr = R.getColumnPointer(j) + i;
        A  dot        = calculateDot(vL.data(), RColumnPtr, m);
        A  prod       = tauLeft * dot;
        for (uint32_t k = 0; k < m; k++)
            RColumnPtr[k] -= vL[k] * prod;
    }

    T(i - col, i - col) = tauLeft;
}
template<typename A>
A calcNormedHouseholder(std::vector<A>& v) {
    A sigma = A(0);
    for (uint32_t i = 1; i < v.size(); ++i)
        sigma += v[i] * v[i];

    if (sigma == A(0))
    {
        return A(0);
    }

    A x0   = v[0];
    A norm = std::sqrt(x0 * x0 + sigma);

    A alpha = (x0 <= A(0)) ? norm : -norm;

    A tau = (alpha - x0) / alpha;

    A inv = A(1) / (x0 - alpha);
    v[0]  = A(1);
    for (uint32_t i = 1; i < v.size(); ++i)
        v[i] *= inv;

    return tau;
}
template<typename A>
void updateTMatrix(HouseholderWS<A>& ws, uint32_t col, uint32_t currBlockSize) {

    // $$T_j = \begin{pmatrix} T_{j-1} & -\tau_j T_{j-1} V_{j-1}^T v_j \\ 0 & \tau_j \end{pmatrix}$$

    Matrix<A>& V   = ws.V;
    Matrix<A>& T   = ws.T;
    Matrix<A>& res = ws.VTMultVStorage;

    const uint32_t effectiveRows = V.rows - col;

    SubMatrixView<A> Vblock(V, col, 0, effectiveRows, currBlockSize);
    transposeMultMatrix(Vblock, Vblock, res);

    for (uint32_t i = 0; i < currBlockSize; i++)
    {
        A* TColumnPtrI = T.getColumnPointer(i);
        for (uint32_t k = 0; k < i; k++)
        {
            A* TColumnPtrK = T.getColumnPointer(k);
            A  factor      = res(k, i);
            for (uint32_t j = 0; j <= k; j++)
            {
                TColumnPtrI[j] += TColumnPtrK[j] * factor;
            }
        }
        A leadingElement = -T(i, i);
        for (uint32_t j = 0; j < i; j++)
        {
            TColumnPtrI[j] *= leadingElement;
        }
    }
}

#endif