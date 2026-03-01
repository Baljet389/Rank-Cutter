#include "testSVD.h"



void testSVD(uint32_t r, uint32_t c, uint32_t precision, bool print) {
    Matrix<double> mat(r, c);
    fillMatrixRandomValues(mat, double(0), double(1));
    SVD<double> svd = calcSVD(mat);

    assertOrthorgonality(svd.U);
    assertOrthorgonality(svd.V);
    assertSVD<double>(svd, mat);

    if (print) {
        std::cout << "---Matrix---\n";
        printMatrixMatlab(mat, precision);
        std::cout << "---U---\n";
        printMatrixMatlab(svd.U, precision);
        std::cout << "---Sigma---\n";
        printMatrixMatlab(svd.S, precision);
        std::cout << "---V---\n";        
        printMatrixMatlab(svd.V, precision);
    }
}
void testQR(uint32_t r, uint32_t c, uint32_t precision, bool print) {
    Matrix<double> mat(r, c);
    fillMatrixRandomValues(mat, double(0), double(1));
    QR<double> qr = calcQRBlocked(mat, true);

    assertOrthorgonality(qr.Q);
    assertQR(qr, mat);

    if (print) {
        std::cout << "---Matrix---\n";
        printMatrixMatlab(mat, precision);
        std::cout << "---Q---\n";
        printMatrixMatlab(qr.Q, precision);
        std::cout << "---R---\n";
        printMatrixMatlab(qr.R, precision);
    }
}
