#include "testSVD.h"


void test::testSVD(uint32_t r, uint32_t c, uint32_t precision, bool print, bool time) {
    Matrix<double> mat(r, c);
    fillMatrixRandomValues(mat, double(0), double(1));
    Timer t;
    t.startTimer();
    if (time)
        std::cout << "Start Timer\n";

    SVD<double> svd = calcSVD(mat);
    if (time)
        std::cout << "SVD took " << t.stopTimer().count() << "\n";


    assertOrthorgonality(svd.U);
    assertOrthorgonality(svd.V);
    assertSVD<double>(svd, mat);

    if (print)
    {
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
void test::testQR(uint32_t r, uint32_t c, uint32_t precision, bool print, bool time) {
    Matrix<double> mat(r, c);
    fillMatrixRandomValues(mat, double(0), double(1));
    Timer t;
    t.startTimer();
    if (time)
        std::cout << "Start Timer\n";
    QR<double> qr = calcQRBlocked(mat, true);
    if (time)
        std::cout << "QR took " << t.stopTimer().count() << "\n";

    assertOrthorgonality(qr.Q);
    assertQR(qr, mat);

    if (print)
    {
        std::cout << "---Matrix---\n";
        printMatrixMatlab(mat, precision);
        std::cout << "---Q---\n";
        printMatrixMatlab(qr.Q, precision);
        std::cout << "---R---\n";
        printMatrixMatlab(qr.R, precision);
    }
}
