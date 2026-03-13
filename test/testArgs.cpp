#include "testArgs.h"


void test::parseArgs(int argc, char* argv[]) {
    argparse::ArgumentParser program("rac_tests");

    program.add_argument("-r", "--rows")
      .help("Number of rows")
      .default_value(uint32_t(5))
      .scan<'u', uint32_t>();

    program.add_argument("-c", "--cols")
      .help("Number of columns")
      .default_value(uint32_t(5))
      .scan<'u', uint32_t>();

    program.add_argument("-p", "--precision")
      .help("Digits after decimal point for printing")
      .default_value(uint32_t(4))
      .scan<'u', uint32_t>();

    program.add_argument("--no-print")
      .help("Disable printing of matrices")
      .default_value(false)
      .implicit_value(true);

    program.add_argument("-t", "--time")
      .help("Stopping the time (using std::chrono)")
      .default_value(false)
      .implicit_value(true);

    argparse::ArgumentParser svd_command("svd");
    svd_command.add_description("Run SVD test");

    argparse::ArgumentParser qr_command("qr");
    qr_command.add_description("Run QR test");

    program.add_subparser(svd_command);
    program.add_subparser(qr_command);

    try
    { program.parse_args(argc, argv); } catch (const std::runtime_error& err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    auto rows  = program.get<uint32_t>("--rows");
    auto cols  = program.get<uint32_t>("--cols");
    auto prec  = program.get<uint32_t>("--precision");
    bool print = !program.get<bool>("--no-print");
    bool time  = program.get<bool>("--time");

    Timer t;
    t.startTimer();

    if (time)
        std::cout << "Start Timer: \n";

    if (program.is_subcommand_used(svd_command))
    {
        std::cout << "Running SVD: " << rows << "x" << cols << ", Precision: " << prec
                  << ", Print: " << (print ? "Yes" : "No") << "\n";
        testSVD(rows, cols, prec, print, time);
        if (time)
        {
            auto duration = t.stopTimer();
            std::cout << "SVD test took: " << duration.count() << "ms\n";
        }
    }
    else if (program.is_subcommand_used(qr_command))
    {
        std::cout << "Running QR: " << rows << "x" << cols << ", Precision: " << prec
                  << ", Print: " << (print ? "Yes" : "No") << "\n";
        testQR(rows, cols, prec, print, time);
        if (time)
        {
            auto duration = t.stopTimer();
            std::cout << "QR test took: " << duration.count() << "ms\n";
        }
    }
    else
    {
        std::cout << program;
        std::exit(1);
    }
}