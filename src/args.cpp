#include "args.h"

void parseArguments(int argc, char* argv[]) {
    argparse::ArgumentParser program(TOSTRING(PROJECT_NAME), TOSTRING(PROJECT_VERSION));

    program.add_argument("input").help("path to the input image file");

    program.add_argument("output").help("path to save the compressed image");

    program.add_argument("-r", "--rank")
      .help("approximation rank for compression")
      .default_value(static_cast<uint32_t>(10))
      .scan<'u', uint32_t>();

    program.add_argument("--rand")
      .help("use randomized SVD algorithm")
      .default_value(false)
      .implicit_value(true);

    try
    { program.parse_args(argc, argv); } catch (const std::runtime_error& err)
    {
        std::cerr << "Error: " << err.what() << std::endl;
        std::cerr << program;
        std::exit(EXIT_FAILURE);
    }

    std::string in         = program.get<std::string>("input");
    std::string out        = program.get<std::string>("output");
    uint32_t    rank       = program.get<uint32_t>("--rank");
    bool        randomized = program.get<bool>("--rand");

    compressImage<double>(in, out, rank, randomized);
}