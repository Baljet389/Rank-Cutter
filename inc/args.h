#ifndef ARGS_H
#define ARGS_H

#include <iostream>
#include <string>
#include <cstdint>
#include <cstdlib>

#include "argparse.hpp"
#include "image.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)


void parseArguments(int argc, char* argv[]);



#endif