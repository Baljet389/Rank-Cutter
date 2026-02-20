# Rank-Cutter ✂️

**Rank-Cutter** (or `rac`) is a command-line utility for image compression using **Singular Value Decomposition (SVD)**. By decomposing an image into its constituent singular values and vectors, Rank-Cutter allows you to "cut" the rank of the image matrix, reducing storage size while maintaining structural integrity (although it will not actually decrease storage size).

---

## 🚀 Features

* **SVD-Based Compression:** Leverage linear algebra to compress images by retaining only the most significant singular values.
* **Randomized SVD:** Includes an optimized randomized algorithm for faster decomposition on high-resolution images.
* **Header-Only Dependencies:** Minimal footprint using `stb_image` and `argparse`.
* **Modern C++:** Written in C++17 for efficiency and clean syntax.

---

## 🛠️ Prerequisites

* **C++17** compatible compiler (GCC 7+, Clang 5+, or MSVC 2017+)
* **CMake** 3.10 or higher
* **Dependencies** (Included as headers):
    * `stb_image.h` / `stb_image_write.h` (Image I/O)
    * `argparse` (Argument parsing)

---

## 🔨 Building the Project

To build the executable `rac`, use the standard CMake workflow:

```bash
cmake -S . -B build
cmake --build build
```

The executable will be generated in `build`.

---

## 🖥️ Usage

Run the executable by providing the input path and the desired output path.

```bash
rac [options] input output
```

### Arguments

| Argument | Description |
| :--- | :--- |
| `input` | Path to the source image file (JPG, PNG, BMP, etc.). |
| `output` | Path where the compressed image will be saved. |

### Optional Flags

* `-h, --help` : Show the help message.
* `-v, --version` : Print version information.
* `-r, --rank` : Set the approximation rank (default: **10**). 
* `--rand` : Use the **Randomized SVD** algorithm for faster processing.

---

### Example

```bash
./rac input.jpg output.png --rank 50 --rand
```

## 📦 Dependencies

* [nothings/stb](https://github.com/nothings/stb)
* [p-ranav/argparse](https://github.com/p-ranav/argparse)
