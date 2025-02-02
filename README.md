# CppLEX: a Primal Simplex Algorithm implementation.

This is a C++ implementation of the Primal Simplex Algorithm, done within the scope of the course Programació Matemàtica from the Bachelor's Degree in Mathematics at UPC-BarcelonaTech.

## Installation
To install this repository, you can run the following command
```sh
git clone https://www.github.com/Atellas23/CppLEX
```
Once you have cloned the repository, you can `cd CppLEX` into the directory, and execute `make` to compile the code. After this steps, you are all set to use the program. If you want to, you can test the functions used in the code by issuing `make test` (NOTE: the results are only in spanish).

## Usage
Once you have the repository installed and the code compiled, you can use the program to find solutions to a Linear Programming problem in standard form, which is
```
(PL)    min z = c'x
        s.t.:
                Ax = b
                x >= 0
```
The problem data shall be submitted to the program in the following format: a plaintext file (the extension is not important, as long as it's a plaintext file) with the matrix A, followed by the line `end of matrix`, followed by a line with the restrictions b, followed by the line `end of restrictions`, followed by a line with the cost vector c, followed by the line `end of costs`. For example, the following (which is the content of the `data_model` file) would be an acceptable format for the input of our implementation:
```
1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1
end of matrix
0 0 0 0
end of restrictions
1 1 1 1
end of costs
```
Then, suppose this is stored in a file named `data.txt`. Then, you can simply execute `./simplex.exe data.txt` and you are good to go. When the algorithm finishes, you will have a report of the results both in the terminal and in a file named like the input data file with the extension `.out`.
