# otPSI-CA
We present a C++ implementation of our otPSI-CA protocol. Informally, over-threshold private set intersection cardinality (otPSI-CA) test is a two party protocol where the parties--a receiver and a sender--each hold a set of at most $n$ elements. The receiver learns only whether $\big|X\cap Y\big|\geq\tau$, while the sender learns nothing.

## Building
Our otPSI-CA protocol requires CMake (>= 3.16) and a C++ compiler (with support for C++17 standard). The protocol also requires NTL and Crypto++ which can be installed downloading the official site and using 
```
make && make install
```

To compile, run the following command from the project root:
```
mkdir build &&  cd build
cmake .. && cmake -- build .
```

## Benchmarking
We provide four benchmarking scripts in the root of project that can reproduce the benchmarking results from the paper.
```
~% otpsica_bench n tau eta
~% gs19_bench n tau eta
```
where $n$: the set size, $\tau$: a threshold, and $\eta$: the number of sets. Note that in GS19, $t$ is set to $n-\tau$.
Therefore, [GS19] targets on the setting where two sets are quite similar to each other.
