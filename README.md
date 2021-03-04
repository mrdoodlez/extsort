# extsort
External sort algorithm demo

gendata - data pattern generator
================================
Generates binary file fiiled with 4-byte unsigned integers in ascending, descending or uniform random order
Please run

./gendata -h

for more info

extsort - extrnal sorting utility
=================================
Performs multithreaded external sorting in ascending order of file specified.
Requires additional sizeof(input file) disk space.
Input file will be overwritten.

Please run

./extsort -h

for more info.

Algorithms used: external sorting & merge, radix counting sort.
Radix counting sort is proved to have const * N performance.
