# kinetics_DNA_Hybridization
Gillespie script to simulate the effects of binding during DNA hybridization

__Compile the code__:

`g++ -std=c++17 -O3 fileName.cpp -o kDNA`

__Run the code__:

seq = 'ACATTTAGAGTAGTCCTTGGAGATTTTATGGAGATG'
stop = 1000

`./kDNA --seq $seq --stop $stop`
