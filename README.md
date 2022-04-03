# kinetics_DNA_Hybridization
Gillespie script to simulate the effects of binding during DNA hybridization

Compile the code:
`g++ -std=c++17 -O3 fileName.cpp -o kDNA`

Run the code:
seq = 'ACATTTAGAGTAGTCCTTGGAGATTTTATGGAGATG'
stop = 1000

`./kDNA --seq $seq --stop $stop`
