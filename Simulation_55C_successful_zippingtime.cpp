//GILLESPIE SIMULATION OF SUCCESSFUL ZIPPING TIME AT 55C

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <tuple>

double getEnergy(int j, int k, int jn, int kn);
std::pair<int, int> getParams(int yL, int xL, int len);
std::tuple<std::string, std::string, std::string, std::string> parseParams(int argc, char* argv[]);

int main(int argc, char* argv[]) {

    if (argc != 9) printf("Error: check input parameters!!!");

    auto [seq, num1, num2, stop] = parseParams(argc, argv);

    int randNum1 = std::stoi(num1);
    int randNum2 = std::stoi(num2);
    int stopCondition = std::stoi(stop);

    // make directory
    // std::system(("mkdir -p run-" + runIndex).c_str());

    int len = size(seq);

    std::vector<int> s1;
    std::vector<int> s2;
//In this simulation, unstructred sequences are coded with capital letters. For example: ACATTTAGAGTAGTCCTTGGAGATTTTATGGAGATG
//In stem-loop sequences, free tails are coded with capital letters, stem-loop regions are coded with lowercase letters. For example: AAGATGGTGAGTgccatcttAAAACTTACTGGAGAT
    for (auto& s : seq) {
        if (s == 'A') {
            s1.push_back(1);
            s2.push_back(2);
        } else if (s == 'T') {
            s1.push_back(2);
            s2.push_back(1);
        } else if (s == 'C') {
            s1.push_back(3);
            s2.push_back(4);
        } else if (s == 'G') {
            s1.push_back(4);
            s2.push_back(3);
        } else if (s == 'a') {
            s1.push_back(11);
            s2.push_back(22);
        } else if (s == 't') {
            s1.push_back(22);
            s2.push_back(11);
        } else if (s == 'c') {
            s1.push_back(33);
            s2.push_back(44);
        } else if (s == 'g') {
            s1.push_back(44);
            s2.push_back(33);
        }
    }

    std::vector<int> xVec;
    std::vector<int> yVec;

    double kForm = pow(10, 9);
    double time = 10000000.0;
    double t = 0.0;
    int xL = 0; int xR = 0;
    int yL = 0; int yR = 0;
    int g = 0;
    int hBonds = 0;
    int success = 0;
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::default_random_engine mt(seed);
    std::uniform_int_distribution<int> distInt(randNum1, randNum2);
    std::uniform_real_distribution<double> dist01(0, 1);

    while (t < time) {
        if (xL == 0 && xR == 0) {
            // int x, y;
            // make sure x == y
            int x = distInt(mt);
            int y = x;

            g = y - x;
            xL = x; xR = x;
            yL = y; yR = y;
            hBonds = 1;
            double r2 = dist01(mt);
            double kTotal = len*len*kForm;
            double tau = (-1.0/kTotal) * log(1.0 - r2);
            t += tau;
        }
        else {
            double randNum = dist01(mt);
            double r2 = dist01(mt);
            double kB1, kB2, kFL, kFR;
            if (xR == xL) {
                int k = s1[xL - 1];
                int kn = s1[xL];
                int j = s2[yL - 1];
                int jn = s2[yL];
                kB1 = kForm*exp(getEnergy(j, k, jn, kn));
                kB2 = 0.0;
            }
            else {
                int k = s1[xL - 1];
                int kn = s1[xL];
                int j = s2[yL - 1];
                int jn = s2[yL];
                kB1 = kForm*exp(getEnergy(j, k, jn, kn));
                int k_1 = s1[xR - 1];
                int kn_1 = s1[xR];
                int j_1 = s2[yR - 1];
                int jn_1 = s2[yR];
                kB2 = kForm*exp(getEnergy(j_1, k_1, jn_1, kn_1));
            }
            int rMax = getParams(yL, xL, len).first;
            int lMin = getParams(yL, xL, len).second;
            if (xL == lMin && xR != rMax) {
                kFL = 0.0;
                kFR = kForm;
            } else if (xR == rMax && xL != lMin) {
                kFR = 0.0;
                kFL = kForm;
            } else if (xR == rMax && xL == lMin) {
                kFL = kFR = 0.0;
            } else {
                kFL = kFR = kForm;
            }
            double kTotal = kB1 + kB2 + kFL + kFR;
            if (randNum >= 0 && randNum <= kB1 / kTotal){
                xL++; yL++; hBonds--;
                if (hBonds == 0) {
                    xR = xL = 0;
                    yR = yL = 0;
                }
                double tau = (-1.0/kTotal) * log(1.0 - r2);
                t += tau;
            } else if (randNum > kB1/kTotal && randNum <= (kB1 + kB2)/kTotal) {
                xR--; yR--; hBonds--;
                if (hBonds == 0) {
                    xR = xL = 0;
                    yR = yL = 0;
                }
                double tau = (-1.0/kTotal) * log(1.0 - r2);
                t += tau;
            } else if (randNum > (kB1 + kB2)/kTotal && randNum <= (kB1 + kB2 + kFL)/kTotal) {
                xL--; yL--; hBonds++;
                double tau = (-1.0/kTotal) * log(1.0 - r2);
                t += tau;
            } else if (randNum > (kB1 + kB2 + kFL)/kTotal && randNum <= 1.0) {
                xR++; yR++; hBonds++;
                double tau = (-1.0/kTotal) * log(1.0 - r2);
                t += tau;
            }
        }
        if (hBonds == 36) {
            printf("%.12f\n", t);
            t = 0;
            xR = xL = 0;
            yR = yL = 0;
            success++;
        }
        if(success == stopCondition) {
            // printf("done!\n");
            break;
        }
    }

    return 0;
}

double getEnergy(int j, int k, int jn, int kn) {
    double en = 0;
	// Nearest-neighbor free energies of SantaLucia et al. at 55C - 0.5M NaCl - unit of kT
    if ((k == 1 && kn == 1 && j == 2 && jn == 2) || (k == 2 && kn == 2 && j == 1 && jn == 1)) en = -0.89;
    else if (k == 1 && kn == 2 && j == 2 && jn == 1) en = -0.71;
    else if (k == 2 && kn == 1 && j == 1 && jn == 2) en = -0.21;
    else if ((k == 3 && kn == 1 && j == 4 && jn == 2) || (k == 2 && kn == 4 && j == 1 && jn == 3)) en = -1.63;
    else if ((k == 4 && kn == 2 && j == 3 && jn == 1) || (k == 1 && kn == 3 && j == 2 && jn == 4)) en = -1.63;
    else if ((k == 3 && kn == 2 && j == 4 && jn == 1) || (k == 1 && kn == 4 && j == 2 && jn == 3)) en = -1.39;
    else if ((k == 4 && kn == 1 && j == 3 && jn == 2) || (k == 2 && kn == 3 && j == 1 && jn == 4)) en = -1.40;
    else if (k == 3 && kn == 4 && j == 4 && jn == 3) en = -2.68;
    else if (k == 4 && kn == 3 && j == 3 && jn == 4) en = -2.89;
    else if ((k == 4 && kn == 4 && j == 3 && jn == 3) || (k == 3 && kn == 3 && j == 4 && jn == 4)) en = -2.34;
    
    // Repulsive nearest-neighbor free energies for stem-loop region at 55C - 0.5M NaCl - unit of kT
    if ((k == 11 && kn == 11 && j == 22 && jn == 22) || (k == 22 && kn == 22 && j == 11 && jn == 11)) en = 0.89;
    else if (k == 11 && kn == 22 && j == 22 && jn == 11) en = 0.71;
    else if (k == 22 && kn == 11 && j == 11 && jn == 22) en = 0.21;
    else if ((k == 33 && kn == 11 && j == 44 && jn == 22) || (k == 22 && kn == 44 && j == 11 && jn == 33)) en = 1.63;
    else if ((k == 44 && kn == 22 && j == 33 && jn == 11) || (k == 11 && kn == 33 && j == 22 && jn == 44)) en = 1.63;
    else if ((k == 33 && kn == 22 && j == 44 && jn == 11) || (k == 11 && kn == 44 && j == 22 && jn == 33)) en = 1.39;
    else if ((k == 44 && kn == 11 && j == 33 && jn == 22) || (k == 22 && kn == 33 && j == 11 && jn == 44)) en = 1.40;
    else if (k == 33 && kn == 44 && j == 44 && jn == 33) en = 2.68;
    else if (k == 44 && kn == 33 && j == 33 && jn == 44) en = 2.89;
    else if ((k == 44 && kn == 44 && j == 33 && jn == 33) || (k == 33 && kn == 33 && j == 44 && jn == 44)) en = 2.34;

    // Santa 1:
    else if ((k == 1  && kn == 1  && j == 2  && jn == 1) || (k == 1  && kn == 2  && j == 1  && jn == 1)) en = 1.23;
    else if ((k == 3  && kn == 1  && j == 4  && jn == 1) || (k == 1  && kn == 4  && j == 1  && jn == 3)) en = 0.95;
    else if ((k == 4  && kn == 1  && j == 3  && jn == 1) || (k == 1  && kn == 3  && j == 1  && jn == 4)) en = 0.67;
    else if ((k == 2  && kn == 1  && j == 1  && jn == 1) || (k == 1  && kn == 1  && j == 1  && jn == 2)) en = 0.93;
    else if ((k == 1  && kn == 3  && j == 2  && jn == 3) || (k == 3  && kn == 2  && j == 3  && jn == 1)) en = 2.58;
    else if ((k == 3  && kn == 3  && j == 4  && jn == 3) || (k == 3  && kn == 4  && j == 3  && jn == 3)) en = 1.60;
    else if ((k == 4  && kn == 3  && j == 3  && jn == 3) || (k == 3  && kn == 3  && j == 3  && jn == 4)) en = 1.29;
    else if ((k == 2  && kn == 3  && j == 1  && jn == 3) || (k == 3  && kn == 1  && j == 3  && jn == 2)) en = 1.35;
    else if ((k == 1  && kn == 4  && j == 2  && jn == 4) || (k == 4  && kn == 2  && j == 4  && jn == 1)) en = 0.17;
    else if ((k == 3  && kn == 4  && j == 4  && jn == 4) || (k == 4  && kn == 4  && j == 4  && jn == 3)) en = 0.35;
    else if ((k == 4  && kn == 4  && j == 3  && jn == 4) || (k == 4  && kn == 3  && j == 4  && jn == 4)) en = -1.23;
    else if ((k == 2  && kn == 4  && j == 1  && jn == 4) || (k == 4  && kn == 1  && j == 4  && jn == 2)) en = 0.85;
    else if ((k == 1  && kn == 2  && j == 2  && jn == 2) || (k == 2  && kn == 2  && j == 2  && jn == 1)) en = 1.57;
    else if ((k == 3  && kn == 2  && j == 4  && jn == 2) || (k == 2  && kn == 4  && j == 2  && jn == 3)) en = 0.45;
    else if ((k == 4  && kn == 2  && j == 3  && jn == 2) || (k == 2  && kn == 3  && j == 2  && jn == 4)) en = 1.08;
    else if ((k == 2  && kn == 2  && j == 1  && jn == 2) || (k == 2  && kn == 1  && j == 2  && jn == 2)) en = 1.31;
    // Santa 2:
    else if ((k == 1  && kn == 1  && j == 2  && jn == 3 ) || (k == 3  && kn == 2  && j == 1  && jn == 1)) en = 1.48;
    else if ((k == 1  && kn == 3  && j == 2  && jn == 1 ) || (k == 1  && kn == 2  && j == 3  && jn == 1)) en = 1.00;
    else if ((k == 3  && kn == 1  && j == 4  && jn == 3 ) || (k == 3  && kn == 4  && j == 1  && jn == 3)) en = 1.30;
    else if ((k == 3  && kn == 3  && j == 4  && jn == 1 ) || (k == 1  && kn == 4  && j == 3  && jn == 3)) en = 1.49;
    else if ((k == 4  && kn == 1  && j == 3  && jn == 3 ) || (k == 3  && kn == 3  && j == 1  && jn == 4)) en = 1.05;
    else if ((k == 4  && kn == 3  && j == 3  && jn == 1 ) || (k == 1  && kn == 3  && j == 3  && jn == 4)) en = 1.06;
    else if ((k == 2  && kn == 1  && j == 1  && jn == 3 ) || (k == 3  && kn == 1  && j == 1  && jn == 2)) en = 1.45;
    else if ((k == 2  && kn == 3  && j == 1  && jn == 1 ) || (k == 1  && kn == 1  && j == 3  && jn == 2)) en = 1.78;
    // Santa 3:
    else if ((k == 1  && kn == 1  && j == 2  && jn == 4 ) || (k == 4  && kn == 2  && j == 1  && jn == 1 )) en = 0.40;
    else if ((k == 1  && kn == 4  && j == 2  && jn == 1 ) || (k == 1  && kn == 2  && j == 4  && jn == 1 )) en = 0.23;
    else if ((k == 3  && kn == 1  && j == 4  && jn == 4 ) || (k == 4  && kn == 4  && j == 1  && jn == 3 )) en = 0.23;
    else if ((k == 3  && kn == 4  && j == 4  && jn == 1 ) || (k == 1  && kn == 4  && j == 4  && jn == 3 )) en = 0.70;
    else if ((k == 4  && kn == 1  && j == 3  && jn == 4 ) || (k == 4  && kn == 3  && j == 1  && jn == 4 )) en = -0.32;
    else if ((k == 4  && kn == 4  && j == 3  && jn == 1 ) || (k == 1  && kn == 3  && j == 4  && jn == 4 )) en = -0.79;
    else if ((k == 2  && kn == 1  && j == 1  && jn == 4 ) || (k == 4  && kn == 1  && j == 1  && jn == 2 )) en = 0.94;
    else if ((k == 2  && kn == 4  && j == 1  && jn == 1 ) || (k == 1  && kn == 1  && j == 4  && jn == 2 )) en = 1.11;
    // Santa 4:
    else if ((k == 1  && kn == 4  && j == 2  && jn == 2 ) || (k == 2  && kn == 2  && j == 4  && jn == 1 )) en = 1.33;
    else if ((k == 1  && kn == 2  && j == 2  && jn == 4 ) || (k == 4  && kn == 2  && j == 2  && jn == 1 )) en = 0.52;
    else if ((k == 3  && kn == 4  && j == 4  && jn == 2 ) || (k == 2  && kn == 4  && j == 4  && jn == 3 )) en = -0.30;
    else if ((k == 3  && kn == 2  && j == 4  && jn == 4 ) || (k == 4  && kn == 4  && j == 2  && jn == 3 )) en = -0.15;
    else if ((k == 4  && kn == 4  && j == 3  && jn == 2 ) || (k == 2  && kn == 3  && j == 4  && jn == 4 )) en = -0.05;
    else if ((k == 4  && kn == 4  && j == 2  && jn == 2 ) || (k == 2  && kn == 2  && j == 4  && jn == 4 )) en = 0.90;
    else if ((k == 4  && kn == 2  && j == 3  && jn == 4 ) || (k == 4  && kn == 3  && j == 2  && jn == 4 )) en = -0.47;
    else if (k == 4   && kn == 2  && j == 2  && jn == 4 ) en = 1.80;
    else if ((k == 2  && kn == 4  && j == 1  && jn == 2 ) || (k == 2  && kn == 1  && j == 4  && jn == 2)) en = 0.91;
    else if (k == 2   && kn == 4  && j == 4  && jn == 2 ) en = 1.21;
    else if ((k == 2  && kn == 2  && j == 1  && jn == 4 ) || (k == 4  && kn == 1  && j == 2  && jn == 2)) en = 0.88;
    // Santa 5:
    else if ((k == 1  && kn == 3  && j == 2  && jn == 2) || (k == 2  && kn == 2  && j == 3  && jn == 1 )) en = 1.21;
    else if ((k == 1  && kn == 2  && j == 2  && jn == 3) || (k == 3  && kn == 2  && j == 2  && jn == 1 )) en = 1.55;
    else if ((k == 3  && kn == 3  && j == 4  && jn == 2) || (k == 2  && kn == 4  && j == 3  && jn == 3 )) en = 1.28;
    else if ((k == 3  && kn == 2  && j == 4  && jn == 3) || (k == 3  && kn == 4  && j == 2  && jn == 3 )) en = 0.99;
    else if ((k == 4  && kn == 3  && j == 3  && jn == 2) || (k == 2  && kn == 3  && j == 3  && jn == 4 )) en = 1.03;
    else if ((k == 4  && kn == 2  && j == 3  && jn == 3) || (k == 3  && kn == 3  && j == 2  && jn == 4 )) en = 1.44;
    else if ((k == 2  && kn == 3  && j == 1  && jn == 2) || (k == 2  && kn == 1  && j == 3  && jn == 2 )) en = 1.78;
    else if ((k == 2  && kn == 2  && j == 1  && jn == 3) || (k == 3  && kn == 1  && j == 2  && jn == 2 )) en = 1.44;

    return en;
}

std::pair<int, int> getParams(int yL, int xL, int len) {
    int registry = yL - xL;
    int rMax, lMin;
    if (registry >= 0) {
        rMax = len - registry;
        lMin = 1;
    } else {
        rMax = len;
        lMin = 1 - registry;
    }
    return std::make_pair(rMax, lMin);
}
std::tuple<std::string, std::string, std::string, std::string> parseParams(int argc, char* argv[]) {
    std::string seq, num1, num2, stop;
    if (argc > 1) {
        for (int i = 0; i < argc; i++){
            std::string temp(argv[i]);
            if (temp == "--seq") seq = std::string (argv[i + 1]);
            if (temp == "--num1") num1 = std::string (argv[i + 1]);
            if (temp == "--num2") num2 = std::string (argv[i + 1]);
            if (temp == "--stop") stop = std::string (argv[i + 1]);
            // if (temp == "--runId") runId = std::string (argv[i + 1]);
        }
    }
    return std::make_tuple(seq, num1, num2, stop);
}
