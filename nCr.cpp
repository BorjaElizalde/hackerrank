//
// Created by Elizalde Borja on 2023-07-22.
//
#include <iostream>
#include <vector>

std::vector<long> factorial_mod_p(long n, int p) {
    std::vector<long> result(n + 1);
    result[0] = 1;
    long j;
    for (int i = 1; i <= n; ++i) {
        j = i;
        while (j % p == 0) {
            j = j / p;
        }
        result[i] = (result[i - 1] * j) % p;
    }

    return result;
}

int find_base(int p_k, std::vector<int> vector_primes) {

    for (int p : vector_primes) {
        if (p_k % p == 0) {
            return p;
        }
    }

    // Return -1 if no prime factor is found
    return -1;
}


int main() {
    // 2147483647
    long n = 9999999;
    int p = 11;
    std::vector<int> vector_primes = {3, 11, 13, 37};
    p = find_base(p, vector_primes);
    std::vector<long> factorial_array = factorial_mod_p(n, p);

    std::cout << "Factorial Array:\n";
    std::cout << "Factorial of " << n << " is: " << factorial_array[n] << std::endl;

    return 0;
}