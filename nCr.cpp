//
// Created by Elizalde Borja on 2023-07-22.
//
#include <iostream>
#include <vector>
#include <cmath>

int calculate_factorial(long, int);
const std::vector<int> vector_factors = {27, 11, 13, 37};
const std::vector<int> vector_primes = {3, 11, 13, 37};
const long N = 38;

int factorial_mod(long n, int index) {
    int prime = vector_primes[index];
    int factor = vector_factors[index];
    int result = 1;
    int j;
    for (int i = 1; i <= n; ++i) {
        j = i;
        while (j%prime == 0) {
            j = j/prime;
        }
        if ((result * (j%factor)) % factor == 0) {
            int a = 3;
        }
        result = (result * (j%factor)) % factor;

    }
    return result;
}

std::vector<long> get_factorial_simple(long n, int index) {
    std::vector<long> result(n + 1);
    result[0] = 1;
    //long j;
    for (int i = 1; i <= n; ++i) {
        if (i % vector_primes[index] != 0) {
            result[i] = (result[i - 1] * i) % vector_factors[index];
        } else {
            result[i] = result[i - 1];
        }
    //    j = i;
    //    while (j % vector_primes[index] == 0) {
    //        j = j / vector_primes[index];
    //    }
    //    result[i] = (result[i - 1] * j) % vector_factors[index];
    }

    return result;
}

std::vector<std::vector<long>> get_factorial_multiple(long n) {
    u_long factors_len = vector_factors.size();
    std::vector<std::vector<long>> res_matrix(factors_len, std::vector<long>(n+1));

    for (int i = 0; i <= factors_len - 1; i++) {
        std::vector<long> row_factorial_mod = get_factorial_simple(n, i);
        std::copy(row_factorial_mod.begin(), row_factorial_mod.end(), res_matrix[i].begin());
    }

    return res_matrix;
}

const std::vector<std::vector<long>> factorial_matrix = get_factorial_multiple(N);


int mod_inverse(int A, int M) {
    for (int X = 1; X < M; X++)
        if (((A % M) * (X % M)) % M == 1)
            return X;
    return 1;
}

int calculate_factorial(long n, int index) {
    int prime = vector_primes[index];
    int factor = vector_factors[index];
    if (n < prime) {
        return factorial_matrix[index][n];
    }
    else if (n < N) {
        return (factorial_matrix[index][n] * calculate_factorial(n/prime, index))%factor;
    } else {
        int factorial_rec = calculate_factorial(n/prime, index);
        //int factorial_sign = ((int) pow(factorial_matrix[index][factor], n/factor)) % factor;
        int factorial_sign = ((int) pow(-1, n/factor));
        int factorial_mod = factorial_matrix[index][n % factor];
        int out = (factorial_sign * factorial_mod * factorial_rec) % factor;
        if (out >= 0) {
            return out;
        } else {
            return factor + out;
        }
    }
}




int p_adic_valuation_factorial(long n, int p) {
    long L = static_cast<long>(std::log(n) / std::log(p));
    std::vector<int> values;

    for (int i = 1; i <= L; ++i) {
        values.push_back(n / static_cast<int>(std::pow(p, i)));
    }

    int sum_values = 0;
    for (int value : values) {
        sum_values += value;
    }

    return sum_values;
}
int get_mod_simple2(long n, int k, int i) {
    int num = factorial_matrix[i][n];
    int num_padic = p_adic_valuation_factorial(n, vector_primes[i]);
    int den1 = factorial_matrix[i][k];
    int den1_padic = p_adic_valuation_factorial(k, vector_primes[i]);
    int den2 = factorial_matrix[i][n-k];
    int den2_padic = p_adic_valuation_factorial(n-k, vector_primes[i]);

    int p_pow = num_padic - (den1_padic + den2_padic);
    int den_inv = mod_inverse((den1*den2) % vector_factors[i], vector_factors[i]);
    int non_inv_mod = pow(vector_primes[i], p_pow);
    int out = (non_inv_mod*num*den_inv) % vector_factors[i];
    return out;

}

int get_mod_simple(long n, int k, int i) {
    int num = calculate_factorial(n, i);
    int num_padic = p_adic_valuation_factorial(n, vector_primes[i]);
    int den1 = calculate_factorial(k, i);
    int den1_padic = p_adic_valuation_factorial(k, vector_primes[i]);
    int den2 = calculate_factorial(n-k, i);
    int den2_padic = p_adic_valuation_factorial(n-k, vector_primes[i]);

    int p_pow = num_padic - (den1_padic + den2_padic);
    int den_inv = mod_inverse((den1*den2) % vector_factors[i], vector_factors[i]);
    int non_inv_mod = pow(vector_primes[i], p_pow);
    non_inv_mod = non_inv_mod%vector_factors[i];
    int out = (non_inv_mod*num*den_inv) % vector_factors[i];
    return out;

}

std::vector<long> get_mod_multiple(long n, int k) {
    u_long factors_len = vector_factors.size();
    std::vector<long> vectors_mod(factors_len);
    for (int i = 0; i <= factors_len - 1; i++) {
        vectors_mod[i] = get_mod_simple(n, k, i);
    }
    return vectors_mod;
}

int CRT(std::vector<long> vectors_mod) {
    u_long factors_len = vectors_mod.size();
    // Compute product of all numbers
    int prod = 1;
    for (int i = 0; i < factors_len; i++)
        prod *= vector_factors[i];

    // Initialize result
    int result = 0;

    // Apply above formula
    for (int i = 0; i < factors_len; i++) {
        int pp = prod / vector_factors[i];
        result += vectors_mod[i] * mod_inverse(pp, vector_factors[i]) * pp;
    }

    return result % prod;
}

int solve(long n, int k) {
    std::vector<long> vectors_mod = get_mod_multiple(n, k);
    int nCk_mod = CRT(vectors_mod);
    return nCk_mod;

}

int main() {
    // 2147483647
    int n1 = 943832402;
    int k1 = 460466991;
    //95238
    //142857
    int out = solve(n1, k1);

    // Calculate res mods primes for all primes
        // Use factorials calculations
        // Calculate padic valuation of factorial
        // Find
    // Apply chinese res theorem to calculate final mod
    // for p in vector_primes:

    std::cout << "Factorial Array:\n";
    std::cout << "Factorial of " << n1 << " is: " << out << std::endl;

    return 0;
}