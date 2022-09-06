#include "Bignum.h"
#include <chrono>
#include <iostream>
#include <random>
#include <chrono>
#include <ctime>

using namespace std;

void bignum_Test()
{
    int M = 1000,
        T = 1000,
        good_count = 0;
    Bignum A,
        B,
        C,
        D;
    srand(time(NULL));

    while (T--)
    {
        int n = rand() % M + 1,
        m = rand() % M + 1;

        Bignum AA(n, 3),
        BB(m, 3);

        if (AA >= BB)
        {
            A = AA;
            B = BB;
        }
        else
        {
            B = AA;
            A = BB;
        }

        C = A / B;
        D = A % B;

        if (A == C * B + D && A - D == C * B && D < B)
        {
            good_count++;
        }
        
        // if (n >= m)
        // {
        //     if (A == C * B + D && A - D == C * B )
        //     {
        //         good_count++;
        //     }
        // }
        // else
        // {
        //     if (A == C * B + D && D < B)
        //     {
        //         good_count++;
        //     }
        // }
        cout << "T " << T << endl;
        cout << "n " << n << endl;
        cout << "m " << m << endl;
    }
    cout << good_count << endl;

    if (good_count == 1000)
    {
        cout << "Correct" << endl;
    }
    else
    {
        cout << "Incorrect" << endl;
    }
}

void fast_sq_Test()
{
    srand(time(NULL));
    int N = 300;
    while (N-- > 0) 
    {
        Bignum f(rand() % 1000 + 1, 3);
        Bignum fast_sq = f.fast_sq();
        Bignum mul = f * f;
        if (fast_sq != mul) 
        {
            cout << "Failed test with: " << endl;
            cout << "f_fast_sq: " << fast_sq << endl;
            cout << "  f_mul_f: " << mul << endl;
            cout << "f: " << f << endl << endl;
            break;
        }
    }
    cout << "Tests correct" << endl;

    int n = 50;
    int M = 5000;
    auto t1 = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n; i++) 
    {
        Bignum f(M, 3);
        Bignum f_sq = f.fast_sq();
    }
    auto t2 = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    std::cout << "fast_sq: time taken for " << n << " numbers to square: " << ms_int.count() << std::endl;

    t1 = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n; i++) 
    {
        Bignum f(M, 3);
        Bignum mul = f * f;
    }
    t2 = chrono::high_resolution_clock::now();
    ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    std::cout << "mul: time taken for " << n << " numbers to square: " << ms_int.count() << std::endl;

    // Bignum f;
    // f.input();
    // Bignum fast_sq = f.fast_sq();
    // Bignum mul = f * f;
    // fast_sq.print();
    // cout << endl;
    // mul.print();
    // cout << endl;
}

void fast_pow()
{
    srand(time(NULL));
    int N = 2;
    while (N-- > 0) 
    {
        Bignum f(1, 3);
        Bignum y(1, 3);
        Bignum pow = f.pow(y);
        Bignum st = f.stupid_pow(y);
        if (pow != st) 
        {
            cout << "Failed test with: " << endl;
            cout << "pow: " << pow << endl;
            cout << "  st: " << st << endl;
            cout << "f: " << f << endl << endl;
            cout << "  st: " << y << endl << endl;
            break;
        }
    }
    cout << "Tests correct" << endl;

    int n = 2;
    int M = 1;
    auto t1 = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n; i++) 
    {
        Bignum f(M, 3);
        Bignum y(M, 3);
        Bignum f_pow = f.pow(y);
    }
    auto t2 = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    std::cout << "pow: time taken for number to pow: " << ms_int.count() << std::endl;

    t1 = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n; i++) 
    {
        Bignum f(M, 3);
        Bignum y(M, 3);
        Bignum f_st = f.stupid_pow(y);
    }
    t2 = chrono::high_resolution_clock::now();
    ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    std::cout << "stupid_pow: time taken for number to pow: " << ms_int.count() << std::endl;

    // Bignum x;
    // Bignum y;
    // x.input();
    // y.input();
    // Bignum z = x.pow(y);
    // z.print();
    // cout << endl;
    // Bignum st = x.stupid_pow(y);
    // st.print();
    // cout << endl;
}

// void barret_z()
// {
//     srand(time(NULL));
//     int M = 2000;
//     int N = 100;
//     int d;
//     int c;
    
//     while (N-- > 0) 
//     {
//         int one = rand() % M + 1;
//         int two = rand() % M + 1;
//         if (one >= two)
//         {
//             c = one;
//             d = two;
//         }
//         else
//         {
//             c = two;
//             d = one;
//         }
//         Bignum a(d, 3);
//         Bignum m(c, 3);
//         Bignum r = a % m;
//         Bignum r_fast = a.barret(m, m.get_z(m));
//         if (r != r_fast) 
//         {
//             cout << "Failed test with: " << endl;
//             cout << "        a: " << a << endl;
//             cout << "        m: " << m << endl;
//             cout << " a%m: " << r << endl;
//             cout << " a%m fast: " << r_fast << endl;
//             break;
//         }
//     }
//     cout << "Tests correct" << endl;

//     N = 2000;
//     Bignum m(800, 3);
//     Bignum z = m.get_z(m);
//     auto t1 = chrono::high_resolution_clock::now();
//     for (size_t i = 0; i < N; i++) 
//     {
//         Bignum x(1000, 3);
//         try{
//             Bignum res_fast = x.barret(m, z);
//         } catch (std::invalid_argument& e) { }
//     }
//     auto t2 = chrono::high_resolution_clock::now();
//     auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
//     std::cout << "barret: time taken for " << N << " numbers to mod: " << ms_int.count() << std::endl;

//     t1 = chrono::high_resolution_clock::now();
//     for (size_t i = 0; i < N; i++) 
//     {
//         Bignum x(1000, 3);
//         Bignum res = x % m;
//     }
//     t2 = chrono::high_resolution_clock::now();
//     ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
//     std::cout << "a%m: time taken for " << N << " numbers to mod: " << ms_int.count() << std::endl;

//     // Bignum x;
//     // Bignum m;
//     // x.input();
//     // m.input();
//     // Bignum z = m.get_z(m);
//     // Bignum r = x.barret(m, z);
//     // r.print();
//     // cout << endl;

//     // Bignum slow = x % m;
//     // slow.print();
//     // cout << endl;
// }

void barret_z()
{
    const ssize_t M = 1000;
    ssize_t N = 100;

    thread_local static mt19937 mt(static_cast<uint32_t>(time(nullptr)));
    while (N --> 0) {
        size_t one = mt() % M + 1;
        size_t two = mt() % M + 1;
        Bignum a(one, 3);
        Bignum b(two, 3);

        Bignum q = a / b;
        Bignum r = a % b;
        try {
            Bignum fast_r = a.barret(b, b.get_z(b));

            if (r != fast_r) {
                cout << "Failed barret_mod with values: " << endl;
                cout << "        a: " << a << endl;
                cout << "        b: " << b << endl;
                cout << " a%b slow: " << r << endl;
                cout << " a%b fast: " << fast_r << endl;
            }

        } catch (std::invalid_argument& e) {
        }

        if (a != q * b + r) {
            cout << "Failed [a == c * b + d] with values: " << endl;
            cout << "        a: " << a << endl;
            cout << "        b: " << b << endl;
            cout << "q = a / b: " << q << endl;
            cout << "r = a % b: " << r << endl;
        }

        if (a - r != q * b) {
            cout << "Failed [a - d == c * b] with values: " << endl;
            cout << "        a: " << a << endl;
            cout << "        b: " << b << endl;
            cout << "q = a / b: " << q << endl;
            cout << "r = a % b: " << r << endl;
        }

        if (r >= b) {
            cout << "Failed [d < b] with values: " << endl;
            cout << "        a: " << a << endl;
            cout << "        b: " << b << endl;
            cout << "q = a / b: " << q << endl;
            cout << "r = a % b: " << r << endl;
        }
    }
    cout << "Tests correct" << endl;

    { 
        const ssize_t M = 1000;
        ssize_t N = 2000;

        Bignum m(M - (M/3), 3);
        Bignum z = m.get_z(m);

        auto t1 = chrono::high_resolution_clock::now();
        for (size_t i = 0; i < N; i++) 
        {
            Bignum x(M, 3);
            try{
                Bignum res_fast = x.barret(m, z);
            } catch (std::invalid_argument& e) { }
        }
        auto t2 = chrono::high_resolution_clock::now();
        auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
        std::cout << "barret: time taken for " << N << " numbers to mod: " << ms_int.count() << std::endl;

        t1 = chrono::high_resolution_clock::now();
        for (size_t i = 0; i < N; i++) 
        {
            Bignum x(M, 3);
            Bignum res = x % m;
        }
        t2 = chrono::high_resolution_clock::now();
        ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
        std::cout << "a%m: time taken for " << N << " numbers to mod: " << ms_int.count() << std::endl;
    }

    {
        // Bignum x;
        // Bignum m;
        // x.input();
        // m.input();
        // Bignum z = m.get_z(m);
        // Bignum r = x.barret(m, z);
        // r.print();
        // cout << endl;

        // Bignum slow = x % m;
        // slow.print();
        // cout << endl;
    }
}

int main()
{
    // fast_sq_Test();
    // fast_pow();
    // barret_z();

    // int t = 10;
    // for (int i = 0; i < 6; i++)
    // {
    // 	Bignum rnd(5, 3);
    // 	rnd.print();
    //     cout << endl;
    // 	if (rnd.ferma_test(t))
    // 	{
    // 		cout << "\tTest Ferma: prime number" << endl;
    // 	}
    // 	else
    // 	{
    // 		cout << "\tTest Ferma: composite number" << endl;
    // 	}
    // 	if (rnd.solovay_strassen(t))
    // 	{
    // 		cout << "\tSolovay_strassen: prime number" << endl;
    // 	}
    // 	else
    // 	{
    // 		cout << "\tSolovay_strassen: composite number" << endl;
    // 	}
    // 	// cout << "\n" << endl;
    // }
    // cout << endl;

    // Bignum number = BASE(565);
    // size_t attempt = 0;
    // while (number.ferma_test(3) == false) 
    // { 
    //     attempt++; 
    // }
    // std::cout << "Ferma prime test lied on " << attempt << " attempt" << std::endl;

    // attempt = 0;
    // while (number.solovay_strassen(2) == false) 
    // { 
    //     attempt++; 
    // }
    // std::cout << "Solovay Strassen lied on " << attempt << " attempt" << std::endl;

    // number = BASE(25889); // is prime number
    // number.print();
    // std::cout << (number.ferma_test(400) == 1 ? "\tPrime number to Ferma" : "\tComposite number to Ferma") << std::endl;
    // std::cout << (number.solovay_strassen(17) == 1 ? "\tPrime number to Solovay" : "\tComposite number to Solovay") << std::endl;

    // number = DBASE(9381); // is not prime number
    // number.print();
    // std::cout << (number.ferma_test(1) == 1 ? "\tPrime number to Ferma" : "\tComposite number to Ferma") << std::endl;
    // std::cout << (number.solovay_strassen(10) == 1 ? "\tPrime number to Solovay" : "\tComposite number to Solovay") << std::endl;

    // number = BASE(561);
    // number.print();
    // std::cout << (number.ferma_test(1) == 1 ? "\tPrime number to Ferma" : "\tComposite number to Ferma") << std::endl;
    // std::cout << (number.solovay_strassen(15) == 1 ? "\tPrime number to Solovay" : "\tComposite number to Solovay") << std::endl;

    // Bignum number;
    // number.input();
    // number.print();
    // cout << endl;
    // std::cout << (number.ferma_test(1) == 1 ? "\tPrime number to Ferma" : "\tComposite number to Ferma") << std::endl;
    // std::cout << (number.solovay_strassen(100) == 1 ? "\tPrime number to Solovay" : "\tComposite number to Solovay") << std::endl;
    // Bignum fi = number.phi();
    // fi.print();
    // cout << endl;


    Bignum strong_prime = strong_prime.gen_strong(5);
    // std::cout << "Generated strong prime: " << strong_prime << std::endl;
    std::cout << "Generated strong prime: ";
    strong_prime.print();
    cout << endl;
    std::cout << (strong_prime.solovay_strassen(100) == 1 ? "Prime number to Solovay" : "Composite number to Solovay") << std::endl;

    return 0;
}
