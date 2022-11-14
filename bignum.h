#ifndef BIG_NUM
#define BIG_NUM
#include <iostream>
#include <random>
#include <cstddef>
#include <ctime>
#include <cstdint>
#include <string>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sstream>
#include <map>

using namespace std;

static std::mt19937& mt() 
{
    // initialize once per thread
    thread_local static std::mt19937 mt(static_cast<uint32_t>(time(nullptr)));
    return mt;
}

typedef unsigned short BASE;
typedef unsigned long long DBASE;
typedef long long DBASE_2;
#define BASE_SIZE (sizeof(BASE) * 8)    //размер типа в битах
const DBASE base = (2LL << (BASE_SIZE - 1));

class Bignum{
private:
    BASE *coef; 
    int len;
    int max_len;
public:
    Bignum();   //конструктор по умолчанию              +
    Bignum(int , int ); //конструктор с параметрами     +
    void len_norm(); //                                 +
    Bignum (const Bignum &); //конструктор копирования  +
    Bignum& operator=(const Bignum &); //               +
    ~Bignum();//                                        +
    Bignum(BASE );//                                    +
    Bignum& operator=(const BASE&); 

    bool operator==(const Bignum&) const;//+
    bool operator!=(const Bignum&);//+
    bool operator>(const Bignum&) const;// +
    bool operator<(const Bignum&) const;// +
    bool operator>=(const Bignum&);//+
    bool operator<=(const Bignum&);//+

    Bignum operator+(const Bignum&);// +
    Bignum operator+=(const Bignum&);//+
    Bignum operator-(const Bignum&);// +
    Bignum operator-=(const Bignum&);//+
    Bignum operator+(BASE );//         +
    Bignum operator*(BASE );//         +
    Bignum operator*(const Bignum&);// +
    Bignum operator*=(const Bignum&);//+
    Bignum operator/(BASE );//         +
    Bignum operator%(BASE );//           +
    Bignum operator/(Bignum&);  
    Bignum operator%(Bignum&);
    Bignum& operator/=(BASE);
    Bignum& operator/=(Bignum&);
    Bignum input();//                  +
    void print();//                    +

    friend istream& operator>>(istream &in, Bignum &other); //16-ый ввод   +
    friend ostream& operator<<(ostream &out, Bignum &other); //16-ый вывод +
    Bignum& size(size_t );
    Bignum fast_sq();
    Bignum pow(const Bignum& n);
    Bignum pow_m(const Bignum& pow, Bignum& mod);
    Bignum stupid_pow(const Bignum& n);
    Bignum barret(const Bignum& m, const Bignum& z);
    Bignum get_z(Bignum& mod);
    Bignum shiftr(size_t k);
    Bignum shiftl(size_t k);
    bool is_even();
    bool is_odd();
    Bignum phi();
    bool ferma_test(size_t t = 100);
    bool solovay_strassen(size_t rounds);
    Bignum gen_strong(size_t half_size);
    Bignum sqrt_num(Bignum num);
    Bignum trial_div(Bignum num);
    void trial_div_method(vector <Bignum> &p, Bignum n);
    Bignum alway(Bignum n);
    void method_ferma(pair<Bignum, Bignum> &res, Bignum num);
    Bignum method_p_pollard(Bignum num);
    Bignum Gelfond(Bignum g, Bignum n, Bignum a);
};

#endif
