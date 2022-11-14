#include "bignum.h"

using namespace std;

Bignum:: Bignum()
{
    max_len = 1;
    len = max_len;
    coef = new BASE [max_len];
    coef[0] = 0;
}

Bignum:: Bignum(int max_len_, int key)
{
    max_len = max_len_;
    len = max_len;
    coef = new BASE [max_len];

    switch(key)
    {
        case 2:
        {
            len = 1;
            for(int i = 0; i < max_len; i++)
            {
                coef[i] = 0;
            }
            break;
        }
        case 3:
        {
            // srand(time(NULL));
            // for(int i = 0; i < max_len; i++)
            // {
            //     coef[i] = rand();
            // }
            thread_local static std::mt19937 mersenne = mt();

            for (size_t i = 0; i < max_len; i++) 
            {
                coef[i] = mersenne();
            }
            len_norm();
            break;
        }
        default:
            cout << "Please enter the correct number 2 or 3" << endl;
    }
}

void Bignum::len_norm()
{
    // len = max_len;
    while ((len > 1) && (coef[len - 1] == 0))
        len--;
}

Bignum:: Bignum (const Bignum& other)
{
    len = other.len;
    max_len = other.max_len;
    coef = new BASE[len];

    for (int i = 0; i < len; i++)
    {
        coef[i] = other.coef[i];
    }
}

Bignum& Bignum::operator=(const BASE& bn) 
{
    len = 1;
    coef[0] = bn;
    return *this;
}

Bignum& Bignum:: operator=(const Bignum &other)
{
    if(this != &other)  //other=other
    {
        delete[] coef;
        len = other.len;
        max_len = other.max_len;
        coef = new BASE [len]; 

        for (size_t i = 0; i < len; ++i) 
        {
            coef[i] = other.coef[i];
        }
    }
    return *this;
}

Bignum:: ~Bignum()
{
    if(coef != NULL)
    { 
        delete [] coef; 
        coef = NULL;
        len = 0;
        max_len = 0;
    } 
}

bool Bignum:: operator==(const Bignum &other) const
{
    if (len != other.len)
        return false;

    for (size_t i = 0; i < len; i++)
    {
        if (coef[i] != other.coef[i])
            return false;
    }

    return true;
}

bool Bignum:: operator!=(const Bignum &other)
{
    return !(*this == other);
}

bool Bignum:: operator>(const Bignum &other) const
{
    if (len > other.len)
        return true;

    if (len < other.len)
        return false;

    for (ssize_t i = ssize_t(len) - 1; i >= 0; i--)
    {
        if (coef[i] > other.coef[i])
            return true;
        if (coef[i] < other.coef[i])
            return false;
    }

    return false;
}

bool Bignum:: operator<(const Bignum &other) const
{
    if (len > other.len)
        return false;

    if (len < other.len)
        return true;

    for (ssize_t i = ssize_t(len) - 1; i >= 0; i--)
    {
        if (coef[i] > other.coef[i])
            return false;
        if (coef[i] < other.coef[i])
            return true;
    }

    return false;   
}

bool Bignum:: operator>=(const Bignum &other)
{
    return !(*this < other);
}

bool Bignum:: operator<=(const Bignum &other)
{
    return !(*this > other);
}

bool Bignum_check(const string input)
{
    for (char i : input)
    {
        if( '0' > i || (i > '9' && 'A' > i)) 
        {
            cout << i << '\n';
            return false;
        }
        else if( 'f' < i || (i < 'a' && 'F' < i)) 
        {
            cout << i << '\n';
            return false;
        }
    }

    return true;
}

istream& operator>>(istream &in, Bignum &other) //16-ый ввод
{
    string input;
    in >> input;
    if(!Bignum_check(input))
    {
        printf("\nIncorrect line entered\n");
        exit(0);
    }

    other.len = (input.length() - 1) / (BASE_SIZE / 4) + 1;
    other.max_len = other.len;
    if (other.coef != NULL)
    {
        delete[] other.coef;
        other.coef = NULL;
    }
    other.coef = new BASE[other.max_len]();
    int temp = 0;
    int j = 0;

    for (int i = input.length() - 1, k = 0; i >= 0; i--, k++) 
    {
        if(k >= BASE_SIZE / 4)
        {
            k = 0;
            j++;
        }
        if( '0' <= input[i] && input[i] <= '9') 
        {
            temp = input[i] - '0';
        }
        if('a' <= input[i] && input[i] <= 'f') 
        {
            temp = input[i] - 'a' + 10;
        }
        if('A' <= input[i] && input[i] <= 'F') 
        {
            temp = input[i] - 'A' + 10;
        }
        other.coef[j] |= (temp << (k * 4));
    }
    other.len_norm();

    return in;
}

ostream& operator<<(ostream &out, Bignum &other) //16-ый вывод
{
    for(int i = other.len - 1; i >= 0; --i)
    {
        cout.width(BASE_SIZE/4);
        cout.fill('0');
        cout << hex << other.coef[i];
    }

    return out; 
}

Bignum Bignum::operator+(BASE b)
{
    DBASE x = b;
    int j = 0;
    Bignum res(len + 1, 2);
    res.len = len + 1;
    for (; j<len; j++){
        res.coef[j]= x+=coef[j];
        x = x >> BASE_SIZE;
    }
    res.coef[j] = x;
    res.len_norm();

    return res;
}

Bignum Bignum:: operator+(const Bignum &other)
{
    Bignum res(max(len, other.len) + 1, 2);
    res.len = max(len, other.len) + 1;
    int j;
    BASE k = 0;
    DBASE temp = 0;

    for (j = 0; j < min(len, other.len); j++)
    {
        temp = DBASE(coef[j]) + DBASE(other.coef[j]) + DBASE(k); 
        res.coef[j] = (BASE)temp; // остаток от деления
        k = temp >> BASE_SIZE; //частное 
    }

    if (j >= len)
    {
        while (j<other.len)
        {
            temp = DBASE(other.coef[j]) + DBASE(k);
            res.coef[j] = temp;
            k = temp >> BASE_SIZE;
            j++;
        }
    }
    else
    {
        while (j<len)
        {
            temp = DBASE(coef[j]) + DBASE(k);
            res.coef[j] = temp;
            k = temp >> BASE_SIZE;
            j++;
        }
    }

    if (k==0)
    {
        res.max_len = res.len = j;
    }
    else
    {
        res.coef[j] = k;
        res.max_len = res.len = j + 1;
    }

    return res;
}

Bignum Bignum:: operator+=(const Bignum &other)
{
    Bignum res;
    res = *this + other;
    *this = res;

    return *this;
}

Bignum Bignum:: operator-(const Bignum &other)
{
    if(*this < other)
    {
        throw std::invalid_argument("first value should be bigger than second to subtract");
    }

    int j = 0;
    BASE k = 0;
    DBASE temp = 0;
    Bignum res(len, 2);
    res.len = len;

    while (j < other.len)
    {
        temp = (DBASE(1) << BASE_SIZE) + DBASE(coef[j]);  // temp = 1*b+coef[j]
        temp = temp - DBASE(other.coef[j]) - DBASE(k);
        res.coef[j] = temp;
        k = !(temp >> BASE_SIZE);
        j++;
    }
    while (j < len)
    {
        temp = (DBASE(1) << BASE_SIZE) + DBASE(coef[j]);
        temp -= DBASE(k);
        res.coef[j] = temp;
        k = !(temp >> BASE_SIZE);
        j++;
    }
    res.len_norm();

    return res;
}

Bignum Bignum:: operator-=(const Bignum &other)
{
    Bignum res;
    res = *this - other;
    *this = res;

    return *this;
}

//умножение на цифру
Bignum Bignum:: operator*(BASE num)
{   
    int j = 0;
    BASE k = 0;
    DBASE temp = 0;
    Bignum res(len + 1, 2);
    res.len = len + 1;
    
    while (j < len)
    {
        temp = DBASE(coef[j]) * DBASE(num) + DBASE(k);
        res.coef[j] = (BASE)temp;
        k = temp >> BASE_SIZE;
        j++;
    }
    res.coef[j] = k;
    res.len_norm();

    return res;
}

Bignum Bignum:: operator*(const Bignum &other)
{
    int i;     // по первому числу
    int j = 0; // по второму

    Bignum res (len + other.len, 2);
    res.len = len + other.len;
    BASE k;
    DBASE temp = 0;

    while (j < other.len)
    {
        //если коэффициент 2-ого множителя равен нулю, то пропуск
        if (other.coef[j] == 0)
        {
            j++;
            continue;
        }
        i = 0;
        k = 0;

        while (i < len)
        {
            temp = DBASE(coef[i]) * DBASE(other.coef[j]) + DBASE(res.coef[i + j]) + k;
            res.coef[i + j] = temp;
            k = temp >> BASE_SIZE;
            i++;
        }
        res.coef[j + len] = k;
        j++;
    }
    res.len_norm();

    return res;
}

Bignum Bignum:: operator*=(const Bignum &other)
{
    Bignum res;
    res = *this * other;
    *this = res;
    return *this;
}

Bignum Bignum:: operator/(BASE num)
{
    if (num == 0) 
    {
        throw std::overflow_error("can't divide by zero");
    }
    int j = 0;
    BASE r = 0;
    DBASE temp = 0;
    Bignum res(len, 2);
    res.len = len;

    while (j < len)
    {
        temp = (DBASE(r) << BASE_SIZE) + DBASE(coef[len - 1 - j]);
        res.coef[len - 1 - j] = temp / DBASE(num);                 
        r = temp % DBASE(num);
        j++;
    }
    res.len_norm();  

    return res;
}

Bignum Bignum:: operator%(BASE other)
{
    if (other == 0) 
    {
        throw std::overflow_error("can't divide by zero");
    }
    int j = 0;
    DBASE tmp, r = 0;

    while (j < len)
    {
        tmp = (DBASE(r) << BASE_SIZE) + coef[len - 1 - j];
        r = tmp % other;
        j++;
    }

    return Bignum(r);
}

void Bignum:: print()
{
    if (*this == Bignum(0)) 
    {
        cout << "0";
        return;
    }
    string res;
    Bignum tmp = *this, b;
    int i;
    BASE num = 0;
    BASE cc = 10;

    while(!(tmp.coef[0] == 0 && tmp.len == 1))
    {
        Bignum temp = tmp % cc;
        res = to_string(temp.coef[0]) + res;
        tmp = tmp / cc;
    }
    // reverse(res.begin(), res.end());
    cout << res;
}

Bignum Bignum::input()
{
    string dec;
    cin >> dec;

    size_t len_s = dec.length();
    if (coef != nullptr)
    {
        delete[] coef;
        coef = nullptr;
    }
    max_len = len = 1;
    coef = new BASE[max_len];
    coef[0] = 0;

    size_t j = 0;
    BASE basis = 10;

    while (j < len_s)
    {
        *this = *this * Bignum(basis) + BASE(dec[j] - '0');
        j++;
    }

    return *this;
}

Bignum:: Bignum (BASE num)
{
    max_len = 1;
    len = max_len;
    coef = new BASE [max_len];
    coef[0] = num;
}

Bignum& Bignum::size(size_t new_max_len) 
{
    auto* new_coef = new BASE[new_max_len];

    size_t new_size = (new_max_len > len) ? len : new_max_len;
    for (size_t i = 0; i < new_size; i++) 
    {
        new_coef[i] = coef[i];
    }

    for (size_t i = new_size; i < new_max_len; i++) 
    {
        new_coef[i] = 0;
    }

    this->~Bignum();
    coef = new_coef;
    max_len = new_max_len;
    len = new_size;

    return *this;
}

Bignum Bignum::operator/(Bignum& bn)
{
    if (*this < bn)
    {
        return Bignum(0);
    }
    if (bn.len == 1)
    {
        return *this / bn.coef[0];
    }
    if (*this == bn) 
    {
        return Bignum(1);
    }

    Bignum res(len - bn.len + 1, 2);
    res.len = len - bn.len + 1;

    BASE d = base / (bn.coef[bn.len - 1] + 1);
    Bignum u = *this * d;
    Bignum v = bn * d;

    if (this->len == u.len) 
    {
        if (u.len == u.max_len) 
        {
            u.size(u.max_len + 1);
        }
        u.coef[this->len] = 0;
        u.len = this->len + 1;
    }

    for (auto j = ssize_t(len - bn.len); j >= 0; j--) 
    {
        DBASE q = ((DBASE(u.coef[j + v.len]) << BASE_SIZE) + u.coef[j + v.len - 1]) / v.coef[v.len - 1];
        DBASE r = ((DBASE(u.coef[j + v.len]) << BASE_SIZE) + u.coef[j + v.len - 1]) % v.coef[v.len - 1];

        if (q == base || q * v.coef[v.len - 2] > base * r + u.coef[j + v.len - 2]) 
        {
            q -= 1;
            r += v.coef[v.len - 1];

            if (r < base && (q == base || q * v.coef[v.len - 2] > base * r + u.coef[j + v.len - 2])) 
            {
                q -= 1;
            }
        }

        Bignum sub = Bignum(q) * v;

        DBASE_2 tmp;
        uint8_t k = 0;
        size_t min_size = (bn.len + 1 > sub.len) ? (sub.len) : (bn.len + 1);
        for (size_t i = 0; i < min_size; i++) 
        {
            tmp = DBASE_2(u.coef[j + i]) - sub.coef[i] - k;
            u.coef[j + i] = tmp; // % base

            k = 0;
            if (tmp < 0) 
            {
                k = 1;
            }
        }

        for (size_t i = min_size; k && i < bn.len + 1; i++) 
        {
            tmp = DBASE_2(u.coef[j + i]) - k;
            u.coef[j + i] = tmp; // % base

            k = 0;
            if (tmp < 0) 
            {
                k = 1;
            }
        }

        if (k != 0) 
        {
            q -= 1;
            BASE k = 0;
            for (size_t i = 0; i < v.len; i++) 
            {
                tmp = u.coef[j + i] + v.coef[i] + k;
                u.coef[j + i] = tmp; // % base
                k = tmp >> BASE_SIZE;
            }
            u.coef[v.len + j] += k;
        }

        res.coef[j] = q;
    }
    len_norm();

    return res;
}


Bignum Bignum::operator%(Bignum& bn)
{
    if (*this < bn) 
    {
        return *this;
    }
    if (*this == bn) 
    {
        return Bignum(0);
    }
    if (bn.len == 1) 
    {
        return *this % bn.coef[0];
    }

    BASE d = base / (bn.coef[bn.len - 1] + 1);
    Bignum u = *this * d;
    Bignum v = bn * d;
    if (this->len == u.len) 
    {
        if (u.len == u.max_len) 
        {
            u.size(u.max_len + 1);
        }
        u.coef[this->len] = 0;
        u.len = this->len + 1;
    }

    for (auto j = ssize_t(this->len - bn.len); j >= 0; j--) 
    {
        DBASE q = ((DBASE(u.coef[j + v.len]) << BASE_SIZE) + u.coef[j + v.len - 1]) / v.coef[v.len - 1];
        DBASE r = ((DBASE(u.coef[j + v.len]) << BASE_SIZE) + u.coef[j + v.len - 1]) % v.coef[v.len - 1];

        if (q == base || q * v.coef[v.len - 2] > base * r + u.coef[j + v.len - 2]) 
        {
            q -= 1;
            r += v.coef[v.len - 1];

            if (r < base && (q == base || q * v.coef[v.len - 2] > base * r + u.coef[j + v.len - 2])) 
            {
                q -= 1;
            }
        }

        Bignum sub = Bignum(q) * v;

        DBASE_2 tmp;
        uint8_t k = 0;
        size_t min_size = (bn.len + 1 > sub.len) ? (sub.len) : (bn.len + 1);
        for (size_t i = 0; i < min_size; i++) 
        {
            tmp = DBASE_2(u.coef[j + i]) - sub.coef[i] - k;
            u.coef[j + i] = tmp; // % base

            k = 0;
            if (tmp < 0) 
            {
                k = 1;
            }
        }

        for (size_t i = min_size; k && i < bn.len + 1; i++) 
        {
            tmp = DBASE_2(u.coef[j + i]) - k;
            u.coef[j + i] = tmp; // % base

            k = 0;
            if (tmp < 0) 
            {
                k = 1;
            }
        }

        if (k != 0) 
        {
            BASE k = 0;
            for (size_t i = 0; i < v.len; i++) 
            {
                tmp = u.coef[j + i] + v.coef[i] + k;
                u.coef[j + i] = tmp; // % base
                k = tmp >> BASE_SIZE;
            }
            u.coef[v.len + j] += k;
        }
    }

    Bignum r = u / d;
    return r;
}

Bignum& Bignum::operator/=(Bignum& bn) 
{
    *this = *this / bn;
    return *this;
}

Bignum& Bignum::operator/=(BASE n) 
{
    *this = *this / n;
    return *this;
}

Bignum Bignum::fast_sq() 
{
	Bignum res (2 * len + 1, 2);
    res.len = 2 * len + 1;
    DBASE tmp = 0;
    for (size_t i = 0; i < len; i++)
    {
        tmp = DBASE(res.coef[2 * i]) + DBASE(coef[i]) * DBASE(coef[i]); //uv
        res.coef[2 * i] = tmp; // v

        for (size_t j = i + 1; j < len; j++)
        {
            tmp = DBASE(res.coef[i + j]) + DBASE(2) * DBASE(coef[i]) * DBASE(coef[j]) + DBASE(tmp >> BASE_SIZE);
            res.coef[i + j] = tmp; 
        }
        tmp >>= BASE_SIZE;
        tmp += DBASE(res.coef[i + len]);
        tmp += DBASE(res.coef[i + len + 1]) << BASE_SIZE;

        res.coef[i + len] = tmp;
        res.coef[i + len + 1] = tmp >> BASE_SIZE;
    }
    res.len_norm();
    return res;
}

Bignum Bignum:: pow(const Bignum& y)
{
    if (y.coef[0] == 0)
    {
        return  Bignum(BASE(1));
    }

    int n = (y.len - 1) * BASE_SIZE;
    int idx = BASE_SIZE - 1;
    while (true)
    {
        if (y.coef[y.len - 1] & (1 << idx))
        {
            n += idx + 1;
            break;
        }
        idx--;
    }

    Bignum z(1);
    z = *this;
    for (int i = n - 2; i >= 0; i--)
    {
        z = z.fast_sq();
        if (y.coef[i / BASE_SIZE] & (1 << (i % BASE_SIZE)))
        {
            z = z * (*this);
        }
    }
    return z;
}

Bignum Bignum:: pow_m(const Bignum& pow, Bignum& mod)
{
    	if (mod.coef[0] == 1 && mod.len == 1)
		{
			Bignum res;
			return res;
		}

		int n = (pow.len - 1) * BASE_SIZE;
		int idx = BASE_SIZE - 1;
		while (true)
		{
			if (pow.coef[pow.len - 1] & (1 << idx))
			{
				n += idx + 1;
				break;
			}
			idx--;
		}

		Bignum z(1);
		z = (*this) % mod;
		for (int i = n - 2; i >= 0; i--)
		{
			z = z.fast_sq() % mod;
			if (pow.coef[i / BASE_SIZE] & (1 << (i % BASE_SIZE)))
			{
				z = (z * (*this)) % mod;
			}
		}
		return z;
}

Bignum Bignum::stupid_pow(const Bignum& n) 
{
    Bignum i((BASE) 0);
    Bignum res((BASE) 1);
    while (i < n) {
        res *= *this;
        i += 1;
    }

    return res;
}

Bignum Bignum::shiftl(size_t k)
{
    size_t new_len = k + len;
    Bignum res(new_len, 2);
    res.len = new_len;
    for (size_t i = k; i < new_len; i++) 
    {
        res.coef[i] = coef[i - k];
    }

    return res;
}

Bignum Bignum::shiftr(size_t k)
{
    size_t new_size = (k > len) ? 0 : len - k;
    Bignum res(new_size, 2);
    res.len = new_size;
    for (size_t i = 0; i < new_size; i++) 
    {
        res.coef[i] = coef[i + k];
    }

    return res;
}

Bignum Bignum::barret(const Bignum& m, const Bignum& z)
{
    if (len > 2 * m.len) 
    {
        throw std::invalid_argument("number is to big for this method to work!");
    }

    Bignum q = shiftr(m.len - 1) * z;
    q = q.shiftr(m.len + 1);

    Bignum r1 = *this;
    if (r1.len > m.len) 
    {
        r1.len = m.len;
    }

    Bignum r2 = q * m;
    if (r2.len > m.len) 
    {
        r2.len = m.len;
    }

    Bignum res(m.len + 2, 2);
    if (r1 >= r2) 
    {
        res = r1 - r2;
    } 
    else 
    {
        res.coef[m.len + 1] = 1;
        res.len = m.len + 2;        //зам
        res += r1 - r2;
    }

    while (res >= m) 
    {
        res -= m;
    }

    return res;
}

Bignum Bignum:: get_z(Bignum& mod) 
{
    Bignum z = Bignum(2 * mod.len + 1, 2);
    z.coef[2 * mod.len] = 1;
    z.len = 2 * mod.len + 1;
    z = z / mod;
    return z;
}

bool Bignum::is_even() 
{
    return coef[0] % 2 == 0;
}

Bignum Bignum::phi() 
{
    Bignum result = *this;
    Bignum n = *this;
    for (Bignum i = Bignum(BASE(2)); i * i <= n; i = i + BASE(1))
    {
        if (n % i == 0) 
        {
            while (n % i == 0)
            {
                n = n / i;
            }
            result -= result / i;
        }
    }
    if (n > 1)
    {
        result -= result / n;
    }

    return result;
}

bool Bignum::ferma_test(size_t t)
{
    if (*this == 3 || *this == 2) 
    {
        return true;
    }

    if (coef[0] % 2 == 0) 
    {
        return false;
    }

    while (t > 0) 
    {
        Bignum a = Bignum(len, 3);
        if (Bignum(BASE(2)) > a || a > *this - 2) 
        {
            continue;
        }

        Bignum r = a.pow_m(*this - 1, *this);
        if (r != 1) 
        {
            return false;
        }

        t -= 1;
    }

    return true;
}

int8_t jacobi(Bignum a, Bignum n) 
{
    if (a == 0) 
    {
        return 0;
    }
    if (a == 1) 
    {
        return 1;
    }

    Bignum b = a;
    Bignum k = 0;
    while (b.is_even()) 
    {
        b /= 2;
        k += 1;
    }

    int8_t s;
    if (k.is_even()) 
    {
        s = 1;
    } 
    else 
    {
        Bignum rem = n % 8;
        if (rem == 1 || rem == 7) 
        {
            s = 1;
        } 
        else 
        {
            s = -1;
        }
    }

    if (n % 4 == 3 && b % 4 == 3) 
    {
        s = -s;
    }

    if (b == 1) 
    {
        return s;
    }

    return s * jacobi(n % b, b);
}

bool Bignum::solovay_strassen(size_t t) 
{
    if (*this == 3 || *this == 2) 
    {
        return true;
    }

    if (coef[0] % 2 == 0) 
    {
        return false;
    }

    Bignum n = *this - 1;
    while (t) 
    {
        Bignum a = Bignum(len, 3);
        if (a < Bignum(BASE(2)) || a > *this - 2) 
        {
            continue;
        }

        Bignum r = a.pow_m(n / 2, *this);
        if (r != 1 && r != n) 
        {
            return false;
        }

        int8_t s = jacobi(a, *this);
        if (s == 1) 
        {
            if (r != 1) 
            {
                return false;
            }
        } 
        else 
        {
            if (r != n) 
            {
                return false;
            }
        }

        t--;
    }

    return true;
}

Bignum gen_prime(size_t len) 
{
    Bignum prime;
    do 
    {
        prime = Bignum(len, 3);
    } while (!prime.solovay_strassen(100));
    
    return prime;
}

Bignum Bignum::gen_strong(size_t len) 
{
    Bignum s = gen_prime(len);
    Bignum t = gen_prime(len);

    Bignum i(1, 3);
    Bignum r = Bignum(BASE(2)) * i * t + 1;
    while (!r.solovay_strassen(100)) 
    {
        i += 1;
        r = Bignum(BASE(2)) * i * t + 1;
    }

    Bignum p0 = Bignum(BASE(2)) * s * s.pow_m(r - 2, r) - 1;
    Bignum j(1, 3);
    Bignum p = Bignum(BASE(2)) * j * r * s + p0;
    while (!p.solovay_strassen(100)) 
    {
        j += 1;
        p = Bignum(BASE(2)) * j * r * s + p0;
    }

    return p;
}

Bignum Bignum::sqrt_num(Bignum num)
{
    if (num == 1)
    {
        return num;
    }
    
    Bignum x = num / 2;
    Bignum x0;
    do
    {
        x0 = x;
        x = (((num / x) + x) / 2);
    } 
    while (x < x0);

    return x0;
} 

void Bignum::trial_div_method(vector <Bignum> &p, Bignum n)
{
    Bignum one = BASE(1);
    Bignum two = BASE(2);
    Bignum three = BASE(3);
    Bignum six = BASE(6);

    p.push_back(one);

    if (n == one)
    {
        return;
    }
    if (n == two)
    {
        p.push_back(two);
        return;
    }
    if (n == three)
    {
        p.push_back(three);
        return;
    }

    while (n % two == 0)
    {
        p.push_back(two);
        n = n / two;
    }

    if (n.solovay_strassen(100))
    {
        p.push_back(n);
        return;
    }

    Bignum d1 = BASE(3); 
    Bignum d2 = BASE(5);
    Bignum d3 = BASE(7);

    while (n != one)
    {
        Bignum r;
        r = n % d1;
        if (r == 0)
        {
            p.push_back(d1);
            n = n / d1;
            continue;
        }
        Bignum q;
        q = n / d1;
        if (q > d1)
        {
            d1 = d2;
            d2 = d3;
            d3 = d1 + six;
        }
        else
        {
            p.push_back(n);
            return;
        }
    }

    return;
}

Bignum Bignum::alway(Bignum n)
{
    Bignum two = BASE(2);
    while (n % two == 0)
    {
        n = n / two;
    }
    
    Bignum d = sqrt_num(n);
    d = (sqrt_num(d) * BASE(2)) + 9;
    
    Bignum r1;
    r1 = n % d;
    Bignum r2;
    Bignum m = (d - 2);
    r2 = n % m;
    Bignum q = ((n / m) - (n / d)) * BASE(4);
    Bignum s = sqrt_num(n);

    while (r1 != 0)
    {
        d = d + 2;
        if (d > s)
        {
            cout << "no divider" << endl;
            return 0;
        }

        Bignum r;
        if ((r1 * 2) >= r2)
        {
            r = ((r1 * 2) - r2 + q);
            r2 = r1;
            r1 = r;
        }

        if ((r1 * 2) < r2)
        {
            r = ((q + d) - (r2 - (r1 * 2)));
            r2 = r1;
            r1 = r;
            q = q + 4;
        }
        // r = ((r1 * 2) - r2 + q);
        // r2 = r1;
        // r1 = r;

        if (r1 < 0)
        {
            r1 = r1 + d;
            q = q + 4;
        }

        while (r1 >= d)
        {
            r1 = r1 - d;
            q = q - 4;
        }
    }

    return d;
}

void Bignum::method_ferma(pair<Bignum, Bignum> &res, Bignum num)
{
    Bignum two = BASE(2);
    while (num % two == 0)
    {
        num = num / two;
    }

    Bignum x = sqrt_num(num);
    if (x.fast_sq() == num)
    {
        res.first = res.second = x;
        res.first.print();
        cout << " * ";
        res.second.print();
        cout << endl;
        return;
    }
    
    Bignum y;
    Bignum z;
    do
    {
        x = x + 1;
        if (x == ((num + 1) / 2))
        {
            cout << "num - prime" << endl;
            return;
        }
        z = x.fast_sq() - num;
        y = sqrt_num(z);
    } while (y.fast_sq() != z);
    
    res.first = x + y;
    res.first.print();
    cout << " * ";
    res.second = x - y;
    res.second.print();
    cout << endl;

    return;
}

Bignum Bignum::method_p_pollard(Bignum num)
{
    Bignum two = BASE(2);
    while (num % two == 0)
    {
        num = num / two;
    }

    Bignum a = BASE(2);
    Bignum b = BASE(2);
    Bignum d;
    Bignum m;
    Bignum n;
    do
    {
        a = (a.fast_sq() + 1) % num;
        b = (b.fast_sq() + 1) % num;
        b = (b.fast_sq() + 1) % num;

        if (a == b)
        {
            cout << "no result" << endl;
            return 0;
        }

        if (a < b)
        {
            m = (b - a);
        }
        else
        {
            m = (a - b);
        }

        n = num;
        while (m != n)
        {
            if (m > n)
            {
                m = m - n;
            }
            else
            {
                n = n - m;
            }

            d = m;
        }
    } while (d == 1);

    return d;
}

Bignum Bignum::Gelfond(Bignum g, Bignum n, Bignum a)
{
    Bignum p = n + 1;
    Bignum h = sqrt_num(n) + 1;
    Bignum b = g.pow_m(h, p);

    map <Bignum, Bignum> table;
    for (Bignum u = 1; u <= h; u = u + 1)
    {
        Bignum tmp = b.pow_m(u, p);
        table[tmp] = u;
    }

    map <Bignum, Bignum> :: iterator it = table.begin();
    // for (int i = 0; it != table.end(); it++, i++) 
    // {  
    //     cout << i << ") Key ";
    //     Bignum temp1 = it->first;
    //     temp1.print();
    //     cout << ", value ";
    //     Bignum temp2 = it->second;
    //     temp2.print();
    //     cout << endl;
    // }
    // it = table.begin();
    map <Bignum, Bignum> :: iterator iter = table.end();
    long long l = 1;
    long long r = table.size();
    long long m = 0;
    Bignum temp = (a * g) % p;
    Bignum v = 1;

    while (v < h)
    {
        while (l < r - 1)
        {
            m = (l + r) / 2;
            it = table.begin();
            advance(it, m);
            // Bignum mid = it->first;
            // mid.print();
            // cout << " = mid" << endl;
            if (it->first < temp)
            {
                l = m;
            }
            else
            {
                r = m;
            }
        }
        it = table.begin();
        advance(it, r);
        if (it->first == temp)
        {
            break;
        }
        l = 0;
        r = table.size();
        temp = (g * temp) % p;
        v = v + 1;
    }
    
    // Bignum coincidence = it->first;
    // coincidence.print();
    // cout << " - coincidence" << endl;
    Bignum u = it->second;
    Bignum x = ((h * u) - v) % n;

    return x;
}
