#include <iostream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>

enum Rounding { ROUND_ABS_DOWN, ROUND_ABS_UP, ROUND_FLOOR, ROUND_CEIL };

template <class Number>
class Interval {

private:
    static uint32_t OUTPUT_PRECISION;
    static uint32_t TRUNC_PRECISION;
    static uint32_t OVERPRECISION_QUOTIENT;

    Number left_;
    Number right_;

    // Технические функции

    static Number NumberMin(std::vector<Number> numbers) {
        Number cmin = numbers[0];
        for (uint32_t i = 1; i < numbers.size(); ++i) {
            if (numbers[i] < cmin) {
                cmin = numbers[i];
            }
        }
        return cmin;
    }

    static Number NumberMax(std::vector<Number> numbers) {
        Number cmax = numbers[0];
        for (uint32_t i = 1; i < numbers.size(); ++i) {
            if (numbers[i] > cmax) {
                cmax = numbers[i];
            }
        }
        return cmax;
    }

    uint32_t NullType() {
        if (left_ != left_ || right_ != right_) {
            return 4;
        }
        int sigleft = left_.Sig();
        int sigright = right_.Sig();
        if (sigleft != sigright) {
            return 3;
        }
        if (sigleft == 1 && left_.IsZero()) {
            return 2;
        }
        if (sigright == -1 && right_.IsZero()) {
            return 1;
        }
        return 0;
    }

    Interval GetNan() {
        return Interval(left_.GetNan(), right_.GetNan());
    }

    Interval GetInf() {
        return Interval(-left_.GetInf(), right_.GetInf());
    }

    static void SetRounding(Rounding type) {
        Number::SetRounding(type);
    }
public:
    // Статические функции

    static void SetPrecision(uint32_t precision) {
        OUTPUT_PRECISION = precision;
        TRUNC_PRECISION = precision;
        Number::SetOutputPrecision(OUTPUT_PRECISION);
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Number::SetupPrecision();
        Number::SetCalcPrecision(TRUNC_PRECISION);
    }

    static void SetOutputPrecision(uint32_t precision) {
        OUTPUT_PRECISION = precision;
        Number::SetOutputPrecision(OUTPUT_PRECISION);
    }

    static void SetTruncPrecision(uint32_t precision) {
        TRUNC_PRECISION = precision;
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Number::SetupPrecision();
        Number::SetCalcPrecision(TRUNC_PRECISION);
    }

    static void SetOverPrecisionQuotient(uint32_t quotient) {
        OVERPRECISION_QUOTIENT = quotient;
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Number::SetupPrecision();
        Number::SetCalcPrecision(TRUNC_PRECISION);
    }

    // Конструкторы

    Interval() {}
    Interval(const Number& left, const Number& right) : left_(left), right_(right) {}
    Interval(const std::string& left, const std::string& right) : left_(left), right_(right) {}
    Interval(const std::vector<Number>& values) : left_(values[0]), right_(values[1]) {}
    Interval(const std::vector<std::string>& values) : left_(values[0]), right_(values[1]) {}
    Interval(const Interval<Number>& values) : left_(values.left_), right_(values.right_) {}
    Interval(const std::pair<Number, Number>& values) : left_(values.first), right_(values.second) {}
    Interval(const std::pair<Number, std::string>& values) : left_(values.first), right_(values.second) {}
    Interval(const std::pair<std::string, Number>& values) : left_(values.first), right_(values.second) {}
    Interval(const std::pair<std::string, std::string>& values) : left_(values.first), right_(values.second) {}

    // Технические функции

    void Normalize() {
        if (right_ < left_) {
            Number temp = left_;
            left_ = right_;
            right_ = temp;
        }
    }

    Number& first() {
        return left_;
    }
    Number& second() {
        return right_;
    }
    Number& left() {
        return left_;
    }
    Number& right() {
        return right_;
    }


    // Операторы

    Interval<Number>& operator=(std::vector<Number> values) {
        left_ = values[0];
        right_ = values[1];
        return this;
    }
    Interval<Number>& operator=(Interval<Number> values) {
        left_ = values.left_;
        right_ = values.right_;
        return *this;
    }
    Interval<Number>& operator=(std::pair<Number, Number> values) {
        left_ = values.first;
        right_ = values.second;
        return *this;
    }

    Number& operator[](uint32_t i) {
        return i ? right_ : left_;
    }

    bool operator==(Interval<Number> other) {
        return (left_ == other.left_ && right_ == other.right_ && NullType() == other.NullType());
    }
    bool operator!=(Interval<Number> other) {
        return (left_ != other.left_ || right_ != other.right_ || NullType() != other.NullType());
    }

    Interval<Number> operator+() {
        return Interval(left_, right_);
    }
    Interval<Number> operator+(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = left_ + other.left_;
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = right_ + other.right_;
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    Interval<Number>& operator+=(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        left_ += other.left_;
        Interval<Number>::SetRounding(ROUND_CEIL);
        right_ += other.right_;
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return *this;
    }

    Interval<Number> operator-() {
        return Interval(-right_, -left_);
    }
    Interval<Number> operator-(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = left_ - other.right_;
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = right_ - other.left_;
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    Interval<Number>& operator-=(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        left_ -= other.right_;
        Interval<Number>::SetRounding(ROUND_CEIL);
        right_ -= other.left_;
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return *this;
    }

    Interval<Number> operator*(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = NumberMin({ left_ * other.left_, left_ * other.right_,
                                     right_ * other.left_, right_ * other.right_ });
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = NumberMax({ left_ * other.left_, left_ * other.right_,
                                      right_ * other.left_, right_ * other.right_ });
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    Interval<Number>& operator*=(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        left_ = NumberMin({ left_ * other.left_, left_ * other.right_,
                            right_ * other.left_, right_ * other.right_ });
        Interval<Number>::SetRounding(ROUND_CEIL);
        right_ = NumberMax({ left_ * other.left_, left_ * other.right_,
                             right_ * other.left_, right_ * other.right_ });
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return *this;
    }

    Interval<Number> operator/(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        uint32_t stype = NullType();
        uint32_t otype = other.NullType();
        if (stype == 4 || otype == 4) {
            Number::SetCalcPrecision(TRUNC_PRECISION);
            return GetNan();
        }
        if ((stype == 3 && otype == 0) || (otype < 3 && stype < 3)) {
            Interval<Number>::SetRounding(ROUND_FLOOR);
            Number newleft = NumberMin({ left_ / other.left_, left_ / other.right_,
                                         right_ / other.left_, right_ / other.right_ });
            Interval<Number>::SetRounding(ROUND_CEIL);
            Number newright = NumberMax({ left_ / other.left_, left_ / other.right_,
                                          right_ / other.left_, right_ / other.right_ });
            Number::SetCalcPrecision(TRUNC_PRECISION);
            return Interval(newleft, newright);
        }
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return GetInf();
    }
    Interval<Number>& operator/=(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        int stype = NullType();
        int otype = other.NullType();
        if (stype == 4 || otype == 4) {
            left_ = left_.GetNan();
            right_ = right_.GetNan();
            Number::SetCalcPrecision(TRUNC_PRECISION);
            return *this;
        }
        if ((stype == 3 && otype == 0) || (otype < 3 && stype < 3)) {
            Interval<Number>::SetRounding(ROUND_FLOOR);
            left_ = NumberMin({ left_ / other.left_, left_ / other.right_,
                                right_ / other.left_, right_ / other.right_ });
            Interval<Number>::SetRounding(ROUND_CEIL);
            right_ = NumberMax({ left_ / other.left_, left_ / other.right_,
                                 right_ / other.left_, right_ / other.right_ });
            Number::SetCalcPrecision(TRUNC_PRECISION);
            return *this;
        }
        left_ = -left_.GetInf();
        right_ = right_.GetInf();
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return *this;
    }


    // Трансцендентные функции

    Interval<Number> Pow(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = NumberMin({ left_.Pow(other.left_), left_.Pow(other.right_),
                                     right_.Pow(other.left_), right_.Pow(other.right_) });
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = NumberMax({ left_.Pow(other.left_), left_.Pow(other.right_),
                                      right_.Pow(other.left_), right_.Pow(other.right_) });
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    
    static Interval<Number> Cos(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number yrd1 = Number::Cos(other.left_);
        Number yrd2 = Number::Cos(other.right_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number yru1 = Number::Cos(other.left_);
        Number yru2 = Number::Cos(other.right_);
        Number pi = Number::GetPi();
        Number twopi = pi << 1u;

        Number tn1 = other.left_ / twopi;
        Number::SetCalcPrecision(0);
        +tn1;
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Number tn2 = other.right_ / twopi;
        Number::SetCalcPrecision(0);
        +tn2;
        Number b;
        Interval<Number>::SetRounding(ROUND_CEIL);
        if (tn1 <= tn2) {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            b = Number(1);
        } else {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            b = NumberMax({ yru1, yru2 });
        }
        tn1 = (other.left_ - pi) / twopi;
        Number::SetCalcPrecision(0);
        +tn1;
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        tn2 = (other.right_ - pi) / twopi;
        Number::SetCalcPrecision(0);
        +tn2;
        Number a;
        Interval<Number>::SetRounding(ROUND_FLOOR);
        if (tn1 <= tn2) {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            a = Number(-1);
        } else {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            a = NumberMin({ yrd1, yrd2 });
        }
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(a, b);
    }
    static Interval<Number> Sin(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number yrd1 = Number::Sin(other.left_);
        Number yrd2 = Number::Sin(other.right_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number yru1 = Number::Sin(other.left_);
        Number yru2 = Number::Sin(other.right_);
        Number pi = Number::GetPi();
        Number halfpi = pi >> 1u;
        Number twopi = pi << 1u;

        Number tn1 = (other.left_ - halfpi) / twopi;
        Number::SetCalcPrecision(0);
        +tn1;
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Number tn2 = (other.right_ - halfpi) / twopi;
        Number::SetCalcPrecision(0);
        +tn2;
        Number b;
        Interval<Number>::SetRounding(ROUND_CEIL);
        if (tn1 <= tn2) {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            b = Number(1);
        } else {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            b = NumberMax({ yru1, yru2 });
        }
        tn1 = (other.left_ + halfpi) / twopi;
        Number::SetCalcPrecision(0);
        +tn1;
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        tn2 = (other.right_ + halfpi) / twopi;
        Number::SetCalcPrecision(0);
        +tn2;
        Number a;
        Interval<Number>::SetRounding(ROUND_FLOOR);
        if (tn1 <= tn2) {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            a = Number(-1);
        } else {
            Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
            a = NumberMin({ yrd1, yrd2 });
        }
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(a, b);
    }
    static Interval<Number> Tan(Interval<Number> other) {
        return Sin(other) / Cos(other);
    }    
    static Interval<Number> Cot(Interval<Number> other) {
        return Cos(other) / Sin(other);
    }

    static Interval<Number> Exp(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = Number::Exp(other.left_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = Number::Exp(other.right_);
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    static Interval<Number> Log(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = Number::Log(other.left_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = Number::Log(other.right_);
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }

    static Interval<Number> ATan(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = Number::ATan(other.left_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = Number::ATan(other.right_);
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    static Interval<Number> ASin(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = Number::ASin(other.left_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = Number::ASin(other.right_);
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    static Interval<Number> ACos(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = Number::ACos(other.right_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = Number::ACos(other.left_);
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }
    static Interval<Number> ACot(Interval<Number> other) {
        Number::SetCalcPrecision(TRUNC_PRECISION * OVERPRECISION_QUOTIENT);
        Interval<Number>::SetRounding(ROUND_FLOOR);
        Number newleft = Number::ACot(other.left_);
        Interval<Number>::SetRounding(ROUND_CEIL);
        Number newright = Number::ACot(other.right_);
        Number::SetCalcPrecision(TRUNC_PRECISION);
        return Interval(newleft, newright);
    }

    // Другое

    template <class SNumber>
    friend std::ostream& operator<<(std::ostream&, Interval<SNumber>);
};

template <class SNumber>
std::ostream& operator<<(std::ostream& os, Interval<SNumber> a) {
    os << "[";
    Interval<SNumber>::SetRounding(ROUND_FLOOR);
    os << a[0] << ", ";
    Interval<SNumber>::SetRounding(ROUND_CEIL);
    os << a[1] << "]";
    return os;
}

class StringLongInteger {
public:
    static const int32_t BASE = 1000000000;
    static const int32_t MAXITEM = 999999999;
    static const int32_t TWOBASE = 2000000000;
    static const int32_t HALFBASE = 500000000;
    static const int32_t QUARTERBASE = 250000000;
    bool sign_ = false;
    std::vector<uint32_t> data_;
    void RLZ() {
        while (data_.size() > 1 && data_.back() == 0) {
            data_.pop_back();
        }
        if (data_.size() == 1 && data_[0] == 0) sign_ = false;
        if (data_.size() == 0) {
            data_.push_back(0);
        }
    }
    StringLongInteger() {}

    StringLongInteger(uint32_t somenumber) : sign_(false) {
        data_.push_back(somenumber % StringLongInteger::BASE);
        if (somenumber >= StringLongInteger::BASE) {
            data_.push_back(somenumber / StringLongInteger::BASE);
        }
        RLZ();
    }
    StringLongInteger(int32_t somenumber) : sign_(false) {
        if (somenumber < 0) {
            sign_ = true;
            data_.push_back(-somenumber % StringLongInteger::BASE);
            if (-somenumber >= StringLongInteger::BASE) {
                data_.push_back(-somenumber / StringLongInteger::BASE);
            }
        } else {
            data_.push_back(somenumber % StringLongInteger::BASE);
            if (somenumber >= StringLongInteger::BASE) {
                data_.push_back(somenumber / StringLongInteger::BASE);
            }
        }
        RLZ();
    }

    StringLongInteger(std::string str) : sign_(false) {
        uint32_t strlen = str.length();
        if (strlen > 0) {
            sign_ = str[0] == '-';
            if (sign_) str = str.substr(1);
        }
        for (int64_t i = str.length(); i > 0; i -= 9) {
            if (i < 9)
                data_.push_back(atoi(str.substr(0, i).c_str()));
            else
                data_.push_back(atoi(str.substr(i - 9, 9).c_str()));
        }
        RLZ();
    }

    StringLongInteger SUB(uint32_t n) {
        StringLongInteger result;
        result.sign_ = sign_;
        if (n == 0) {
            result.data_.resize(1);
            return result;
        }
        result.data_.resize(n);
        uint32_t sz = data_.size();
        for (uint32_t i = 0; i < n; ++i) {
            if (i >= sz) {
                result.data_[i] = 0;
            } else {
                result.data_[i] = data_[i];
            }
        }
        return result;
    }
    uint32_t Get10N(uint32_t n) {
        switch (n) {
        case 0:
            return 1;
        case 1:
            return 10;
        case 2:
            return 100;
        case 3:
            return 1000;
        case 4:
            return 10000;
        case 5:
            return 100000;
        case 6:
            return 1000000;
        case 7:
            return 10000000;
        case 8:
            return 100000000;
        default:
            return 1000000000;
        }
    }
    bool SLHard(uint32_t n) {
        bool res = SL(n / 9);
        n = n % 9;
        uint32_t sz = data_.size();
        if (n && (data_[0] != 0 || sz != 1)) {
            uint32_t cutpow = Get10N(n);
            uint32_t shiftpow = Get10N(9 - n);
            uint32_t shift = data_[0] % cutpow;
            data_[0] /= cutpow;
            res |= shift;
            for (uint32_t i = 1; i < sz; ++i) {
                uint32_t shift = data_[i] % cutpow;
                data_[i - 1] += shiftpow * shift;
                data_[i] /= cutpow;
            }
            if (data_[sz - 1] == 0 && sz > 1) {
                data_.pop_back();
            }
        }
        return res;
    }
    bool SL(uint32_t n = 1) {
        bool res = false;
        if (data_.size() <= n) {
            if (data_.size() == 0) {
                data_ = { 0 };
                sign_ = false;
                return false;
            }
            if (data_.size() > 1 || data_[0] != 0) {
                res = true;
            }
            sign_ = false;
            data_ = { 0 };
            return res;
        }
        uint32_t sz = data_.size();
        for (uint32_t i = 0; i < n; ++i) {
            if (data_[i]) {
                res = true;
            }
        }
        for (uint32_t i = n; i < sz; ++i) {
            data_[i - n] = data_[i];
        }
        for (uint32_t i = 0; i < n; ++i)
            data_.pop_back();
        return res;
    }
    void SRHard(uint32_t n) {
        if (n % 9) {
            SR(n / 9 + 1);
            SLHard(9 - (n % 9));
        } else {
            SR(n / 9);
        }
    }
    void SR(uint32_t n = 1) {
        if (data_.size() == 0) {
            data_.push_back(0);
            return;
        }
        uint32_t sz = data_.size();
        data_.resize(sz + n);
        for (uint32_t i = data_.size(); i > n; --i)
            data_[i - 1] = data_[i - 1 - n];
        for (uint32_t i = 0; i < n; ++i)
            data_[i] = 0;
    }

    void Mul2() {
        RLZ();
        bool carry = 0;
        for (uint32_t i = 0; i < data_.size(); ++i) {
            data_[i] <<= 1;
            if (carry) {
                data_[i] += 1;
            }
            if (data_[i] > BASE) {
                carry = true;
                data_[i] -= BASE;
            } else {
                carry = false;
            }
        }
        if (carry) {
            data_.push_back(1);
        }
    }
    void Mul4() {
        RLZ();
        uint32_t carry = 0;
        for (uint32_t i = 0; i < data_.size(); ++i) {
            data_[i] <<= 2;
            if (carry) {
                data_[i] += carry;
            }
            carry = 0;
            if (data_[i] > TWOBASE) {
                carry += 2;
                data_[i] -= TWOBASE;
            }
            if (data_[i] > BASE) {
                carry += 1;
                data_[i] -= BASE;
            }
        }
        if (carry) {
            data_.push_back(carry);
        }
    }
    bool Div2() {
        bool carry = false;
        for (uint32_t i = data_.size(); i > 0; --i) {
            bool nextcarry = data_[i - 1] & 1;
            data_[i - 1] >>= 1;
            if (carry) {
                data_[i - 1] += HALFBASE;
            }
            carry = nextcarry;
        }
        RLZ();
        return carry;
    }
    bool Div4() {
        uint32_t carry = 0;
        for (uint32_t i = data_.size(); i > 0; --i) {
            uint32_t nextcarry = data_[i - 1] & 3;
            data_[i - 1] >>= 2;
            if (carry & 1) {
                data_[i - 1] += QUARTERBASE;
            }
            if (carry & 2) {
                data_[i - 1] += HALFBASE;
            }
            carry = nextcarry;
        }
        RLZ();
        return carry;
    }
    void Inc() {
        if (data_[0] != MAXITEM) {
            data_[0] += 1;
        } else {
            *this += 1;
        }
    }
    void Dec() {
        if (data_[0] != 0) {
            data_[0] -= 1;
        }
        else {
            *this -= 1;
        }
    }
    void IncOrDec(int32_t a) {
        if (a >= 0) {
            Inc();
        } else {
            Dec();
        }
    }

    const StringLongInteger operator -() const {
        StringLongInteger copy(*this);
        copy.sign_ = !copy.sign_;
        return copy;
    }

    bool operator ==(const StringLongInteger& o) const {
        if (sign_ != o.sign_) return o.data_[0] == 0 && data_[0] == 0 && o.data_.size() == 1 && data_.size() == 1;
        if (data_.empty()) {
            if (o.data_.empty() || (o.data_.size() == 1 && o.data_[0] == 0)) return true;
            else return false;
        }

        if (o.data_.empty()) {
            if (data_.size() == 1 && data_[0] == 0) return true;
            else return false;
        }
        if (data_.size() != o.data_.size()) return false;
        for (uint32_t i = 0; i < data_.size(); ++i) if (data_[i] != o.data_[i]) return false;

        return true;
    }
    bool operator !=(const StringLongInteger& o) const {
        return !((*this) == o);
    }
    bool operator <(const StringLongInteger& o) const {
        if ((*this) == o) 
            return false;
        if (sign_) {
            if (o.sign_) return ((-o) < (-(*this)));
            else return true;
        }
        else if (o.sign_) return false;
        else {
            if (data_.size() != o.data_.size()) {
                return data_.size() < o.data_.size();
            }
            else {
                for (uint32_t i = data_.size(); i > 0; --i) {
                    if (data_[i - 1] != o.data_[i - 1]) return data_[i - 1] < o.data_[i - 1];
                }

                return false;
            }
        }
    }
    bool operator >(const StringLongInteger& o) const {
        return o < (*this);
    }
    bool operator <=(const StringLongInteger& o) const {
        return (*this) < o || o == (*this);
    }
    bool operator >=(const StringLongInteger& o) const {
        return o < (*this) || o == (*this);
    }
    StringLongInteger& operator+=(StringLongInteger o) {
        RLZ();
        o.RLZ();
        if (sign_ != o.sign_) {
            if (sign_) {
                StringLongInteger temp = o;
                temp -= (-*this);
                (*this) = temp;
                return *this;
            } else {
                *this -= (-o);
                return *this;
            }
        }
        int32_t carry = 0;

        uint32_t max_size = std::max(o.data_.size(), data_.size());
        for (uint32_t i = 0; i < max_size || carry != 0; ++i) {
            if (i == data_.size()) data_.push_back(0);
            data_[i] += carry + (i < o.data_.size() ? o.data_[i] : 0);
            carry = data_[i] >= StringLongInteger::BASE;
            if (carry != 0) data_[i] -= StringLongInteger::BASE;
        }

        RLZ();
        return *this;
    }
    StringLongInteger operator+(const StringLongInteger& o) const {
        StringLongInteger temp = (*this);
        return temp += o;
    }
    StringLongInteger& operator-=(StringLongInteger o) {
        RLZ();
        o.RLZ();
        if (o.sign_) { 
            *this += (-o);
            return *this;
        }
        else if ((*this).sign_) {
            StringLongInteger temp = -(*this);
            temp += o;
            (*this) = -temp;
            return *this;
        }
        else if ((*this) < o) {
            StringLongInteger temp = o;
            temp -= (*this);
            (*this) = -temp;
            return *this;
        }
        int32_t carry = 0;

        uint32_t max_size = o.data_.size();
        for (uint32_t i = 0; i < max_size || carry != 0; ++i) {
            uint32_t sub = carry + (i < o.data_.size() ? o.data_[i] : 0);

            carry = data_[i] < sub;
            if (carry != 0) data_[i] += StringLongInteger::BASE;
            data_[i] -= sub;
        }

        RLZ();
        return *this;
    }
    StringLongInteger operator-(const StringLongInteger& o) const {
        StringLongInteger temp = (*this);
        return temp -= o;
    }
    StringLongInteger& operator*=(StringLongInteger o) {
        RLZ();
        o.RLZ();
        if (o == 2) {
            Mul2();
            return *this;
        }
        if (o == 4) {
            Mul4();
            return *this;
        }
        
        StringLongInteger result;
        result.data_.resize(data_.size() + o.data_.size());

        uint32_t max_size = data_.size();
        for (uint32_t i = 0; i < max_size; ++i) {
            int32_t carry = 0;
            uint32_t o_max_size = o.data_.size();
            for (uint32_t j = 0; j < o_max_size || carry != 0; ++j) {
                int64_t cur = result.data_[i + j] +
                    data_[i] * 1LL * (j < o.data_.size() ? o.data_[j] : 0) + carry;
                result.data_[i + j] = static_cast<int32_t>(cur % StringLongInteger::BASE);
                carry = static_cast<int32_t>(cur / StringLongInteger::BASE);
            }
        }
        result.sign_ = sign_ != o.sign_;
        *this = result;
        RLZ();
        return *this;
    }
    StringLongInteger operator*(const StringLongInteger& o) const {
        StringLongInteger temp = (*this);
        return temp *= o;
    }

    std::pair<StringLongInteger, StringLongInteger> DDiv(const StringLongInteger& o) const {
        if (o == 0) throw std::runtime_error("Zero division");
        StringLongInteger b = o;
        b.sign_ = false;
        StringLongInteger result, current;
        result.data_.resize(data_.size());
        for (int64_t i = static_cast<int64_t>(data_.size()) - 1; i >= 0; --i) {
            current.SR();
            current.data_[0] = data_[i];
            current.RLZ();
            int32_t x = 0, l = 0, r = StringLongInteger::BASE;
            while (l <= r) {
                int32_t m = (l + r) / 2;
                StringLongInteger t = b;
                t *= m;
                if (t <= current) {
                    x = m;
                    l = m + 1;
                }
                else r = m - 1;
            }

            result.data_[i] = x;
            StringLongInteger t = b;
            t *= x;
            current -= t;
        }

        result.sign_ = sign_ != o.sign_;
        result.RLZ(); 
        StringLongInteger mod = (*this) - result * o;
        if (mod.sign_) mod += o;
        return { result, mod };
    }

    StringLongInteger& operator/=(StringLongInteger o) {
        RLZ();
        o.RLZ();
        if (o == 0) throw std::runtime_error("Zero division");
        if (o == 2) {
            Div2();
            return *this;
        }
        if (o == 4) {
            Div4();
            return *this;
        }
        StringLongInteger b = o;
        b.sign_ = false;
        StringLongInteger result, current;
        result.data_.resize(data_.size());
        for (int64_t i = static_cast<int64_t>(data_.size()) - 1; i >= 0; --i) {
            current.SR();
            current.data_[0] = data_[i];
            current.RLZ();
            int32_t x = 0, l = 0, r = StringLongInteger::BASE;
            while (l <= r) {
                int32_t m = (l + r) / 2;
                StringLongInteger t = b;
                t *= m;
                if (t <= current) {
                    x = m;
                    l = m + 1;
                }
                else r = m - 1;
            }

            result.data_[i] = x;
            StringLongInteger t = b;
            t *= x;
            current -= t;
        }

        result.sign_ = (sign_ ^ o.sign_);
        *this = result;
        RLZ();
        return *this;
    }
    StringLongInteger operator/(const StringLongInteger& o) const {
        StringLongInteger temp = (*this);
        return temp /= o;
    }

    bool IsOdd() const {
        if (data_.size() == 0) return false;
        return data_[0] & 1;
    }

    const StringLongInteger pow(StringLongInteger n) const {
        StringLongInteger a(*this), result(1);
        while (n != 0) {
            if (n.IsOdd()) result *= a;
            a *= a;
            n /= 2;
        }

        return result;
    }
};
const int StringLongInteger::BASE;

std::ostream& operator<<(std::ostream& os, const StringLongInteger& n) {
    if (n.data_.size() == 0) {
        os << 0;
    } else {
        for (uint32_t i = n.data_.size(); i > 0; --i) {
            if (i == n.data_.size()) {
                os << n.data_[i - 1];
            } else {
                for (uint32_t j = 1000000000; j >= 10; j /= 10) {
                    os << n.data_[i - 1] % j / (j / 10);
                }
            }
        }
    }
    return os;
}

class FixedWrapper {

public:
    static uint32_t CALCPRECISION;
    static uint32_t OUTPRECISION;
    static Rounding ROUNDING;

    uint32_t localfixed_;
    StringLongInteger base_;
    bool spec_;
    StringLongInteger ABS(StringLongInteger a) const {
        return a.sign_ ? -a : a;
    }
    FixedWrapper abs() {
        UpdateFixedPoint();
        return Sig() == -1 ? -(*this) : (*this);
    }
    void CheckPrecision() {
        if (CALCPRECISION % 9) {
            CALCPRECISION = CALCPRECISION + (9 - CALCPRECISION % 9);
        }
    }

    virtual void PARSE(std::string s) {
        CheckPrecision();
        if (s == "INF" || s == "inf" || s == "infty" || s == "+INF" || s == "+inf" || s == "+infty") {
            spec_ = true;
            localfixed_ = 0;
            base_ = StringLongInteger("1");
            return;
        }        
        if (s == "-INF" || s == "-inf" || s == "-infty") {
            spec_ = true;
            localfixed_ = 0;
            base_ = StringLongInteger("-1");
            return;
        }
        if (s == "-0") {
            spec_ = true;
            localfixed_ = 0;
            base_ = StringLongInteger("0");
            return;
        }
        if (s == "nan" || s == "NaN" || s == "NAN" || s == "nil" || s == "NIL" || s == "Nan" || s == "Nil") {
            spec_ = true;
            localfixed_ = 0;
            base_ = StringLongInteger("2");
            return;
        }

        bool sign = s[0] == '-';
        if (sign) s = s.substr(1);
        uint32_t fixedpos = s.size() + 1;
        uint32_t sz = fixedpos - 1;
        std::string news = "";
        news.reserve(sz + 8);
        for (uint32_t i = 0; i < sz; ++i) {
            if (s[i] == '.' || s[i] == ',') {
                if (fixedpos != sz + 1) {
                    spec_ = true;
                    localfixed_ = 0;
                    base_ = StringLongInteger("2");
                    return;
                }
                fixedpos = i;
            } else {
                if (s[i] < '0' || s[i] > '9') {
                    spec_ = true;
                    localfixed_ = 0;
                    base_ = StringLongInteger("2");
                    return;
                }
                news += s[i];
            }
        }
        int32_t prec = static_cast<int32_t>(sz) - fixedpos - 1;
        if (prec == -2) {
            prec = 0;
        }
        for (; prec % 9; ++prec) {
            news += '0';
        }
        StringLongInteger res(news);
        res.sign_ = sign;
        if (prec < CALCPRECISION) {
            uint32_t delta = (CALCPRECISION - prec);
            res.SRHard(delta);
        } else {
            uint32_t delta = (prec - CALCPRECISION);
            bool carry = res.SLHard(delta);
            switch (ROUNDING) {
            case ROUND_ABS_UP:
                if (carry) {
                    res.IncOrDec(Sig());
                }
                break;
            case ROUND_ABS_DOWN:
                break;
            case ROUND_CEIL:
                if (Sig() == 1 && carry) {
                    res.Inc();
                }
                break;
            case ROUND_FLOOR:
                if (Sig() == -1 && carry) {
                    res.Dec();
                }
                break;
            default:
                break;
            }
        }
        res.RLZ();
        base_ = res;
        return;
    }
    FixedWrapper() : base_(0), spec_(false) {}
    FixedWrapper(StringLongInteger base, bool spec = false, uint32_t localfixed = CALCPRECISION % 9 ? CALCPRECISION + 9 - CALCPRECISION % 9 : CALCPRECISION)
        : base_(base), spec_(spec), localfixed_(localfixed) {}
    FixedWrapper(std::string base, uint32_t localfixed = CALCPRECISION % 9 ? CALCPRECISION + 9 - CALCPRECISION % 9 : CALCPRECISION)
        : spec_(false), localfixed_(localfixed) {
        PARSE(base);
    }
    FixedWrapper(int32_t base, bool spec = false)
        : base_(base), spec_(spec), localfixed_(0) {}
    static void SetCalcPrecision(uint32_t prec) {
        CALCPRECISION = (prec / 9) * 9 + (prec % 9 ? 1 : 0);
    }
    static void SetOutputPrecision(uint32_t prec) {
        OUTPRECISION = prec;
    }
    static void SetRounding(Rounding type) {
        ROUNDING = type;
    }
    static void FlipRounding() {
        ROUNDING = Rounding(ROUNDING - ROUNDING % 2 + (ROUNDING + 1) % 2);
    }

    void UpdateFixedPoint() {
        if (spec_)
            return;
        CheckPrecision();
        if (localfixed_ != CALCPRECISION) {
            if (localfixed_ < CALCPRECISION) {
                uint32_t delta = (CALCPRECISION - localfixed_);
                base_.SRHard(delta);
                localfixed_ = CALCPRECISION;
            } else {
                uint32_t delta = (localfixed_ - CALCPRECISION);
                bool carry = base_.SLHard(delta);
                switch (ROUNDING) {
                case ROUND_ABS_UP:
                    if (carry) {
                        base_.IncOrDec(Sig());
                    }
                    break;
                case ROUND_ABS_DOWN:
                    break;
                case ROUND_CEIL:
                    if (Sig() == 1 && carry) {
                        base_.Inc();
                    }
                    break;
                case ROUND_FLOOR:
                    if (Sig() == -1 && carry) {
                        base_.Dec();
                    }
                    break;
                default:
                    break;
                }
                localfixed_ = CALCPRECISION;
            }
        }
        base_.RLZ();
    }
    FixedWrapper GetOutput() {
        UpdateFixedPoint();
        if (spec_) {
            return *this;
        }
        FixedWrapper res = *this;
        uint32_t outprec = OUTPRECISION;
        if (res.localfixed_ > outprec) {
            uint32_t delta = (res.localfixed_ - outprec);
            bool carry = res.base_.SLHard(delta);
            switch (ROUNDING) {
            case ROUND_ABS_UP:
                if (carry) {
                    res.base_.IncOrDec(Sig());
                }
                res.localfixed_ = outprec;
                break;
            case ROUND_ABS_DOWN:
                res.localfixed_ = outprec;
                break;
            case ROUND_CEIL:
                if (Sig() == 1 && carry) {
                    res.base_.Inc();
                }
                res.localfixed_ = outprec;
                break;
            case ROUND_FLOOR:
                if (Sig() == -1 && carry) {
                    res.base_.Dec();
                }
                res.localfixed_ = outprec;
                break;
            default:
                res.localfixed_ = outprec;
                break;
            }
        }
        return res;
    }
    void FastShift(uint32_t a, bool left = true) {
        if (spec_)
            return;
        if (left) {
            localfixed_ += a;
        } else {
            if (localfixed_ > a) {
                localfixed_ -= a;
            } else {
                a -= localfixed_;
                localfixed_ = 0;
                base_.SRHard(a);
            }
        }
    }
    void FastMul2() {
        UpdateFixedPoint();
        if (spec_) {
            return;
        }
        base_.Mul2();
    }
    void FastDiv2() {
        UpdateFixedPoint();
        if (spec_) {
            return;
        }
        bool carry = base_.Div2();
        switch (ROUNDING) {
        case ROUND_ABS_UP:
            if (carry) {
                base_.IncOrDec(Sig());
            }
            break;
        case ROUND_ABS_DOWN:
            break;
        case ROUND_CEIL:
            if (Sig() == 1 && carry) {
                base_.Inc();
            }
            break;
        case ROUND_FLOOR:
            if (Sig() == -1 && carry) {
                base_.Dec();
            }
            break;
        default:
            break;
        }
    }
    void FastDiv4() {
        UpdateFixedPoint();
        if (spec_) {
            return;
        }
        bool carry = base_.Div4();
        switch (ROUNDING) {
        case ROUND_ABS_UP:
            if (carry) {
                base_.IncOrDec(Sig());
            }
            break;
        case ROUND_ABS_DOWN:
            break;
        case ROUND_CEIL:
            if (Sig() == 1 && carry) {
                base_.Inc();
            }
            break;
        case ROUND_FLOOR:
            if (Sig() == -1 && carry) {
                base_.Dec();
            }
            break;
        default:
            break;
        }
    }

    FixedWrapper operator>>(uint32_t n) {
        FixedWrapper res = *this;
        for (uint32_t i = 0; i < n; ++i) {
            res.FastDiv2();
        }
        return res;
    }
    FixedWrapper operator<<(uint32_t n) {
        FixedWrapper res = *this;
        for (uint32_t i = 0; i < n; ++i) {
            res.FastMul2();
        }
        return res;
    }

    int Sig() const {
        if (!spec_) {
            return base_.sign_ ? -1 : 1;
        }
        if (base_ == 0) {
            return -1;
        }
        if (IsNan()) {
            return 0;
        }
        return base_.sign_ ? -1 : 1;
    }
    bool IsInf() const {
        return spec_ && ABS(base_) == 1;
    }
    bool IsPInf() const {
        return spec_ && base_ == 1;
    }
    bool IsNInf() const {
        return spec_ && base_ == -1;
    }
    bool IsZero() const {
        return base_ == 0;
    }
    bool IsPZero() const {
        return !spec_ && base_ == 0;
    }
    bool IsNZero() const {
        return spec_ && base_ == 0;
    }
    bool IsNan() const {
        return spec_ && base_ >= 2;
    }
    static FixedWrapper GetNan() {
        return FixedWrapper(2, true);
    }
    static FixedWrapper GetInf() {
        return FixedWrapper(1, true);
    }

    FixedWrapper& operator=(FixedWrapper other) {
        base_ = other.base_;
        spec_ = other.spec_;
        localfixed_ = other.localfixed_;
        UpdateFixedPoint();
        return *this;
    }

    bool operator==(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return false;
        }
        return base_ == o.base_;
    }
    bool operator!=(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return true;
        }
        return base_ != o.base_;
    }
    bool operator<(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return false;
        }
        if (IsPInf()) {
            return false;
        }
        if (o.IsPInf()) {
            return true;
        }
        if (IsNInf() && !o.IsNInf()) {
            return true;
        }
        if (o.IsNInf()) {
            return false;
        }
        return base_ < o.base_;
    }
    bool operator<=(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return false;
        }
        if (o.IsPInf()) {
            return true;
        }
        if (IsPInf()) {
            return false;
        }
        if (IsNInf() && !o.IsNInf()) {
            return true;
        }
        return base_ <= o.base_;
    }
    bool operator>(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return false;
        }
        if (o.IsPInf()) {
            return false;
        }
        if (IsPInf()) {
            return true;
        }
        if (o.IsNInf() && !IsNInf()) {
            return true;
        }
        return base_ > o.base_;
    }
    bool operator>=(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return false;
        }
        if (IsPInf()) {
            return true;
        }
        if (o.IsPInf()) {
            return false;
        }
        if (o.IsNInf() && !IsNInf()) {
            return true;
        }
        return base_ >= o.base_;
    }
    FixedWrapper operator+() {
        UpdateFixedPoint();
        return FixedWrapper(base_);
    }
    FixedWrapper operator-() {
        UpdateFixedPoint();
        if (base_ == 0) {
            return FixedWrapper(base_, !spec_);
        }
        return FixedWrapper(-base_, spec_);
    }
    FixedWrapper operator+(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return GetNan();
        }
        if (IsPInf()) {
            if (o.IsNInf()) return FixedWrapper(0);
            return FixedWrapper(1, true);
        }
        if (o.IsNInf()) return FixedWrapper(-1, true);
        if (IsNZero() && o.IsZero()) {
            return FixedWrapper(0, true);
        }
        return FixedWrapper(base_ + o.base_);
    }
    FixedWrapper& operator+=(FixedWrapper o) {
        FixedWrapper res = (*this) + o;
        base_ = res.base_;
        spec_ = res.spec_;
        return *this;
    }
    FixedWrapper operator-(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return GetNan();
        }
        if (IsPInf()) {
            if (o.IsPInf()) return FixedWrapper(0);
            return FixedWrapper(1, true);
        }
        if (o.IsNInf()) return FixedWrapper(1, true);
        if (IsNZero() && o.IsZero()) {
            return FixedWrapper(0, true);
        }
        return FixedWrapper(base_ - o.base_);
    }
    FixedWrapper& operator-=(FixedWrapper o) {
        FixedWrapper res = (*this) - o;
        base_ = res.base_;
        spec_ = res.spec_;
        return *this;
    }
    FixedWrapper operator*(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return GetNan();
        }
        if (IsInf() || o.IsInf()) {
            return FixedWrapper(o.Sig() * Sig(), o.base_ != 0 && base_ != 0);
        }
        if (IsZero() || o.IsZero()) {
            return FixedWrapper(0, Sig() * o.Sig() < 0);
        }
        FixedWrapper res = FixedWrapper(base_ * o.base_, false, CALCPRECISION * 2);
        res.UpdateFixedPoint();
        return res;
    }
    FixedWrapper& operator*=(FixedWrapper o) {
        FixedWrapper res = (*this) * o;
        base_ = res.base_;
        spec_ = res.spec_;
        return *this;
    }

    FixedWrapper div(FixedWrapper o) {
        FixedWrapper res = (*this) / o;
        if (res.spec_) {
            return res;
        }
        uint32_t oldprec = CALCPRECISION;
        bool carry = false;
        if (res.base_.SUB(CALCPRECISION / 9) != 0) {
            carry = true;
        }
        res.base_.SL(CALCPRECISION / 9);
        res.localfixed_ = 0;
        return res;
    }

    FixedWrapper operator/(FixedWrapper o) {
        UpdateFixedPoint();
        o.UpdateFixedPoint();
        if (IsNan() || o.IsNan()) {
            return GetNan();
        }
        if (IsInf()) {
            return FixedWrapper(o.Sig() * Sig(), !o.IsInf());
        }
        if (o.IsZero()) {
            return FixedWrapper(o.Sig() * Sig(), !IsZero());
        }
        if (o.IsInf() || IsZero()) {
            return FixedWrapper(0, o.Sig() * Sig() < 0);
        }
        StringLongInteger bs = base_;
        bs.SR(CALCPRECISION / 9);
        FixedWrapper res = FixedWrapper(bs / o.base_, false);
        return res;
    }
    FixedWrapper& operator/=(FixedWrapper o) {
        FixedWrapper res = (*this) / o;
        base_ = res.base_;
        spec_ = res.spec_;
        return *this;
    }

    FixedWrapper clamp(FixedWrapper low = FixedWrapper(-1, true), FixedWrapper high = FixedWrapper(1, true)) {
        UpdateFixedPoint();
        if (*this < low) {
            return low;
        }
        if (*this > high) {
            return high;
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream&, FixedWrapper);

    FixedWrapper MSPow(FixedWrapper);
    static FixedWrapper CalcPi(uint32_t);
    static FixedWrapper GetPi(uint32_t);
    static FixedWrapper MSExpUnchecked(FixedWrapper);
    static FixedWrapper CalcE(uint32_t);
    static FixedWrapper GetE(uint32_t);
    static FixedWrapper MSFasterBinExpUnchecked(FixedWrapper);
    static FixedWrapper MSLogUnchecked(FixedWrapper);
    static FixedWrapper GetLN1p5();
    static FixedWrapper MSFasterBinLogUnchecked(FixedWrapper);
    static FixedWrapper MSSinUnchecked(FixedWrapper);
    static FixedWrapper MSCosUnchecked(FixedWrapper);
    static FixedWrapper MSArctgUnchecked(FixedWrapper);
    static FixedWrapper MSSqrt(FixedWrapper);
    static void SetupCORDIC();
    static void SetupPrecision();
    static std::pair<FixedWrapper, FixedWrapper> CORDIC10(FixedWrapper);
    static std::pair<FixedWrapper, FixedWrapper> CORDIC2(FixedWrapper);
    static FixedWrapper A_CORDIC2(std::pair<FixedWrapper, FixedWrapper>);
    static FixedWrapper CORDIC2_Sin(FixedWrapper);
    static FixedWrapper CORDIC2_Cos(FixedWrapper);
    static FixedWrapper CORDIC2_Tan(FixedWrapper);
    static FixedWrapper CORDIC2_Cot(FixedWrapper);
    static FixedWrapper CORDIC2_ASin(FixedWrapper);
    static FixedWrapper CORDIC2_ACos(FixedWrapper);
    static FixedWrapper CORDIC2_ATan(FixedWrapper);
    static FixedWrapper CORDIC2_ACot(FixedWrapper);
    
    FixedWrapper Pow(FixedWrapper a) {
        return MSPow(a);
    }
    static FixedWrapper Exp(FixedWrapper a) {
        return MSFasterBinExpUnchecked(a);
    }
    static FixedWrapper Log(FixedWrapper a) {
        return MSFasterBinLogUnchecked(a);
    }
    static FixedWrapper Sin(FixedWrapper a) {
        return CORDIC2_Sin(a);
    }
    static FixedWrapper Cos(FixedWrapper a) {
        return CORDIC2_Cos(a);
    }
    static FixedWrapper Tan(FixedWrapper a) {
        return CORDIC2_Tan(a);
    }
    static FixedWrapper Cot(FixedWrapper a) {
        return CORDIC2_Cot(a);
    }
    static FixedWrapper ASin(FixedWrapper a) {
        return CORDIC2_ASin(a);
    }
    static FixedWrapper ACos(FixedWrapper a) {
        return CORDIC2_ACos(a);
    }
    static FixedWrapper ATan(FixedWrapper a) {
        return CORDIC2_ATan(a);
    }
    static FixedWrapper ACot(FixedWrapper a) {
        return CORDIC2_ACot(a);
    }
};
Rounding FixedWrapper::ROUNDING = ROUND_ABS_DOWN;
uint32_t FixedWrapper::OUTPRECISION = 18;
uint32_t FixedWrapper::CALCPRECISION = 12;

uint32_t Interval<FixedWrapper>::OUTPUT_PRECISION = 12;
uint32_t Interval<FixedWrapper>::TRUNC_PRECISION = 12;
uint32_t Interval<FixedWrapper>::OVERPRECISION_QUOTIENT = 2;

std::ostream& operator<<(std::ostream& os, FixedWrapper smt) {
    FixedWrapper n = smt.GetOutput();
    uint32_t outprec = FixedWrapper::OUTPRECISION;
    uint32_t outage = (9 - outprec % 9);
    outage = outage == 9 ? 0 : outage;
    n.base_ *= n.base_.Get10N(outage);
    n.localfixed_ += outage;
    uint32_t sz = n.base_.data_.size();
    uint32_t curprec = n.localfixed_ / 9;
    uint32_t prec = 0;
    bool sign = false;
    if (n.spec_) {
        if (n.IsNInf()) {
            os << "-INF";
            return os;
        }
        if (n.IsPInf()) {
            os << "+INF";
            return os;
        }
        if (n.IsNan()) {
            os << "nan";
            return os;
        }
    }
    if (n.Sig() == -1) {
        os << "-";
        sign = true;
    }
    n.spec_ = false;
    n.base_.sign_ = false;
    if (sz <= curprec) {
        if (outprec == 0) {
            os << "0";
            return os;
        }
        os << "0.";
        for (uint32_t i = curprec * 9; i > sz * 9; i--) {
            os << "0";
            prec += 1;
            if (prec == outprec) {
                break;
            }
        }
        for (uint32_t i = sz; i > 0; --i) {
            for (uint32_t j = 1000000000; j >= 10; j /= 10) {
                os << n.base_.data_[i - 1] % j / (j / 10);
                prec++;
                if (prec == outprec) {
                    return os;
                }
            }
        }
        for (; prec < outprec; ++prec) {
            os << "0";
        }
    } else {
        for (uint32_t i = sz; i > curprec; --i) {
            if (i == sz) {
                os << n.base_.data_[i - 1];
            } else {
                for (uint32_t j = 1000000000; j >= 10; j /= 10) {
                    os << n.base_.data_[i - 1] % j / (j / 10);
                }
            }
        }
        if (outprec > 0) {
            os << ".";
        }
        for (uint32_t i = curprec; i > 0; --i) {
            for (uint32_t j = 1000000000; j >= 10; j /= 10) {
                os << n.base_.data_[i - 1] % j / (j / 10);
                prec++;
                if (prec == outprec) {
                    break;
                }
            }
        }
        for (; prec < outprec; ++prec) {
            os << "0";
        }
    }
    return os;
}

class ConstantsHolder {
public:

    std::vector<std::pair<FixedWrapper, uint32_t>> BestPi = { {3, 0}, {4, 0}, {3, 0}, {4, 0} };
    std::vector<std::pair<FixedWrapper, uint32_t>> BestE = { {2, 0}, {3, 0}, {2, 0}, {3, 0} };
    FixedWrapper mone = FixedWrapper("-1");
    FixedWrapper zero = FixedWrapper("0");
    FixedWrapper onetenth = FixedWrapper("0.1");
    FixedWrapper fivetenth = FixedWrapper("0.5");
    FixedWrapper one = FixedWrapper("1");
    FixedWrapper oneandonetenth = FixedWrapper("1.1");
    FixedWrapper oneandfivetenth = FixedWrapper("1.5");
    FixedWrapper two = FixedWrapper("2");
    FixedWrapper three = FixedWrapper("3");
    FixedWrapper four = FixedWrapper("4");
    FixedWrapper five = FixedWrapper("5");
    FixedWrapper eight = FixedWrapper("8");
    FixedWrapper ten = FixedWrapper("10");
    FixedWrapper twentyfour = FixedWrapper("24");
    FixedWrapper thirtytwo = FixedWrapper("32");
    std::vector<std::pair<FixedWrapper, uint32_t>> BestLn1p5 =
        {{FixedWrapper("0.4"), 1}, {FixedWrapper("0.5"), 1},
         {FixedWrapper("0.4"), 1}, {FixedWrapper("0.5"), 1}};
    uint32_t CORDIC_MAX_PRECISION = 0;
    std::vector<std::vector<FixedWrapper>> CORDIC_DK10 = { {}, {}, {}, {} };
    std::vector<std::vector<FixedWrapper>> CORDIC_ATAN10 = { {}, {}, {}, {} };
    std::vector<std::vector<FixedWrapper>> CORDIC_P10i = { {}, {}, {}, {} };
    std::vector<std::vector<FixedWrapper>> CORDIC_K2 = { {}, {}, {}, {} };
    std::vector<std::vector<FixedWrapper>> CORDIC_DK2 = { {}, {}, {}, {} };
    std::vector<std::vector<FixedWrapper>> CORDIC_ATAN2 = { {}, {}, {}, {} };
    std::vector<std::vector<FixedWrapper>> CORDIC_P2i = { {}, {}, {}, {} };
    
};
ConstantsHolder CONSTANTS;

FixedWrapper FixedWrapper::MSPow(FixedWrapper x) {
    return MSFasterBinExpUnchecked(x * MSFasterBinLogUnchecked(*this));
}

FixedWrapper FixedWrapper::CalcPi(uint32_t prec = FixedWrapper::CALCPRECISION) {
    uint32_t originalprecision = FixedWrapper::CALCPRECISION;
    FixedWrapper::SetCalcPrecision(prec * 11 / 10 + 9);
    FixedWrapper lasts = CONSTANTS.zero;
    FixedWrapper t = CONSTANTS.three;
    FixedWrapper s = CONSTANTS.three;
    FixedWrapper n = CONSTANTS.one;
    FixedWrapper na = CONSTANTS.zero;
    FixedWrapper d = CONSTANTS.zero;
    FixedWrapper da = CONSTANTS.twentyfour;
    while (s != lasts) {
        lasts = s;
        n += na;
        na += CONSTANTS.eight;
        d += da;
        da += CONSTANTS.thirtytwo;
        t = t * n / d;
        s += t;
    } 

    if (prec > CONSTANTS.BestPi[FixedWrapper::ROUNDING].second) {
        CONSTANTS.BestPi[FixedWrapper::ROUNDING].first = s;
        CONSTANTS.BestPi[FixedWrapper::ROUNDING].second = prec;
    }
    FixedWrapper::SetCalcPrecision(originalprecision);
    +s;
    switch (FixedWrapper::ROUNDING) {
    case ROUND_ABS_DOWN:
        s.base_ -= s.Sig();
        break;
    case ROUND_ABS_UP:
        s.base_ += s.Sig();
        break;
    case ROUND_CEIL:
        s.base_ += 1;
        break;
    case ROUND_FLOOR:
        s.base_ -= 1;
        break;
    default:
        break;
    }

    return s;
}

FixedWrapper FixedWrapper::GetPi(uint32_t prec = FixedWrapper::CALCPRECISION) {
    if (prec > CONSTANTS.BestPi[FixedWrapper::ROUNDING].second) {
        CalcPi(prec);
    }
    return +CONSTANTS.BestPi[FixedWrapper::ROUNDING].first;
}

FixedWrapper FixedWrapper::MSExpUnchecked(FixedWrapper x) {
    if (x.IsNan()) {
        return x;
    }
    if (x.IsNInf()) {
        return CONSTANTS.zero;
    }
    uint32_t originalprecision = FixedWrapper::CALCPRECISION;
    FixedWrapper::SetCalcPrecision(originalprecision * 2);
    FixedWrapper i = 0;
    FixedWrapper lasts = 0;
    FixedWrapper s = 1;
    FixedWrapper fact = 1;
    FixedWrapper num = 1;
    FixedWrapper delta = 0;
    while (s != lasts) {
        lasts = s;
        i += CONSTANTS.one;
        fact *= i;
        num *= x;
        delta = num / fact;
        s += delta;
    }
    FixedWrapper::SetCalcPrecision(originalprecision);
    +s;
    switch (FixedWrapper::ROUNDING) {
    case ROUND_ABS_DOWN:
        s.base_ -= s.Sig();
        break;
    case ROUND_ABS_UP:
        s.base_ += s.Sig();
        break;
    case ROUND_CEIL:
        s.base_ += 1;
        break;
    case ROUND_FLOOR:
        s.base_ -= 1;
        break;
    default:
        break;
    }
    return s.clamp(CONSTANTS.zero);
}

FixedWrapper FixedWrapper::CalcE(uint32_t prec = FixedWrapper::CALCPRECISION) {
    uint32_t originalprecision = FixedWrapper::CALCPRECISION;
    FixedWrapper::SetCalcPrecision(prec * 11 / 10 + 9);

    FixedWrapper s = MSExpUnchecked(1);

    if (prec > CONSTANTS.BestE[FixedWrapper::ROUNDING].second) {
        CONSTANTS.BestE[FixedWrapper::ROUNDING].first = s;
        CONSTANTS.BestE[FixedWrapper::ROUNDING].second = prec;
    }
    FixedWrapper::SetCalcPrecision(originalprecision);
    return +s;
}

FixedWrapper FixedWrapper::GetE(uint32_t prec = FixedWrapper::CALCPRECISION) {
    if (prec > CONSTANTS.BestE[FixedWrapper::ROUNDING].second) {
        CalcE(prec);
    }
    return +CONSTANTS.BestE[FixedWrapper::ROUNDING].first;
}

FixedWrapper FixedWrapper::MSFasterBinExpUnchecked(FixedWrapper x) {
    if (x.IsPInf() || x.IsNan()) {
        return x;
    }
    if (x.IsNInf()) {
        return CONSTANTS.zero;
    }
    if (x > CONSTANTS.onetenth || x < -CONSTANTS.onetenth) {
        x.FastShift(1);
        FixedWrapper res = MSFasterBinExpUnchecked(x);
        FixedWrapper res2 = res * res;
        FixedWrapper res4 = res2 * res2;
        return (res4 * res4 * res2);
    } else {
        return MSExpUnchecked(x);
    }
}

FixedWrapper FixedWrapper::MSLogUnchecked(FixedWrapper x) {
    if (x <= CONSTANTS.zero) {
        return -FixedWrapper::GetInf();
    }
    x -= CONSTANTS.one;
    uint32_t originalprecision = FixedWrapper::CALCPRECISION;
    FixedWrapper::SetCalcPrecision(originalprecision * 11 / 10 + 9);
    FixedWrapper i = 1;
    FixedWrapper lasts = 0;
    FixedWrapper s = x;
    FixedWrapper fact = 1;
    FixedWrapper num = x;
    FixedWrapper antinum = x;
    FixedWrapper delta = 0;
    bool sign = false;
    while (s != lasts) {
        lasts = s;
        i += CONSTANTS.one;
        num *= x;
        FixedWrapper::FlipRounding();
        antinum *= x;
        FixedWrapper::FlipRounding();

        sign ^= true;
        if (sign) {
            FixedWrapper::FlipRounding();
            delta = antinum / i;
            FixedWrapper::FlipRounding();
            s -= delta;
        } else {
            delta = num / i;
            s += delta;
        }
    }
    FixedWrapper::SetCalcPrecision(originalprecision);
    +s;
    switch (FixedWrapper::ROUNDING) {
    case ROUND_ABS_DOWN:
        s.base_ -= s.Sig();
        break;
    case ROUND_ABS_UP:
        s.base_ += s.Sig();
        break;
    case ROUND_CEIL:
        s.base_ += 1;
        break;
    case ROUND_FLOOR:
        s.base_ -= 1;
        break;
    default:
        break;
    }
    return s;
}

FixedWrapper FixedWrapper::GetLN1p5() {
    if (FixedWrapper::CALCPRECISION > CONSTANTS.BestLn1p5[FixedWrapper::ROUNDING].second) {
        FixedWrapper s = MSLogUnchecked(CONSTANTS.oneandfivetenth);
        CONSTANTS.BestLn1p5[FixedWrapper::ROUNDING].first = s;
        CONSTANTS.BestLn1p5[FixedWrapper::ROUNDING].second = FixedWrapper::CALCPRECISION;

    }
    return +CONSTANTS.BestLn1p5[FixedWrapper::ROUNDING].first;
}

FixedWrapper FixedWrapper::MSFasterBinLogUnchecked(FixedWrapper x) {
    if (x.IsPInf() || x.IsNan()) {
        return x;
    }
    if (x <= CONSTANTS.zero) {
        return -FixedWrapper::GetInf();
    }
    FixedWrapper e = GetE();
    if (x > e) {
        FixedWrapper res = MSFasterBinLogUnchecked(x / e);
        return res + CONSTANTS.one;
    } else if (x > CONSTANTS.oneandfivetenth) {
        FixedWrapper res = MSFasterBinLogUnchecked(x / CONSTANTS.oneandfivetenth);
        FixedWrapper ln1p5 = GetLN1p5();
        return res + ln1p5;
    } else {
        if (x < CONSTANTS.one) {
            return -MSFasterBinLogUnchecked(CONSTANTS.one / x);
        }
        return MSLogUnchecked(x);
    }
}

FixedWrapper FixedWrapper::MSSinUnchecked(FixedWrapper x) {
    FixedWrapper pi = GetPi(); // 1
    FixedWrapper twopi = pi * CONSTANTS.two;
    FixedWrapper k = x.div(twopi);
    x -= twopi * k;
    if (x > pi) {
        x -= twopi;
    }
    if (x < -pi) {
        x += twopi;
    }
    FixedWrapper xs = x * x;
    FixedWrapper i = 1;
    FixedWrapper lasts = 0;
    FixedWrapper s = x;
    FixedWrapper fact = 1;
    FixedWrapper num = x;
    FixedWrapper delta = 0;
    bool sign = false;
    while (s != lasts) {
        lasts = s;
        i += CONSTANTS.two;
        fact *= i * (i + CONSTANTS.mone);
        num *= xs;
        delta = num / fact;
        sign ^= true;
        if (sign) {
            s -= delta;
        } else {
            s += delta;
        }
    }
    switch (FixedWrapper::ROUNDING) {
    case ROUND_ABS_DOWN:
        s.base_ -= s.Sig();
        break;
    case ROUND_ABS_UP:
        s.base_ += s.Sig();
        break;
    case ROUND_CEIL:
        s.base_ += 1;
        break;
    case ROUND_FLOOR:
        s.base_ -= 1;
        break;
    default:
        break;
    }
    return +s.clamp(CONSTANTS.mone, CONSTANTS.one);
}

FixedWrapper FixedWrapper::MSCosUnchecked(FixedWrapper x) {
    FixedWrapper pi = GetPi();
    FixedWrapper twopi = pi * CONSTANTS.two;
    FixedWrapper k = x.div(twopi);
    x -= twopi * k;
    if (x > pi) {
        x -= twopi;
    }
    if (x < -pi) {
        x += twopi;
    }
    FixedWrapper xs = x * x;
    FixedWrapper i = 0;
    FixedWrapper lasts = 0;
    FixedWrapper s = 1;
    FixedWrapper fact = 1;
    FixedWrapper num = 1;
    FixedWrapper delta = 0;
    bool sign = false;
    while (s != lasts) {
        lasts = s;
        i += CONSTANTS.two;
        fact *= i * (i + CONSTANTS.mone);
        num *= xs;
        delta = num / fact;
        sign ^= true;
        if (sign) {
            s -= delta;
        } else {
            s += delta;
        }
    }
    switch (FixedWrapper::ROUNDING) {
    case ROUND_ABS_DOWN:
        s.base_ -= s.Sig();
        break;
    case ROUND_ABS_UP:
        s.base_ += s.Sig();
        break;
    case ROUND_CEIL:
        s.base_ += 1;
        break;
    case ROUND_FLOOR:
        s.base_ -= 1;
        break;
    default:
        break;
    }
    return +s.clamp(CONSTANTS.mone, CONSTANTS.one);
}

FixedWrapper FixedWrapper::MSArctgUnchecked(FixedWrapper x) {
    x = x.clamp(CONSTANTS.mone, CONSTANTS.one);
    FixedWrapper pi = GetPi();
    FixedWrapper mpi = -pi;
    FixedWrapper s = x;
    if (x == CONSTANTS.one) {
        s = pi / CONSTANTS.four;
    } else if (x == CONSTANTS.mone) {
        s = mpi / CONSTANTS.four;
    } else {
        FixedWrapper xs = x * x;
        FixedWrapper i = 1;
        FixedWrapper lasts = 0;
        FixedWrapper num = x;
        FixedWrapper delta = 0;
        bool sign = false;
        while (s != lasts) {
            lasts = s;
            i += CONSTANTS.two;
            num *= xs;
            delta = num / i;
            sign ^= true;
            if (sign) {
                s -= delta;
            } else {
                s += delta;
            }
        }
    }
    switch (FixedWrapper::ROUNDING) {
    case ROUND_ABS_DOWN:
        s.base_ -= s.Sig();
        break;
    case ROUND_ABS_UP:
        s.base_ += s.Sig();
        break;
    case ROUND_CEIL:
        s.base_ += 1;
        break;
    case ROUND_FLOOR:
        s.base_ -= 1;
        break;
    default:
        break;
    }
    return +s.clamp(mpi, pi);
}

FixedWrapper FixedWrapper::MSSqrt(FixedWrapper x) {
    return x.MSPow(CONSTANTS.fivetenth);
}

void FixedWrapper::SetupCORDIC() {
    Rounding rndg = FixedWrapper::ROUNDING;
    uint32_t precision = FixedWrapper::CALCPRECISION;
    for (uint32_t i = 0; i < 4; ++i) {
        // Calculate K
        FixedWrapper::SetRounding(Rounding(i));
        std::vector<FixedWrapper> DK = CONSTANTS.CORDIC_DK10[i];
        std::vector<FixedWrapper> K2 = CONSTANTS.CORDIC_K2[i];
        std::vector<FixedWrapper> DK2 = CONSTANTS.CORDIC_DK2[i];
        if (K2.size() == 0 || CONSTANTS.CORDIC_MAX_PRECISION < precision) {
            CONSTANTS.CORDIC_ATAN10[i] = {};
            CONSTANTS.CORDIC_ATAN2[i] = {};
            FixedWrapper inversetwopow = CONSTANTS.one;
            FixedWrapper doubleinversetwopow = CONSTANTS.one;
            FixedWrapper inversetenpow = CONSTANTS.one;
            FixedWrapper doubleinversetenpow = CONSTANTS.one;
            FixedWrapper DKK = CONSTANTS.one;
            FixedWrapper KK2 = CONSTANTS.one;
            FixedWrapper DKK2 = CONSTANTS.one;
            
            // CORDIC10
            for (uint32_t j = 0; j < precision; ++j) {
                FixedWrapper::FlipRounding();
                DKK = MSSqrt(doubleinversetenpow + CONSTANTS.one);
                doubleinversetenpow.FastShift(2);
                FixedWrapper::FlipRounding();
                DK.push_back(CONSTANTS.one / DKK);

                CONSTANTS.CORDIC_P10i[i].push_back(inversetenpow);
                CONSTANTS.CORDIC_ATAN10[i].push_back(MSArctgUnchecked(inversetenpow));
                inversetenpow.FastShift(1);
            }
            //
            // CORDIC2
            uint32_t cordic2_max_iter = precision * 10 / 3 + 1;
            for (uint32_t j = 0; j < cordic2_max_iter; ++j) {
                FixedWrapper::FlipRounding();
                DKK2 = MSSqrt(doubleinversetwopow + CONSTANTS.one);
                KK2 *= DKK2;
                doubleinversetwopow.FastDiv4();
                FixedWrapper::FlipRounding();
                K2.push_back(CONSTANTS.one / KK2);
                DK2.push_back(CONSTANTS.one / DKK2);

                CONSTANTS.CORDIC_P2i[i].push_back(inversetwopow);
                CONSTANTS.CORDIC_ATAN2[i].push_back(MSArctgUnchecked(inversetwopow));
                inversetwopow.FastDiv2();
            }
            //
            CONSTANTS.CORDIC_DK10[i] = DK;
            CONSTANTS.CORDIC_K2[i] = K2;
            CONSTANTS.CORDIC_DK2[i] = DK2;
        }
    }
    CONSTANTS.CORDIC_MAX_PRECISION = precision;
    FixedWrapper::SetRounding(rndg);
}

std::pair<FixedWrapper, FixedWrapper> FixedWrapper::CORDIC10(FixedWrapper a) {
    if (a.IsInf() || a.IsNan()) {
        return { a, a };
    }

    uint32_t precision = FixedWrapper::CALCPRECISION;
    Rounding rndg = FixedWrapper::ROUNDING;
    Rounding antirndg = Rounding(rndg - rndg % 2 + (rndg + 1) % 2);
    FixedWrapper pi = GetPi(); // 1
    FixedWrapper twopi = pi;
    twopi.FastMul2();
    FixedWrapper halfpi = pi;
    halfpi.FastDiv2();
    FixedWrapper k = a.div(twopi);
    a -= twopi * k;
    if (a > pi) {
        a -= twopi;
    }
    if (a < -pi) {
        a += twopi;
    }
    bool mirror = false;
    if (a > halfpi) {
        a = pi - a;
        mirror = true;
    }
    if (a < -halfpi) {
        a = -pi - a;
        mirror = true;
    }
    FixedWrapper K = CONSTANTS.one;
    FixedWrapper theta = CONSTANTS.zero;
    FixedWrapper x = CONSTANTS.one;
    FixedWrapper y = CONSTANTS.zero;
    FixedWrapper antix = CONSTANTS.one;
    FixedWrapper antiy = CONSTANTS.zero;
    uint32_t i = 0;
    uint32_t last = 0;
    while (true) {
        ++i;
        FixedWrapper sigma = theta < a ? CONSTANTS.one : CONSTANTS.mone;
        FixedWrapper delta = (a - theta).abs();
        uint32_t ll = last, rr = precision + 1, mm;
        while (ll + 1 < rr) {
            mm = (ll + rr) / 2;
            if (CONSTANTS.CORDIC_ATAN10[rndg][mm - 1] > delta) {
                ll = mm;
            } else {
                rr = mm;
            }
        }
        if (rr == precision + 1)
            break;
        
        if (rr > 1) {
            FixedWrapper::FlipRounding();
            FixedWrapper newtheta1 = theta + sigma * CONSTANTS.CORDIC_ATAN10[antirndg][rr - 1];
            FixedWrapper newtheta2 = theta + sigma * CONSTANTS.CORDIC_ATAN10[antirndg][rr - 2];
            if ((a - newtheta2).abs() < (a - newtheta1).abs()) {
                --rr;
            }
            FixedWrapper::FlipRounding();
        }

        last = rr - 1;
        K *= CONSTANTS.CORDIC_DK10[rndg][last];
        FixedWrapper::FlipRounding();
        theta += sigma * CONSTANTS.CORDIC_ATAN10[antirndg][last];
        FixedWrapper antidltx = sigma * antiy * CONSTANTS.CORDIC_P10i[antirndg][last];
        FixedWrapper antidlty = sigma * antix * CONSTANTS.CORDIC_P10i[antirndg][last];
        FixedWrapper::FlipRounding();

        FixedWrapper dltx = sigma * y * CONSTANTS.CORDIC_P10i[antirndg][last];
        FixedWrapper dlty = sigma * x * CONSTANTS.CORDIC_P10i[antirndg][last];
        x -= antidltx;
        y += dlty;

        FixedWrapper::FlipRounding();
        antix -= dltx;
        antiy += antidlty;
        FixedWrapper::FlipRounding();
    }
    if (mirror)
        return { (-x * K).clamp(CONSTANTS.mone, CONSTANTS.one), (y * K).clamp(CONSTANTS.mone, CONSTANTS.one) };
    else
        return { (x * K).clamp(CONSTANTS.mone, CONSTANTS.one), (y * K).clamp(CONSTANTS.mone, CONSTANTS.one) };
}

std::pair<FixedWrapper, FixedWrapper> FixedWrapper::CORDIC2(FixedWrapper a) {
    if (a.IsInf() || a.IsNan()) {
        return { a, a };
    }

    uint32_t precision = FixedWrapper::CALCPRECISION;
    Rounding rndg = FixedWrapper::ROUNDING;
    Rounding antirndg = Rounding(rndg - rndg % 2 + (rndg + 1) % 2);
    FixedWrapper pi = GetPi(); // 1
    FixedWrapper twopi = pi;
    twopi.FastMul2();
    FixedWrapper halfpi = pi;
    halfpi.FastDiv2();
    FixedWrapper k = a.div(twopi);
    a -= twopi * k;
    if (a > pi) {
        a -= twopi;
    }
    if (a < -pi) {
        a += twopi;
    }
    bool mirror = false;
    if (a > halfpi) {
        a = pi - a;
        mirror = true;
    }
    if (a < -halfpi) {
        a = -pi - a;
        mirror = true;
    }
    FixedWrapper K = CONSTANTS.one;
    FixedWrapper theta = CONSTANTS.zero;
    FixedWrapper x = CONSTANTS.one;
    FixedWrapper y = CONSTANTS.zero;
    FixedWrapper antix = CONSTANTS.one;
    FixedWrapper antiy = CONSTANTS.zero;

    uint32_t cordic2_max_iter = precision * 10 / 3 + 1;
    for (uint32_t i = 0; i < cordic2_max_iter; ++i) {
        FixedWrapper sigma = theta < a ? CONSTANTS.one : CONSTANTS.mone;

        FixedWrapper::FlipRounding();
        theta += sigma * CONSTANTS.CORDIC_ATAN2[antirndg][i];
        FixedWrapper antidltx = sigma * antiy * CONSTANTS.CORDIC_P2i[rndg][i];
        FixedWrapper antidlty = sigma * antix * CONSTANTS.CORDIC_P2i[rndg][i];
        FixedWrapper::FlipRounding();

        FixedWrapper dltx = sigma * y * CONSTANTS.CORDIC_P2i[rndg][i];
        FixedWrapper dlty = sigma * x * CONSTANTS.CORDIC_P2i[rndg][i];
        x -= antidltx;
        y += dlty;

        FixedWrapper::FlipRounding();
        antix -= dltx;
        antiy += antidlty;
        FixedWrapper::FlipRounding();
    }
    K = CONSTANTS.CORDIC_K2[antirndg][cordic2_max_iter - 1];
    if (mirror)
        return { (-x * K).clamp(CONSTANTS.mone, CONSTANTS.one), (y * K).clamp(CONSTANTS.mone, CONSTANTS.one) };
    else
        return { (x * K).clamp(CONSTANTS.mone, CONSTANTS.one), (y * K).clamp(CONSTANTS.mone, CONSTANTS.one) };
}

FixedWrapper FixedWrapper::A_CORDIC2(std::pair<FixedWrapper, FixedWrapper> a) {
    if (a.first.IsInf() || a.first.IsNan()) {
        return a.first;
    }
    if (a.second.IsInf() || a.second.IsNan()) {
        return a.second;
    }
    bool mirror = false;
    if (a.first < 0) {
        mirror = true;
        a.first = -a.first;
    }
    uint32_t precision = FixedWrapper::CALCPRECISION;
    Rounding rndg = FixedWrapper::ROUNDING;
    Rounding antirndg = Rounding(rndg - rndg % 2 + (rndg + 1) % 2);

    FixedWrapper cos = a.first;
    FixedWrapper sin = a.second;

    FixedWrapper K = CONSTANTS.one;
    FixedWrapper theta = CONSTANTS.zero;
    FixedWrapper x = CONSTANTS.one;
    FixedWrapper y = CONSTANTS.zero;
    FixedWrapper antix = CONSTANTS.one;
    FixedWrapper antiy = CONSTANTS.zero;

    uint32_t cordic2_max_iter = precision * 10 / 3 + 1;
    for (uint32_t i = 0; i < cordic2_max_iter; ++i) {
        K = i == 0 ? 1 : CONSTANTS.CORDIC_K2[rndg][i - 1];
        FixedWrapper sigma = y * K < sin ? CONSTANTS.one : CONSTANTS.mone;

        FixedWrapper::FlipRounding();
        theta += sigma * CONSTANTS.CORDIC_ATAN2[antirndg][i];
        FixedWrapper antidltx = sigma * antiy * CONSTANTS.CORDIC_P2i[rndg][i];
        FixedWrapper antidlty = sigma * antix * CONSTANTS.CORDIC_P2i[rndg][i];
        FixedWrapper::FlipRounding();

        FixedWrapper dltx = sigma * y * CONSTANTS.CORDIC_P2i[rndg][i];
        FixedWrapper dlty = sigma * x * CONSTANTS.CORDIC_P2i[rndg][i];
        x -= antidltx;
        y += dlty;

        FixedWrapper::FlipRounding();
        antix -= dltx;
        antiy += antidlty;
        FixedWrapper::FlipRounding();
    }
    FixedWrapper pi = GetPi();
    if (mirror)
        theta = pi - theta;
    FixedWrapper halfpi = pi;
    halfpi.FastDiv2();
    return theta.clamp(-halfpi, halfpi);
}

FixedWrapper FixedWrapper::CORDIC2_Cos(FixedWrapper a) {
    return CORDIC2(a).first;
}
FixedWrapper FixedWrapper::CORDIC2_Sin(FixedWrapper a) {
    return CORDIC2(a).second;
}
FixedWrapper FixedWrapper::CORDIC2_Tan(FixedWrapper a) {
    std::pair<FixedWrapper, FixedWrapper> res = CORDIC2(a);
    return res.second / res.first;
}
FixedWrapper FixedWrapper::CORDIC2_Cot(FixedWrapper a) {
    std::pair<FixedWrapper, FixedWrapper> res = CORDIC2(a);
    return res.first / res.second;
}
FixedWrapper FixedWrapper::CORDIC2_ACos(FixedWrapper a) {
    return A_CORDIC2({ a, MSSqrt(CONSTANTS.one - a * a) });
}
FixedWrapper FixedWrapper::CORDIC2_ASin(FixedWrapper a) {
    return A_CORDIC2({ MSSqrt(CONSTANTS.one - a * a), a });
}
FixedWrapper FixedWrapper::CORDIC2_ATan(FixedWrapper a) {
    FixedWrapper cos = CONSTANTS.one / MSSqrt(CONSTANTS.one + a * a);
    return A_CORDIC2({ cos, a * cos });
}
FixedWrapper FixedWrapper::CORDIC2_ACot(FixedWrapper a) {
    FixedWrapper sin = CONSTANTS.one / MSSqrt(CONSTANTS.one + a * a);
    return A_CORDIC2({ a * sin, sin });
}

void FixedWrapper::SetupPrecision() {
    CONSTANTS.oneandfivetenth = FixedWrapper("1.5");
    CONSTANTS.oneandonetenth = FixedWrapper("1.1");
    CONSTANTS.fivetenth = FixedWrapper("0.5");
    CONSTANTS.onetenth = FixedWrapper("0.1");
    SetupCORDIC();
}

#include <ctime>

class SimpleTimer {
public:
    static uint32_t TIMERS;
    clock_t tStart;
    SimpleTimer() : tStart(clock()) {}
    ~SimpleTimer() {
        std::cout << "Elapsed time: " << (double)(clock() - tStart) * 1000.0 / CLOCKS_PER_SEC << " ms" << std::endl;
    }
};
uint32_t SimpleTimer::TIMERS = 0;

#define TIMED SimpleTimer __simple_timer_timed__;

std::string GenString(uint32_t num) {
    std::string a;
    for (uint32_t i = 0; i < num; ++i) {
        a += (rand() % 10) + '0';
    }
    return a;
}
