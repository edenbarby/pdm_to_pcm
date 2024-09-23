#include "ChebyshevPolynomial.hpp"
#include "gtest/gtest.h"
#include <cmath>
#include <numbers>

namespace
{
template <typename T>
class LinearSpacing
{
    T const start_;
    T const delta_;
    size_t const count_;

public:
    class Iterator
    {
        T const start_;
        T const delta_;
        size_t state_;

    public:
        Iterator(T start, T delta, size_t state = 0)
            : start_{start}
            , delta_{delta}
            , state_{state}
        {
        }

        ~Iterator() = default;
        Iterator(Iterator const&) = delete;
        Iterator(Iterator&&) = delete;
        auto operator=(Iterator const&) -> Iterator& = delete;
        auto operator=(Iterator&&) -> Iterator& = delete;

        auto operator!=(Iterator const& other) const -> bool
        {
            return (this->start_ != other.start_) || (this->delta_ != other.delta_)
                || (this->state_ != other.state_);
        }

        auto operator*() const -> T
        {
            return start_ + delta_ * static_cast<T>(state_);
        }

        auto operator++() -> Iterator&
        {
            state_ += 1;
            return *this;
        }
    };

    LinearSpacing(T start, T end, size_t count)
        : start_{start}
        , delta_{(end - start) / static_cast<T>(count - 1)}
        , count_{count}
    {
    }

    ~LinearSpacing() = default;
    LinearSpacing(LinearSpacing const&) = delete;
    LinearSpacing(LinearSpacing&&) = delete;
    auto operator=(LinearSpacing const&) -> LinearSpacing& = delete;
    auto operator=(LinearSpacing&&) -> LinearSpacing& = delete;

    [[nodiscard]] auto begin() const -> LinearSpacing::Iterator
    {
        return {start_, delta_, 0};
    }

    [[nodiscard]] auto end() const -> LinearSpacing::Iterator
    {
        return {start_, delta_, count_};
    }
};
} // namespace

TEST(ChebyshevPolynomialTests, PolynomialDegreeZero)
{
    constexpr std::array<double, 1> coefficients = {1.0};
    for (auto const x : LinearSpacing(-1.0, 1.0, 10000))
    {
        ASSERT_DOUBLE_EQ(ChebyshevPolynomial::Evaluate(coefficients, x), 1.0)
            << "Evaluated at x = " << x;
    }
}

TEST(ChebyshevPolynomialTests, PolynomialDegreeOne)
{
    constexpr std::array<double, 2> coefficients = {0.0, 1.0};
    for (auto const x : LinearSpacing(-1.0, 1.0, 10000))
    {
        ASSERT_DOUBLE_EQ(ChebyshevPolynomial::Evaluate(coefficients, x), x)
            << "Evaluated at x = " << x;
    }
}

TEST(ChebyshevPolynomialTests, PolynomialDegreeTwo)
{
    constexpr std::array<double, 3> coefficients = {0.0, 0.0, 1.0};
    for (auto const x : LinearSpacing(-1.0, 1.0, 10000))
    {
        ASSERT_DOUBLE_EQ(ChebyshevPolynomial::Evaluate(coefficients, x), 2 * x * x - 1)
            << "Evaluated at x = " << x;
    }
}

TEST(ChebyshevPolynomialTests, PolynomialDegreeThree)
{
    constexpr std::array<double, 4> coefficients = {0.0, 0.0, 0.0, 1.0};
    for (auto const x : LinearSpacing(-1.0, 1.0, 10000))
    {
        ASSERT_NEAR(ChebyshevPolynomial::Evaluate(coefficients, x), 4 * x * x * x - 3 * x, 1e-15)
            << "Evaluated at x = " << x;
    }
}

TEST(ChebyshevPolynomialTests, PolynomialDegreeFour)
{
    constexpr std::array<double, 5> coefficients = {0.0, 0.0, 0.0, 0.0, 1.0};
    auto polynomial_function = [](double x) {
        auto const t0 = 1.0;
        auto const t1 = x;
        auto const t2 = 2.0 * x * t1 - t0;
        auto const t3 = 2.0 * x * t2 - t1;
        auto const t4 = 2.0 * x * t3 - t2;
        return t4;
    };

    for (auto const x : LinearSpacing(-1.0, 1.0, 10000))
    {
        ASSERT_NEAR(ChebyshevPolynomial::Evaluate(coefficients, x), polynomial_function(x), 1e-15)
            << "Evaluated at x = " << x;
    }
}

TEST(ChebyshevPolynomialTests, PolynomialDegreeFive)
{
    // Some not particualarly accurate coefficients for the sine function from
    // www.embeddedrelated.com. Just wanted to evaluate something a bit less trival.
    constexpr std::array<float, 6> sinChebyshevCoefficients
        = {0.0F, 1.1336F, 0.0F, -0.13807F, 0.0F, 0.0045584F};
    auto sin_function = [](float x) {
        constexpr auto scale = 0.5F * std::numbers::pi_v<float>;
        return std::sin(scale * x);
    };

    for (auto const x : LinearSpacing(-1.0F, 1.0F, 10000))
    {
        ASSERT_NEAR(
            ChebyshevPolynomial::Evaluate(sinChebyshevCoefficients, x),
            sin_function(x),
            1e-3F)
            << "Evaluated at x = " << x;
    }
}