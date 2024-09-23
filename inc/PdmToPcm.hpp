#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace PdmToPcm
{
class Filter
{
private:
    static constexpr std::size_t decimationRate_ = 72;
    static constexpr std::size_t filterLength_ = 1368;
    static constexpr double kaiserWindowBeta_ = 5.4;

    static constexpr std::size_t samplesPerLookupTableStep_ = 8;

    static_assert(filterLength_ % decimationRate_ == 0);
    static_assert(filterLength_ % samplesPerLookupTableStep_ == 0);
    static_assert(decimationRate_ % samplesPerLookupTableStep_ == 0);

    static constexpr std::size_t numberOfLookupTableSteps_ = filterLength_ / 8;

    // The lookup table exploits the shallowness of the PDM quantisation by packing 8 samples (1 bit
    // per sample) into each lookup step. This completely removes the need to perform expensive
    // multiplications and reduces the number of additions by a factor of 8.
    static constexpr std::size_t lookupTableInnerDimension = (1U << samplesPerLookupTableStep_);
    static constexpr std::size_t lookupTableOuterDimension = numberOfLookupTableSteps_;
    using LookupTableIntegerType = int16_t;
    using LookupTableInnerType = std::array<LookupTableIntegerType, lookupTableInnerDimension>;
    using LookupTableType = std::array<LookupTableInnerType, lookupTableOuterDimension>;
    static constexpr std::size_t lookupTableSize = sizeof(LookupTableType);
    inline static LookupTableType lookupTable_ /* = */;

    struct Bank
    {
        LookupTableType::const_iterator step_;
        int32_t accumulator_;
    };

    static constexpr std::size_t numberOfFilterBanks_ = filterLength_ / decimationRate_;
    std::array<Bank, numberOfFilterBanks_> banks_;

public:
    Filter();
    ~Filter() = default;
    Filter(Filter const&) = delete;
    Filter(Filter&&) = delete;
    void operator=(Filter const&) = delete;
    void operator=(Filter&&) = delete;
    auto Apply(std::span<uint8_t> dataIn, std::span<int32_t> dataOut) -> std::size_t;
};
} // namespace PdmToPcm

// #pragma once

// #include "KaiserWindow.hpp"

// namespace PdmToPcm
// {

// // template <class OutputType, size_t filterLengthDiv8>
// // constexpr auto FindScaleToFitOutputType(
// //     std::array<double, filterLengthDiv8 * 8> const& filterKernel) -> double
// // {
// //     double maxLutStepSum = 0.0;
// //     for (size_t lutStep = 0; lutStep < filterLengthDiv8; ++lutStep)
// //     {
// //         double lutStepSum = 0.0;
// //         for (size_t shift = 0; shift < 8; ++shift)
// //         {
// //             constexpr double dcOffset = 0.5;
// //             constexpr double maxFilterInput = 1 - dcOffset;

// //             size_t const filterStep = 8 * lutStep + shift;

// //             lutStepSum += maxFilterInput * filterKernel[filterStep];
// //         }

// //         if (lutStepSum > maxLutStepSum)
// //         {
// //             maxLutStepSum = lutStepSum;
// //         }
// //     }
// //     constexpr auto maxOutputValue =
// static_cast<double>(std::numeric_limits<OutputType>::max());

// //     return maxOutputValue / maxLutStepSum;
// // }

// template <class OutputType, size_t filterLengthDiv8>
// constexpr auto CreateKernelLookupTable(double const beta, bool const msb = true)
//     -> std::array<std::array<OutputType, 256>, filterLengthDiv8>
// {
//     static_assert(std::is_integral_v<OutputType>);
//     static_assert(std::is_signed_v<OutputType>);

//     constexpr size_t filterLength = filterLengthDiv8 * 8;
//     auto const filterKernel = CreateKernel<filterLength>(beta);

//     double maxLutStepSum = 0.0;
//     for (size_t lutStep = 0; lutStep < filterLengthDiv8; ++lutStep)
//     {
//         double lutStepSum = 0.0;
//         for (size_t shift = 0; shift < 8; ++shift)
//         {
//             constexpr double dcOffset = 0.5;
//             constexpr double maxFilterInput = 1 - dcOffset;

//             size_t const filterStep = 8 * lutStep + shift;

//             lutStepSum += maxFilterInput * filterKernel[filterStep];
//         }

//         if (lutStepSum > maxLutStepSum)
//         {
//             maxLutStepSum = lutStepSum;
//         }
//     }
//     constexpr auto maxOutputValue = static_cast<double>(std::numeric_limits<OutputType>::max());
//     double const scaleToFitOutputType = maxOutputValue / maxLutStepSum;

//     // double const scaleToFitOutputType
//     //     = FindScaleToFitOutputType<OutputType, filterLengthDiv8>(filterKernel);

//     std::array<std::array<OutputType, 256>, filterLengthDiv8> lut{};
//     for (size_t lutStep = 0; lutStep < filterLengthDiv8; ++lutStep)
//     {
//         for (size_t pdm = 0; pdm < 256; ++pdm)
//         {
//             double lutStepSum = 0.0;
//             for (size_t shift = 0; shift < 8; ++shift)
//             {
//                 size_t const bitIdx = msb ? (7 - shift) : shift;

//                 constexpr double dcOffset = 0.5;
//                 double const filterInput = ((pdm >> bitIdx) & 1) - dcOffset;

//                 size_t const filterStep = 8 * lutStep + shift;

//                 lutStepSum += filterInput * filterKernel[filterStep];
//             }

//             lut[lutStep][pdm] = static_cast<OutputType>(std::clamp(
//                 lutStepSum * scaleToFitOutputType,
//                 static_cast<double>(std::numeric_limits<OutputType>::min()),
//                 static_cast<double>(std::numeric_limits<OutputType>::max())));
//         }
//     }

//     return lut;
// }

// // static constexpr double Rolloff3db(double const sampleFrequency)
// // {
// //     return KaiserWindow::RollOff3dBNormalised(numberOfTaps_, beta_) * sampleFrequency;
// // }

// // static constexpr double RolloffNull(double const sampleFrequency)
// // {
// //     return KaiserWindow::RollOffNullNormalised(numberOfTaps_, beta_) * sampleFrequency;
// // }

// template <
//     typename OutputIntegerType,
//     std::size_t numberOfTaps_,
//     std::size_t decimationRate_,
//     std::size_t outputLength_,
//     int beta1000_>
// class Filter
// {
// private:
//     static constexpr double beta_ = static_cast<double>(beta1000_) / 1000.0;

//     static constexpr size_t numberOfFilterBanks_{numberOfTaps_ / decimationRate_};
//     static_assert(
//         numberOfFilterBanks_ * decimationRate_ == numberOfTaps_,
//         "numberOfTaps_ must be a multiple of decimationRate!");

//     struct Bank
//     {
//         size_t step_;
//         OutputIntegerType accumulator_;
//     };

//     std::array<Bank, numberOfFilterBanks_> banks_{};

//     static constexpr size_t tapsPerStep_{8U};
//     static constexpr size_t numberOfSteps_{numberOfTaps_ / tapsPerStep_};
//     static_assert((1 << tapsPerStep_) == 256);
//     static_assert(numberOfTaps_ % tapsPerStep_ == 0, "numberOfTaps_ must be a multiple of 8!");
//     static_assert(decimationRate_ % tapsPerStep_ == 0, "decimationRate must be a multiple of
//     8!");

//     using LutIntegerType = int16_t;
//     using LutType = std::array<std::array<LutIntegerType, 256>, numberOfSteps_>;

//     inline static LutType lut_
//         = KaiserWindow::CreateKernelLookupTable<LutIntegerType, numberOfSteps_>(beta_);

// public:
//     Filter()
//     {
//         constexpr size_t filterBankStepStagger = numberOfSteps_ / numberOfFilterBanks_;

//         for (size_t i = 0; i < numberOfFilterBanks_; ++i)
//         {
//             banks_[i].step_ = (numberOfSteps_ - (i * filterBankStepStagger)) % numberOfSteps_;
//             banks_[i].accumulator_ = 0;
//         }
//     }

//     Filter(Filter const&) = delete;
//     Filter(Filter&&) = delete;
//     void operator=(Filter const&) = delete;
//     void operator=(Filter&&) = delete;

//     static constexpr size_t inputLength_ = (outputLength_ * decimationRate_) / tapsPerStep_;
//     using InputType = std::array<uint8_t, inputLength_>;
//     using OutputType = std::array<OutputIntegerType, outputLength_>;

//     void ApplySlow(InputType const& dataIn, OutputType& dataOut)
//     {
//         // So we don't need to worry about going off the end of dataOut.
//         static_assert((inputLength_ * numberOfFilterBanks_) / numberOfSteps_ == outputLength_);

//         auto outIter = dataOut.begin();

//         for (auto const in : dataIn)
//         {
//             for (auto& bank : banks_)
//             {
//                 bank.accumulator_ += lut_[bank.step_][in];
//                 ++bank.step_;
//                 if (bank.step_ == numberOfSteps_)
//                 {
//                     (*outIter) = bank.accumulator_;
//                     ++outIter;
//                     bank.step_ = 0;
//                     bank.accumulator_ = 0;
//                 }
//             }
//         }
//     }

//     void Apply(InputType const& dataIn, OutputType& dataOut)
//     {
//         constexpr size_t loopUnrollMultiple = 5;
//         static_assert(numberOfFilterBanks_ == 4);
//         static_assert(numberOfSteps_ % loopUnrollMultiple == 0);
//         static_assert(inputLength_ > loopUnrollMultiple);

//         size_t inIdx = 0;
//         size_t outIdx = 0;
//         auto step0 = banks_[0].step_;
//         auto step1 = banks_[1].step_;
//         auto step2 = banks_[2].step_;
//         auto step3 = banks_[3].step_;
//         auto sum0 = banks_[0].accumulator_;
//         auto sum1 = banks_[1].accumulator_;
//         auto sum2 = banks_[2].accumulator_;
//         auto sum3 = banks_[3].accumulator_;

//         while (inIdx < inputLength_)
//         {
//             auto const data0 = dataIn[inIdx];
//             auto const data1 = dataIn[inIdx + 1];
//             auto const data2 = dataIn[inIdx + 2];
//             auto const data3 = dataIn[inIdx + 3];
//             auto const data4 = dataIn[inIdx + 4];
//             inIdx += loopUnrollMultiple;

//             sum0 += lut_[step0][data0] + lut_[step0 + 1][data1] + lut_[step0 + 2][data2]
//                   + lut_[step0 + 3][data3] + lut_[step0 + 4][data4];
//             sum1 += lut_[step1][data0] + lut_[step1 + 1][data1] + lut_[step1 + 2][data2]
//                   + lut_[step1 + 3][data3] + lut_[step1 + 4][data4];
//             sum2 += lut_[step2][data0] + lut_[step2 + 1][data1] + lut_[step2 + 2][data2]
//                   + lut_[step2 + 3][data3] + lut_[step2 + 4][data4];
//             sum3 += lut_[step3][data0] + lut_[step3 + 1][data1] + lut_[step3 + 2][data2]
//                   + lut_[step3 + 3][data3] + lut_[step3 + 4][data4];

//             step0 += loopUnrollMultiple;
//             step1 += loopUnrollMultiple;
//             step2 += loopUnrollMultiple;
//             step3 += loopUnrollMultiple;

//             if (step0 == numberOfSteps_)
//             {
//                 dataOut[outIdx++] = sum0;
//                 step0 = 0;
//                 sum0 = 0;
//             }
//             if (step1 == numberOfSteps_)
//             {
//                 dataOut[outIdx++] = sum1;
//                 step1 = 0;
//                 sum1 = 0;
//             }
//             if (step2 == numberOfSteps_)
//             {
//                 dataOut[outIdx++] = sum2;
//                 step2 = 0;
//                 sum2 = 0;
//             }
//             if (step3 == numberOfSteps_)
//             {
//                 dataOut[outIdx++] = sum3;
//                 step3 = 0;
//                 sum3 = 0;
//             }
//         }

//         banks_[0].step_ = step0;
//         banks_[1].step_ = step1;
//         banks_[2].step_ = step2;
//         banks_[3].step_ = step3;
//         banks_[0].accumulator_ = sum0;
//         banks_[1].accumulator_ = sum1;
//         banks_[2].accumulator_ = sum2;
//         banks_[3].accumulator_ = sum3;
//     }
// };
// } // namespace PdmToPcm