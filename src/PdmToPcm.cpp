#include "PdmToPcm.hpp"
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <span>

namespace PdmToPcm
{
Filter::Filter()
{
    constexpr size_t filterBankStepStagger = numberOfLookupTableSteps_ / numberOfFilterBanks_;

    for (size_t i = 0; i < numberOfFilterBanks_; ++i)
    {
        banks_.at(i).step_
            = (numberOfLookupTableSteps_ - (i * filterBankStepStagger)) % numberOfLookupTableSteps_;
        banks_.at(i).accumulator_ = 0;
    }
}

auto Filter::Apply(std::span<uint8_t> const dataIn, std::span<int32_t> dataOut) -> std::size_t
{
    auto outIter = dataOut.begin();

    if (dataOut.size() > 0)
    {
        for (auto const in : dataIn)
        {
            for (auto& bank : banks_)
            {
                bank.accumulator_ += bank.step_->at(in);
                ++bank.step_;
                if (bank.step_ == lookupTable_.end())
                {
                    (*outIter) = bank.accumulator_;
                    ++outIter;
                    bank.step_ = lookupTable_.begin();
                    bank.accumulator_ = 0;

                    if (outIter == dataOut.end())
                    {
                        return dataOut.size();
                    }
                }
            }
        }
    }

    return std::distance(dataOut.begin(), outIter);
}
} // namespace PdmToPcm