#ifndef UBSMEAR_UBSMEAR_OBJECTS_UBXSECMETA
#define UBSMEAR_UBSMEAR_OBJECTS_UBXSECMETA

namespace ubsmear
{

/*
* @brief A class to hold metadata about a 1D differential cross-section measurement
*/
class UBXSecMeta
{
    public:

        /**
        * @brief Constructor
        *
        * @param nBins the number of bins (including any underflow or overflow bins)
        * @param hasUnderflow if there is an underflow bin
        * @param hasOverflow if there is an overflow bin
        */
        UBXSecMeta(const size_t nBins, const bool hasUnderflow, const bool hasOverflow);

        /**
        * @brief Get the number of bins
        *
        * @return the number of bins
        */
        size_t GetNBins() const;

        /**
        * @brief Get if there is an underflow bin
        *
        * @return boolean, true if there is an underflow bin
        */
        bool HasUnderflow() const;

        /**
        * @brief Get if there is an overflow bin
        *
        * @return boolean, true if there is an overflow bin
        */
        bool HasOverflow() const;

        /**
        * @brief Determine if an input bin index (counting from zero) is an underflow or overflow bin
        *
        * @param iBin the bin index
        *
        * @return boolean, true if the bin index is underflow or overflow, false otherwise
        */
        bool IsUnderOverflowBin(const size_t iBin) const;

    private:
        size_t m_nBins;        ///< The number of bins including any underflow or overflow bins
        bool   m_hasUnderflow; ///< If there is an underflow bin
        bool   m_hasOverflow;  ///< If there is an overflow bin
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBXSecMeta::UBXSecMeta(const size_t nBins, const bool hasUnderflow, const bool hasOverflow) :
    m_nBins(nBins),
    m_hasUnderflow(hasUnderflow),
    m_hasOverflow(hasOverflow)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline size_t UBXSecMeta::GetNBins() const
{
    return m_nBins;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline bool UBXSecMeta::HasUnderflow() const
{
    return m_hasUnderflow;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline bool UBXSecMeta::HasOverflow() const
{
    return m_hasOverflow;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline bool UBXSecMeta::IsUnderOverflowBin(const size_t iBin) const
{
    if (iBin >= m_nBins)
        throw std::out_of_range("UBXSecMeta::IsUnderOverflowBin - Input bin index is out of range");

    // Check if underflow bin
    if (m_hasUnderflow && iBin == 0)
        return true;

    // Check if overflow bin
    if (m_hasOverflow && iBin == m_nBins - 1)
        return true;

    return false;
}

} // namespace ubsmear

#endif
