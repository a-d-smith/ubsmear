#ifndef UBSMEAR_UBSMEAR_OBJECTS_UBXSECMETA
#define UBSMEAR_UBSMEAR_OBJECTS_UBXSECMETA

namespace ubsmear
{

/**
* @brief A class to hold metadata about a 1D differential cross-section measurement
*/
class UBXSecMeta
{
    public:

        /**
        * @brief Constructor
        *
        * @param binEdges the edges of all bins (including any underflow or overflow bins)
        * @param hasUnderflow if there is an underflow bin
        * @param hasOverflow if there is an overflow bin
        * @param scaleByBinWidth if the cross-section corrsponds to a an event rate that's scaled by the bin width
        */
        UBXSecMeta(const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth);

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

        /**
        * @brief Get if bin width scaling is used
        *
        * @return boolean, true if bin width scaling is used
        */
        bool IsScaledByBinWidth() const;

        /**
        *  @brief  Get the bin edges
        *
        *  @return the bin edges
        */
        std::vector<float> GetBinEdges() const;

        /**
        *  @brief  Get the bin widths as a column vector that were used in the cross-section calculation. Note, if scaleByBinWidth is false,
        *          then all bin widths are set to unity.
        *
        *  @return return the bin widths
        */
        UBMatrix GetBinWidths() const;

    private:
        size_t             m_nBins;           ///< The number of bins including any underflow or overflow bins
        std::vector<float> m_binEdges;        ///< The bin edges
        bool               m_hasUnderflow;    ///< If there is an underflow bin
        bool               m_hasOverflow;     ///< If there is an overflow bin
        bool               m_scaleByBinWidth; ///< If the cross-section bins have been scaled by their width
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBXSecMeta::UBXSecMeta(const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth) :
    m_binEdges(binEdges),
    m_hasUnderflow(hasUnderflow),
    m_hasOverflow(hasOverflow),
    m_scaleByBinWidth(scaleByBinWidth)
{
    // Check if the number of bin edges is enough for the underflow, overflow and at least one more bin
    if (static_cast<int>(binEdges.size()) - 1 <= (hasUnderflow ? 1 : 0) + (hasOverflow ? 1 : 0))
        throw std::invalid_argument("UBXSecMeta::UBXSecMeta - insufficient bin edges supplied");

    // Store the number of bins
    m_nBins = binEdges.size() - 1;

    // Check that the bin edges are valid
    for (size_t iBin = 0; iBin < m_nBins; ++iBin)
    {
        const auto lowerEdge = binEdges.at(iBin);
        const auto upperEdge = binEdges.at(iBin + 1);
        const auto width = upperEdge - lowerEdge;

        if (width <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("UBXSecMeta::UBXSecMeta - supplied bin " + std::to_string(iBin) + " has width " + std::to_string(width));
    }
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

// -----------------------------------------------------------------------------------------------------------------------------------------

inline bool UBXSecMeta::IsScaledByBinWidth() const
{
    return m_scaleByBinWidth;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<float> UBXSecMeta::GetBinEdges() const
{
    return m_binEdges;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBXSecMeta::GetBinWidths() const
{
    if (!m_scaleByBinWidth)
    {
        // If the cross-section wasn't scaled by bin width, then just return 1.f for the width of each bin
        return ubsmear::UBMatrix(std::vector<float>(m_nBins, 1.f), m_nBins, 1);
    }

    std::vector<float> elements;
    for (unsigned int iBin = 1; iBin < m_binEdges.size(); ++iBin)
    {
        const auto lowerEdge = m_binEdges.at(iBin - 1);
        const auto upperEdge = m_binEdges.at(iBin);
        const auto width = upperEdge - lowerEdge;

        if (width <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("UBXSecMeta::GetBinWidth - Found a bin with zero or negative width");

        elements.push_back(width);
    }

    return ubsmear::UBMatrix(elements, m_nBins, 1);
}

} // namespace ubsmear

#endif
