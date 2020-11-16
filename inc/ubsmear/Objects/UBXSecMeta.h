#ifndef UBSMEAR_UBSMEAR_OBJECTS_UBXSECMETA
#define UBSMEAR_UBSMEAR_OBJECTS_UBXSECMETA

namespace ubsmear
{

/*
* @brief A structure to hold metadata about a 1D differential cross-section measurement
*/
struct UBXSecMeta
{
    size_t nBins;        ///< The number of bins including any underflow or overflow bins
    bool   hasUnderflow; ///< If there is an underflow bin
    bool   hasOverflow;  ///< If there is an overflow bin
};

} // namespace ubsmear

#endif
