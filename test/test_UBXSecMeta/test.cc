#include "ubsmear.h"
#include "ubsmear_test_macros.h"

#include <assert.h>

using namespace ubsmear;

int main()
{
    // Try to construct some UBXSecMeta with invalid input parameters
    // There needs to be at least one bin that's not underflow/overflow
    UB_ASSERT_THROW( UBXSecMeta(0, false, false) );
    UB_ASSERT_THROW( UBXSecMeta(0, false, true) );
    UB_ASSERT_THROW( UBXSecMeta(1, true, false) );
    UB_ASSERT_THROW( UBXSecMeta(1, true, true) );
    UB_ASSERT_THROW( UBXSecMeta(2, true, true) );

    return 0;
}

