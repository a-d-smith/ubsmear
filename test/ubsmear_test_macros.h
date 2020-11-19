#ifndef UBSMEAR_TEST_MACROS
#define UBSMEAR_TEST_MACROS

#include <stdexcept>

// Try to execute an expression which we want to raise an exception - if no exception is raised, then throw!
#define UB_ASSERT_THROW(expression) {                                                                                                      \
    {                                                                                                                                      \
        bool didThrow = false;                                                                                                             \
        try                                                                                                                                \
        {                                                                                                                                  \
            expression;                                                                                                                    \
        }                                                                                                                                  \
        catch (...)                                                                                                                        \
        {                                                                                                                                  \
            didThrow = true;                                                                                                               \
        }                                                                                                                                  \
        if (!didThrow)                                                                                                                     \
        {                                                                                                                                  \
            throw std::logic_error("No throw for: "#expression);                                                                           \
        }                                                                                                                                  \
    }                                                                                                                                      \
}                                                                                                                                          \

#endif
