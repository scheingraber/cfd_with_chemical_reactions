#include "Range2.h"

Range2::Range2 ()
{
}

Range2::Range2 (int ilow, int ihigh, int jlow, int jhigh)
    : i(Range(ilow, ihigh)), j(Range(jlow, jhigh))
{
}

Range2::Range2 (Range i, Range j)
    : i(i), j(j)
{
}

