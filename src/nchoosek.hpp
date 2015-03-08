#pragma once


namespace ifgt
{


inline int nchoosek(int n, int k)
{
    int n_k = n - k;
    if (k < n_k)
    {
        k = n_k;
        n_k = n - k;
    }

    int nchsk = 1;
    for (int i = 1; i <= n_k; ++i)
    {
        nchsk *= ++k;
        nchsk /= i;
    }

    return nchsk;
}
}
