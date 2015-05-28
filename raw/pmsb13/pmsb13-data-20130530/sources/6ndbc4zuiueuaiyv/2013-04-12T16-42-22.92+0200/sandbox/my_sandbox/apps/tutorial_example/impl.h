#ifndef IMPL_H_
#define IMPL_H_

double quadrat(double x)
{
    return x*x;
}

void iota(seqan::String<int> & result, int begin, int end)
{
    resize(result, end - begin, 0);
    for (int i = begin, k = 0; i < end; ++k, ++i)
        result[k] = i;
}

#endif

