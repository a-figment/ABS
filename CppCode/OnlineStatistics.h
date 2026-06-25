#ifndef ONLINESTATISTICS_H_INCLUDED
#define ONLINESTATISTICS_H_INCLUDED

template <typename T>
struct OnlineStatistics
{
        OnlineStatistics()
        {
                Clear();
        }

        void Clear()
        {
                n = 0;
                Mean = Variance = static_cast<T>(0);
        }

        void Update(const T &value)
        {
                long int n1 = n;
                n++;
                Mean += (value - Mean) / n;
                Variance += ((value - Mean) * (value - Mean) / n) * n1;
        }

        T getMean() const
        {
                return Mean;
        }
        T getVariance() const
        {
                return (Variance / (n - static_cast<T>(1)));
        }
        long int getNumElems() const
        {
                return n;
        }

        long int n;
        T Mean, Variance;
};

#endif
