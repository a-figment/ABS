#ifndef DATAKEYS_H_INCLUDED
#define DATAKEYS_H_INCLUDED

#include <string>
#include <vector>
#include <iostream>

struct DataKeys
{

    DataKeys(const std::vector<std::string> &s, const std::vector<uint> &maxs)
    {
        assert(s.size() == maxs.size());
        DataKey = std::make_pair(s, maxs);
        nKeys = s.size();
        Meta = std::vector<std::string>(s.size(), "null");
        CurrentKey = "";
    }

    // @brief Include @meta description of the @s
    DataKeys(const std::vector<std::string> &s, const std::vector<uint> &maxs, const std::vector<std::string> &meta)
    {
        assert(s.size() == maxs.size());
        assert(s.size() == meta.size());
        DataKey = std::make_pair(s, maxs);
        nKeys = s.size();
        Meta = meta;
        CurrentKey = "";
    }

    std::string MakeKey(const std::vector<int> &idx)
    {
        std::cout << "Making key from indices: ";
        for (auto j : idx)
            std::cout << j << " ";
        std::cout << std::endl;

        if (idx.size() != DataKey.first.size())
        {
            std::cout << "Number of keys (" << idx.size() << ") does not expected number of keys (" << nKeys << "). Check command line arguments. Exiting." << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string key = "";
        for (size_t i = 0; i < DataKey.first.size(); ++i)
        {
            if ((size_t)idx[i] > DataKey.second[i] - 1)
            {
                std::cout << "Key index (" << DataKey.first[i] << idx[i] << ") exceeds maximum key index (" << DataKey.first[i] << DataKey.second[i] << "). Exiting." << std::endl;
                exit(EXIT_FAILURE);
            }
            key += (DataKey.first[i] + std::to_string(idx[i])) + "-";
        }
        key += "-";
        CurrentKey = key;
        return key;
    }

    void Print()
    {
        std::cout << "\n========== Data key summary ==========" << std::endl;
        for (size_t i = 0; i < nKeys; ++i)
        {
            if (DataKey.second[i] == 1)
                std::cout << DataKey.first[i] << "[1] : " << Meta[i] << "." << std::endl;
            else
                std::cout << DataKey.first[i] << "[1-" << DataKey.second[i] << "] : " << Meta[i] << "." << std::endl;
        }
        std::cout << "Current key = " << CurrentKey << std::endl;
        std::cout << "======================================" << std::endl;
    }

    std::pair<std::vector<std::string>, std::vector<uint>> DataKey;
    std::string CurrentKey;
    std::vector<std::string> Meta;
    size_t nKeys;
};

#endif
