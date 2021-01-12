#pragma once

using namespace std;

#include <vector>
#include <string>

// ----------------------------------------------------------------------------
//  La classe DataULS stocke les données de l'ULS
// ----------------------------------------------------------------------------
class Data {

public:

    int sizeBin;
    int numberOfItems; // number of items
    vector<int> sizeObj; // size of the obect
    vector<int> deviationObj;
    int deviation;

    Data(string const & filePath);

};


class Solution {

public:
    vector<vector<double>> bins;
    double opt;
    Solution();

};

