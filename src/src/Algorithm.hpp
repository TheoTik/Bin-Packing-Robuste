
#include "Data.hpp"
#include "First_Fit.hpp"

class Algorithm {
public:
	int numeroMethod;

public:
    virtual Solution solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit Ffq)=0;

};


class Rbpp: public Algorithm {

public:

    virtual Solution solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit Ffq);

};

class RelaxRbpp: public Algorithm {

public:

    virtual Solution solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit Ffq);

};


class GColSameDev: public Algorithm {

public:

    virtual Solution solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit Ffq);

};

class GColNotSameDev: public Algorithm {

public:

    virtual Solution solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit Ffq);

};

class BranchPrice: public Algorithm {

public:

    virtual Solution solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit Ffq);

};






