

#include "Data.hpp"

class SubProblemGCol {
public :
	double optSP;
	IloEnv env;
	IloObjective objective;

public:
    virtual vector<int> solveSubProblemSameDev(Data const &data, vector<double> const &pi, int const &maxDeviations, int const &numeroSP, IloNumVarArray x, IloExpr obj, IloModel modelSP)=0;
    virtual vector<int> solveSubProblemNotSameDev(Data const &data, vector<double> const &pi, int const &maxDeviations, IloNumVarArray x, IloNumVarArray lambda, IloNumVar mu, IloExpr obj, IloModel modelSP)=0;

};


class SubProblem: public SubProblemGCol {

public:

    virtual vector<int> solveSubProblemSameDev(Data const &data, vector<double> const &pi, int const &maxDeviations, int const &numeroSP, IloNumVarArray x, IloExpr obj, IloModel modelSP);
    virtual vector<int> solveSubProblemNotSameDev(Data const &data, vector<double> const &pi, int const &maxDeviations, IloNumVarArray x, IloNumVarArray lambda, IloNumVar mu, IloExpr obj, IloModel modelSP);
};




