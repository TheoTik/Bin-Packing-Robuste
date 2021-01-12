

#include "Data.hpp"

class MasterGCol {
public :
	double solMaster;
	vector<double>pi;
	IloEnv env;
	IloObjective objective;

public:
	virtual Solution solveMasterNotSameDev(Data const &data, vector<vector<int>> const &P, int const &maxDeviation, IloNumVarArray lambda, IloExpr objMaster,IloRangeArray constraint, IloModel modelMaster)=0;
	virtual Solution solveMasterSameDev(Data const &data, vector<vector<int>> const &P, vector<vector<int>> const &P1, int const &maxDeviation, IloNumVarArray lambdap, IloNumVarArray lambdaq, IloExpr objMaster, IloRangeArray constraint, IloModel modelMaster)=0;

};


class Master: public MasterGCol {

public:
	virtual Solution solveMasterNotSameDev(Data const &data, vector<vector<int>> const &P, int const &maxDeviation, IloNumVarArray lambda, IloExpr objMaster,IloRangeArray constraint, IloModel modelMaster);
	virtual Solution solveMasterSameDev(Data const &data, vector<vector<int>> const &P, vector<vector<int>> const &P1, int const &maxDeviation, IloNumVarArray lambdap, IloNumVarArray lambdaq, IloExpr objMaster, IloRangeArray constraint, IloModel modelMaster);

};



