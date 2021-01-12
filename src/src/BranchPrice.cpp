#include <iostream>
#include <ilcplex/ilocplex.h>
#include <cmath>;

#include "Algorithm.hpp"
#include "Master.hpp"
#include "SubProblem.hpp"


vector<vector<double>> transpos(vector<vector<double> > &b)
{

	vector<vector<double> > trans_vec(b[0].size(), vector<double>());

	if (b.size() == 0)
		return trans_vec;

	for (int i = 0; i < b.size(); i++){
		for (int j = 0; j < b[i].size(); j++){
			trans_vec[j].push_back(b[i][j]);
		}
	}

	return trans_vec;
}

vector<vector<double>> keepSol1(IloRangeArray constraint, vector<double> solution, IloNumVarArray lambda){
	vector<vector<double>> finalSol;
	vector<double> vect(constraint.getSize());
	for(int i = 0; i < solution.size(); i++){
		finalSol.push_back(vect);
	}


	for(int i = 0; i < constraint.getSize(); i++){

		for (IloExpr::LinearIterator it = IloExpr(constraint[i].getExpr()).getLinearIterator(); it.ok();++it){

			int posVarInLambda = lambda.find(it.getVar()); // find the position in the list of variables

			if(solution[posVarInLambda] != 0){

				finalSol[posVarInLambda][i] = solution[posVarInLambda];
			}
		}
	}

	vector<vector<double>> finalSol1 = transpos(finalSol);
	vector<vector<double>> finalSol2;

	for(int i = 0; i < constraint.getSize(); i++){
		if(!finalSol1[i].empty()){
			finalSol2.push_back(finalSol1[i]);
		}
	}

	return transpos(finalSol2);
}

void affiche(vector<vector<int>> P){
	for(int i = 0; i < P.size(); i++){
		for(int j = 0; j < P[0].size(); j++){
			cout << P[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void affiche(vector<double> gamma){
	for(int i = 0; i < gamma.size(); i++){
		cout << gamma[i] << " ";
	}
	cout << endl;
}

void affiche(IloModel model){
	cout << model << endl;
}

void launchBP(MasterGCol* m, Data const &data, Solution solution, First_Fit const &Ff, int const &maxDeviation, IloNumVarArray lambda, IloExpr objMaster,IloRangeArray constraint, IloModel modelMaster, int level, SubProblemGCol* sp, IloNumVarArray x, IloNumVarArray delta, IloNumVar mu, IloExpr objSP, IloModel modelSP, vector<int> solSP, IloNumArray a, IloNumArray b, int countColumns){

	// add column to the master
	IloNumColumnArray col(m->env,1);
	col[0] = m->objective(1);
	for(IloInt i = 0; i < constraint.getSize(); i++){
		col[0] += constraint[i](solSP[i]);
	}
	solSP.clear();

	IloNumVarArray lambdavar(m->env, col, a, b);
	lambda.add(lambdavar);
	countColumns++; // count number of columns

	// solve the master and keep duals in m->vals
	solution.bins.clear();
	IloCplex solverMaster(modelMaster);
	solverMaster.solve();

	int item = -1;
	for(int i = 0; i < lambda.getSize(); i++){
		if(solverMaster.getValue(lambda[i]) > 0 && solverMaster.getValue(lambda[i]) < 1){
			item = i;
			break;
		}
	}

	if(item != -1){
		IloNumArray vals(m->env);

		solverMaster.getDuals(vals, constraint);
		cout << "Duals : " << vals << endl;

		int solMaster = solverMaster.getObjValue();
		cout << "Solution master : " << solMaster << endl;

		// solve subproblem and modify pi
		objSP.setLinearCoefs(x, vals);
		IloCplex solverSP(modelSP);
		solverSP.solve();

		// check solution
		double currentObjSP = solverSP.getObjValue();
		cout << "Solution subproblem : " << currentObjSP << endl;

		if(1 - currentObjSP > -0.00001){
			cout << "Convergence is reached" << endl;
		}
		else{
			for(int i = 0; i < data.numberOfItems; i++){
				solSP.push_back(round(solverSP.getValue(x[i])));
			}
		}
		solverMaster.end();
		solverSP.end();
		launchBP(m, data, solution, Ff, maxDeviation, lambda, objMaster, constraint, modelMaster, level, sp, x, delta, mu, objSP, modelSP, solSP, a, b, countColumns);
	}
	else{

		// item is a double -> launch separation







		exit(-1);
	}



	exit(-1);
}


void branch_and_price(MasterGCol* m, Data const &data, Solution solution, First_Fit &Ff, int const &maxDeviation, IloNumVarArray lambda, IloExpr objMaster,IloRangeArray constraint, IloModel modelMaster, int level){

	SubProblemGCol* sp = new SubProblem();
	IloNumVarArray x(sp->env, data.numberOfItems, 0, 1, ILOINT);
	IloNumVarArray delta(sp->env, data.numberOfItems, 0, IloInfinity);
	IloNumVar mu(sp->env, 0, IloInfinity);
	IloExpr objSP(sp->env);
	IloModel modelSP(sp->env);
	vector<int> solSP;

	// solve subproblem
	solSP = sp->solveSubProblemNotSameDev(data, m->pi, maxDeviation, x, delta, mu, objSP, modelSP);
	cout << "Solution subproblem : " << sp->optSP << endl;

	IloNumArray a(m->env);
	IloNumArray b(m->env);
	a.add(0);
	b.add(IloInfinity);
	double solMaster = 0;
	int countColumns = Ff.nb_bins;

	launchBP(m, data, solution, Ff, maxDeviation, lambda, objMaster, constraint, modelMaster, level, sp, x, delta, mu, objSP, modelSP, solSP, a, b, countColumns);

}




// exchange between MasterNotSameDev and SpNotSameDev
Solution BranchPrice::solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit inutile){
	Solution solution; // empty solution

	// master initialization
	MasterGCol* m = new Master();
	IloNumVarArray lambda(m->env, Ff.nb_bins, 0, IloInfinity);
	IloExpr objMaster(m->env);
	IloRangeArray constraint(m->env);
	IloModel modelMaster(m->env);
	solution = m->solveMasterNotSameDev(data, Ff.sol_init, maxDeviation, lambda, objMaster, constraint, modelMaster); // solve the master
	cout << "Solution master : " << m->solMaster << endl;

	branch_and_price(m, data, solution, Ff, maxDeviation, lambda, objMaster, constraint, modelMaster, 0);


	return solution;
}

































