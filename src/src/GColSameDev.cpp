#include <iostream>
#include <ilcplex/ilocplex.h>
#include <cmath>;

#include "Algorithm.hpp"
#include "Master.hpp"
#include "SubProblem.hpp"

/*

void affiche(vector<vector<int>> P){
	for(int i = 0; i < P.size(); i++){
		for(int j = 0; j < P[0].size(); j++){
			cout << P[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void affiche(vector<int> gamma){
	for(int i = 0; i < gamma.size(); i++){
		cout << gamma[i] << " ";
	}
	cout << endl;
}

void affiche(IloModel model){
	cout << model << endl;
}*/

vector<vector<double>> transpo(vector<vector<double> > &b)
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


// exchange between MasterSameDev and SpSameDev
Solution GColSameDev::solve(Data const &data, int const &maxDeviation, First_Fit Ffp, First_Fit Ffq){

	Solution solution; // empty solution

	// master initialization
	MasterGCol* m = new Master(); 	// MasterSameDev a renommé
	IloNumVarArray lambdap(m->env, Ffp.nb_bins, 0, IloInfinity);
	IloNumVarArray lambdaq(m->env, Ffq.nb_bins, 0, IloInfinity);
	IloExpr objMaster(m->env);
	IloRangeArray constraint(m->env);
	IloModel modelMaster(m->env);
	solution = m->solveMasterSameDev(data, Ffp.sol_init, Ffq.sol_init, maxDeviation, lambdap, lambdaq, objMaster, constraint, modelMaster);
	cout << "Solution master : " << m->solMaster << endl;

	// subproblems initialization
	SubProblemGCol* sp0 = new SubProblem();
	SubProblemGCol* sp1 = new SubProblem();

	// initialize solution of the subproblems
	vector<int> solSP1;
	vector<int> solSP0;

	// we have 2 types of bins, we have to solve the 2 subproblems and add associated constraint
	IloNumVarArray xSP0(sp0->env, data.numberOfItems, 0, 1, ILOINT);
	IloNumVarArray xSP1(sp1->env, data.numberOfItems, 0, 1, ILOINT);

	// creation of the objective subproblems functions
	IloExpr objSP0(sp0->env);
	IloExpr objSP1(sp1->env);

	// model of the subproblems
	IloModel modelSP0(sp0->env);
	IloModel modelSP1(sp1->env);

	// solve the subproblems
	if(maxDeviation != 0){
		solSP0 = sp0->solveSubProblemSameDev(data, m->pi, maxDeviation, 0, xSP0, objSP0, modelSP0);
		cout << "Solution subproblem SP0 : " << sp0->optSP << endl;

	}
	solSP1 = sp1->solveSubProblemSameDev(data, m->pi, maxDeviation, 1, xSP1, objSP1, modelSP1);
	cout << "Solution subproblem SP1 : " << sp1->optSP << endl;

	//affiche(solSP0);
	//affiche(solSP1);

	// Columns generation
	bool convergence = false;
	double currentObjSP0 = 0;
	double currentObjSP1 = 0;
	int cpt = 0;
	IloNumArray a(m->env);
	IloNumArray b(m->env);
	a.add(0);
	b.add(IloInfinity);
	double solMaster = 0;
	int countColumnsSP0 = Ffp.nb_bins;
	int countColumnsSP1 = Ffq.nb_bins;

	while(!convergence){

		if(maxDeviation != 0){
			// add a new variable for SP0 if Γ != 0
			IloNumColumnArray col(m->env,1);
			col[0] = m->objective(1);
			for(IloInt i = 0; i < constraint.getSize(); i++){
				col[0] += constraint[i](solSP0[i]);
			}
			solSP0.clear();

			IloNumVarArray lambdavar(m->env, col, a, b);
			lambdap.add(lambdavar);
			countColumnsSP0++;
		}

		// add a new variable for SP1
		IloNumColumnArray col(m->env,1);
		col[0] = m->objective(1);
		for(IloInt i = 0; i < constraint.getSize(); i++){
			col[0] += constraint[i](solSP1[i]);
		}
		solSP1.clear();

		IloNumVarArray lambdavar(m->env, col, a, b);
		lambdaq.add(lambdavar);
		countColumnsSP1++;

		// solve the master and keep duals in m->vals
		solution.bins.clear();
		IloCplex solverMaster(modelMaster);
		solverMaster.solve();
		IloNumArray vals(m->env);

		solverMaster.getDuals(vals, constraint);
		cout << "Duals : " << vals << endl;

		solMaster = solverMaster.getObjValue();
		cout << "Solution master : " << solMaster << endl;

		// solve subproblems and change pi
		if(maxDeviation != 0){
			objSP0.setLinearCoefs(xSP0, vals);
		}
		IloCplex solverSP0(modelSP0);
		if(maxDeviation != 0){
			solverSP0.solve();
			// check solution
			currentObjSP0 = solverSP0.getObjValue();
			cout << "Solution subproblem SP0 : " << currentObjSP0 << endl;
		}

		// same for the 2nd subproblem
		objSP1.setLinearCoefs(xSP1, vals);
		IloCplex solverSP1(modelSP1);
		solverSP1.solve();

		// check solution
		currentObjSP1 = solverSP1.getObjValue();
		cout << "Solution subproblem SP1 : " << currentObjSP1 << endl;


		if(1.0 - currentObjSP0 > -0.0001 && 1.0 - currentObjSP1 > -0.0001){ // check convergence
			convergence = true;
			cout << "Convergence is reached" << endl;

			// build solution
			vector<double> solLambdap;
			vector<double> solLambdaq;

			if(maxDeviation != 0){
				for(int i = 0; i < countColumnsSP0; i++){
					solLambdap.push_back(solverMaster.getValue(lambdap[i]));
				}
			}

			for(int i = 0; i < countColumnsSP1; i++){
				solLambdaq.push_back(solverMaster.getValue(lambdaq[i]));
			}


			vector<vector<double>> finalSol;
			vector<double> vect(constraint.getSize());
			for(int i = 0; i < solLambdap.size(); i++){
				finalSol.push_back(vect);
			}
			for(int i = 0; i < solLambdaq.size(); i++){
				finalSol.push_back(vect);
			}

			for(int i = 0; i < constraint.getSize(); i++){

				for (IloExpr::LinearIterator it = IloExpr(constraint[i].getExpr()).getLinearIterator(); it.ok();++it){

					int posVarInLambda = lambdaq.find(it.getVar()); // find the position in the list of variables

					if(posVarInLambda != -1){

						if(solLambdaq[posVarInLambda] != 0){
							finalSol[posVarInLambda][i] = solLambdaq[posVarInLambda];
						}
					}
					else{
						posVarInLambda = lambdap.find(it.getVar());

						if(maxDeviation != 0 && solLambdap[posVarInLambda] != 0){
							finalSol[posVarInLambda][i] = solLambdap[posVarInLambda];
						}
					}
				}
			}

			vector<vector<double>> finalSol1 = transpo(finalSol);
			vector<vector<double>> finalSol2;

			for(int i = 0; i < constraint.getSize(); i++){
				if(!finalSol1[i].empty()){
					finalSol2.push_back(finalSol1[i]);
				}
			}

			solution.bins = transpo(finalSol2);
			solution.opt = solMaster;

		}
		else{
			if(maxDeviation != 0){
				for(int i = 0; i < data.numberOfItems; i++){
					solSP0.push_back(round(solverSP0.getValue(xSP0[i])));
				}
			}

			for(int i = 0; i < data.numberOfItems; i++){
				solSP1.push_back(round(solverSP1.getValue(xSP1[i])));
			}

			//affiche(solSP0);
			//affiche(solSP1);
		}
		solverMaster.end();
		solverSP0.end();
		solverSP1.end();
	}

	/*for(int i = 0; i < solution.bins.size(); i++){
		for(int j = 0; j < solution.bins[0].size(); j++){
			cout << solution.bins[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;*/

	return solution;
}





