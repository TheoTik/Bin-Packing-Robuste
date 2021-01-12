#include <iostream>
#include <ilcplex/ilocplex.h>
#include <cmath>

#include "Algorithm.hpp"
#include "Master.hpp"
#include "SubProblem.hpp"


vector<vector<double>> transp(vector<vector<double> > &b)
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

vector<vector<double>> keepSol(IloRangeArray constraint, vector<double> solution, IloNumVarArray lambda){
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

	vector<vector<double>> finalSol1 = transp(finalSol);
	vector<vector<double>> finalSol2;

	for(int i = 0; i < constraint.getSize(); i++){
		if(!finalSol1[i].empty()){
			finalSol2.push_back(finalSol1[i]);
		}
	}

	return transp(finalSol2);
}



ILOBRANCHCALLBACK7(callback, Data const&, data, int, item, IloNumVarArray, lambda, double, solMaster, vector<double>, solLambda, vector<double>, coefObjMaster, double, maxobj)
{

	if (getBranchType() != BranchOnVariable)
		return;

	/*IntegerFeasibilityArray feas;
	getFeasibilities(feas, lambda);

	double maxinf = 0.0;

	// we choose the clother than 0 or 1
	for(int i = 0; i < lambda.getSize(); i++){
		if (feas[i] == Infeasible) { // it exists one item not integer
			double inf = solLambda[i] - floor(solLambda[i]);
			if(inf > 0.5){
				inf = 1.0 - inf;
			}
			if(inf >= maxinf && (inf > maxinf || abs(coefObjMaster[i]) >= maxobj)){
				item = i;
				maxinf = inf;
				maxobj = abs(coefObjMaster[i]);
			}
		}
	}*/

	if(item >= 0){
		//makeBranch(lambda[item], solLambda[item], IloCplex::BranchUp, solMaster);
		//makeBranch(lambda[item], solLambda[item], IloCplex::BranchDown, solMaster);

		makeBranch(lambda[item], 0, IloCplex::BranchUp, solMaster);
	}

	return;
}


ILOUSERCUTCALLBACK7(usercut, int, b, int, item, IloNumVarArray, lambda, double, solMaster, vector<double>, solLambda, vector<double>, coefObjMaster, double, maxobj){

	if(b == 1){
		b=0;
		abortCutLoop();
	}

	if (!isAfterCutLoop())
		return;

	if(item != -1){
		add(lambda[item] <= 0);
	}

}

void afficher(vector<vector<double>> P){
	for(int i = 0; i < P.size(); i++){
		for(int j = 0; j < P[0].size(); j++){
			cout << P[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

// exchange between MasterNotSameDev and SpNotSameDev
Solution GColNotSameDev::solve(Data const &data, int const &maxDeviation, First_Fit Ff, First_Fit inutile){
	Solution solution; // empty solution
	Solution solutionfinal;

	// number of bins during the test
	int currentMin = 0;
	int bestMin = 0;

	// master initialization
	MasterGCol* m = new Master();
	IloNumVarArray lambda(m->env, Ff.nb_bins, 0, IloInfinity);
	IloExpr objMaster(m->env);
	IloRangeArray constraint(m->env);
	IloModel modelMaster(m->env);
	solution = m->solveMasterNotSameDev(data, Ff.sol_init, maxDeviation, lambda, objMaster, constraint, modelMaster); // solve the master
	cout << "Solution master : " << m->solMaster << endl;

	int tailleConstraint = constraint.getSize();


	// subproblem initialization
	SubProblemGCol* sp = new SubProblem();
	IloNumVarArray x(sp->env, data.numberOfItems, 0, 1, ILOINT);
	IloNumVarArray delta(sp->env, data.numberOfItems, 0, IloInfinity);
	IloNumVar mu(sp->env, 0, IloInfinity);
	IloExpr objSP(sp->env);
	IloModel modelSP(sp->env);
	vector<int> solSP = sp->solveSubProblemNotSameDev(data, m->pi, maxDeviation, x, delta, mu, objSP, modelSP);
	cout << "Solution subproblem : " << sp->optSP << endl;


	// Columns generation
	bool convergence = false;

	double currentObjSP = 0;
	int cpt = 0;
	IloNumArray a(m->env);
	IloNumArray b(m->env);
	a.add(0);
	b.add(IloInfinity);
	IloNumArray a1(sp->env);
	IloNumArray b1(sp->env);
	a1.add(0);
	b1.add(1);
	double solMaster = 0;
	int countColumns = Ff.nb_bins;
	vector<double> solLambda;
	vector<int> itemiChosenForBranch;
	vector<int> itemjChosenForBranch;

	IloExprArray previous(m->env);

	bool endALL = false;

	while(!endALL){

		while(!convergence){

			// add a new variable (does not work with IloNumColumn) we use an Array of size 1
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
			IloNumArray vals(m->env);


			solverMaster.getDuals(vals, constraint);
			cout << "Duals : " << vals << endl;

			solMaster = solverMaster.getObjValue();
			cout << "Solution master : " << solMaster << endl;

			// solve subproblem and modify pi
			objSP.setLinearCoefs(x, vals);
			IloCplex solverSP(modelSP);
			solverSP.solve();

			// check solution
			currentObjSP = solverSP.getObjValue();
			cout << "Solution subproblem : " << currentObjSP << endl;

			if(1.0 - currentObjSP > -0.00001){ // check convergence
				convergence = true;
				cout << "Convergence is reached" << endl;


				// build solution
				for(int i = 0; i < countColumns; i++){
					solLambda.push_back(solverMaster.getValue(lambda[i]));
				}


				// a partir de la liste des configurations récupérer les solutions
				solution.bins = keepSol(constraint, solLambda, lambda);
				solution.opt = solMaster;

			}
			else{

				for(int i = 0; i < data.numberOfItems; i++){
					solSP.push_back(round(solverSP.getValue(x[i])));
				}

			}
			solverMaster.end();
			solverSP.end();
		}


		//afficher(solution.bins);
		//cout << solution.opt<< endl;


		///////// INTEGER //////// Problem with impact on subproblem objetive

		if(numeroMethod == 1){ // integer method -> add a branch in the master

			// solve Master with new branch

			// find an item which is a double
			int itemi = -1; // position of the item
			int itemj = -1; // position of the item
			double mininf = 1;
			double inf = 0.0;

			// we choose the clother to 0
			for(int i = 0; i < solution.bins.size(); i++){
				for(int j = 0; j < solution.bins[i].size(); j++){
					inf = solution.bins[i][j];

					if (inf > 0 && inf < 1) { // it exists one item not integer
						if(inf < mininf){
							itemi = i;
							itemj = j;
							mininf = inf;
						}
					}
				}
			}


			// if no item double than end !
			if(itemi != -1){

				// first branch
				IloExpr c(m->env);
				c += lambda[itemj]*solution.bins[itemi][itemj];

				// add first branch to the model
				if(c.getLinearIterator().ok()){
					constraint.add(c <= 0);
					modelMaster.add(constraint[constraint.getSize()-1]);
				}
				c.end();


				// solve master
				IloCplex solverMaster(modelMaster);
				solverMaster.solve();

				IloNumArray vals(m->env); // set of duals

				// vector which contains all the previous items chosen
				itemiChosenForBranch.push_back(itemi);
				itemjChosenForBranch.push_back(itemj);


				// if master is infeasible
				if(solverMaster.getStatus() == IloAlgorithm::Infeasible){
					bool nofinish = false;

					// return to the previous node where the master is feasible
					while(!nofinish){

						// return to the previous node by deleting previous constraint and adding the other branch
						modelMaster.remove(constraint[constraint.getSize()-1]);

						constraint[constraint.getSize()-1] = constraint[constraint.getSize()-1].getExpr() >= 1;
						modelMaster.add(constraint[constraint.getSize()-1]);


						// we put back the impact on the subproblem without the previous branch
						if(previous.getSize() != 0){
							objSP -= previous[previous.getSize()-1];

							// we can delete the previous impact on the subproblem
							previous.setSize(previous.getSize()-1);
						}

						// solve the other branch
						IloCplex solverMaster2(modelMaster);
						solverMaster2.solve();

						if(solverMaster2.getStatus( )== IloAlgorithm::Infeasible){

							// node infeasible !!
							modelMaster.remove(constraint[constraint.getSize()-1]); // remove the last branch
							itemiChosenForBranch.pop_back(); // remove the last node

							itemjChosenForBranch.pop_back(); // remove the last node
							constraint.setSize(constraint.getSize()-1); // remove the last constraint

						}
						else{ // master feasible we can continue to branch on this node
							nofinish = true;

							// keep dual of the master
							solverMaster2.getDuals(vals, constraint);
							cout << "Duals : " << vals << endl;

							solMaster = solverMaster2.getObjValue();
							cout << "Solution master : " << solMaster << endl;

						}
					}

					// impact on the subproblem
					if(previous.getSize() != 0){
						previous.add(vals[vals.getSize()-1]*x[itemjChosenForBranch[itemjChosenForBranch.size()-1]]);
						objSP -= previous[previous.getSize()-1];
					}

				}
				else{
					// keep dual of the master
					solverMaster.getDuals(vals, constraint);
					cout << "Duals : " << vals << endl;

					solMaster = solverMaster.getObjValue();
					cout << "Solution master : " << solMaster << endl;

					// impact on the subproblem
					if(previous.getSize() != 0){
						previous.add(vals[vals.getSize()-1]*x[itemjChosenForBranch[itemjChosenForBranch.size()-1]]);
						objSP += previous[previous.getSize()-1];
					}
				}

				// solve subproblem and change pi
				objSP.setLinearCoefs(x, vals);
				IloCplex solverSP(modelSP);
				solverSP.solve();

				// check solution
				currentObjSP = solverSP.getObjValue();
				cout << "Solution subproblem : " << currentObjSP << endl;

				// solution of the subproblem
				for(int i = 0; i < data.numberOfItems; i++){
					solSP.push_back(round(solverSP.getValue(x[i])));
				}

			}
			else{	// continue to check the optimal

				currentMin = m->solMaster;

				if(currentMin > bestMin){
					bestMin = currentMin;
					solution.bins = keepSol(constraint, solLambda, lambda);
					cpt = 0;
				}
				else{
					cout << "fail " << endl;
					cout << "best relax : " << bestMin << endl;
					cpt++;
				}

				// does not evolve
				if(cpt == 1000){
					exit(-1);
				}

				cout << bestMin << endl;

				// continue branch and price

				// break the last variables
				// node complete !! remove until constraint = 0
				bool infea = false;
				while(!infea){
					bool nofinish = false;
					while(!nofinish){

						if(constraint.getSize() == tailleConstraint){
							cout << "Error !" <<endl;
							cout << "best solution found : " << bestMin << endl;
							exit(-1);
						}


						// check when we can go to the right
						if(constraint[constraint.getSize()-1].getLB() == 1){

							modelMaster.remove(constraint[constraint.getSize()-1]); // remove the last branch
							itemiChosenForBranch.pop_back(); // remove the last node

							itemjChosenForBranch.pop_back(); // remove the last node
							constraint.setSize(constraint.getSize()-1); // remove the last constraint

							if(previous.getSize() != 0){
								objSP += previous[previous.getSize()-1]; // put back impact on the subproblem
								previous.setSize(previous.getSize()-1); // delete node
							}


						}
						else{
							nofinish = true;
						}
					}


					// delete last constraint and add new branch
					modelMaster.remove(constraint[constraint.getSize()-1]);

					IloExpr expr = constraint[constraint.getSize()-1].getExpr();
					constraint[constraint.getSize()-1] = expr >= 1;

					modelMaster.add(constraint[constraint.getSize()-1]);

					if(previous.getSize() != 0){
						objSP -= previous[previous.getSize()-1]; // put back impact on the subproblem
						previous.setSize(previous.getSize()-1);
					}

					IloCplex solverModel(modelMaster);
					solverModel.solve();


					if(solverModel.getStatus() != IloAlgorithm::Infeasible){
						infea = true;
					}
					else{


						bool nofinish = false;

						IloNumArray vals(m->env);


						// return to the previous node where the master is feasible
						while(!nofinish){
							// return to the previous node by deleting previous constraint and adding the other branch
							modelMaster.remove(constraint[constraint.getSize()-1]);

							constraint[constraint.getSize()-1] = constraint[constraint.getSize()-1].getExpr() >= 1;
							modelMaster.add(constraint[constraint.getSize()-1]);


							// we put back the impact on the subproblem without the previous branch
							if(previous.getSize() != 0){
								objSP -= previous[previous.getSize()-1];
								// we can delete the previous impact on the subproblem
								previous.setSize(previous.getSize()-1);
							}

							// solve the other branch
							IloCplex solverMaster2(modelMaster);
							solverMaster2.solve();

							if(solverMaster2.getStatus( )== IloAlgorithm::Infeasible){

								// node infeasible !!
								modelMaster.remove(constraint[constraint.getSize()-1]); // remove the last branch
								itemiChosenForBranch.pop_back(); // remove the last node

								itemjChosenForBranch.pop_back(); // remove the last node
								constraint.setSize(constraint.getSize()-1); // remove the last constraint

							}
							else{ // master feasible we can continue to branch on this node
								nofinish = true;

								// solLambda must become integer, we put the associated coef to 1
								solution.bins[itemiChosenForBranch[itemiChosenForBranch.size()-1]][itemjChosenForBranch[itemjChosenForBranch.size()-1]] = 1;

								// keep dual of the master
								solverMaster2.getDuals(vals, constraint);
								cout << "Duals : " << vals << endl;

								solMaster = solverMaster2.getObjValue();
								cout << "Solution master : " << solMaster << endl;

							}
						}

						// impact on the subproblem
						previous.add(vals[vals.getSize()-1]*x[itemjChosenForBranch[itemjChosenForBranch.size()-1]]);
						objSP -= previous[previous.getSize()-1];



					}
				}

				convergence = false;

			}
		}
		else{
			endALL = true;
		}
	}


	return solution;
}

