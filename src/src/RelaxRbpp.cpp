#include <iostream>
#include <ilcplex/ilocplex.h>

#include "Algorithm.hpp"


// ----------------------------------------------------------------------------
//  La méthode solve résout la formulation classique de l'ULS
// ----------------------------------------------------------------------------
Solution RelaxRbpp::solve(Data const &data, int const &maxDeviations, First_Fit Ff, First_Fit inutile) {


	// upper bound on the number of bins (to change)
	int U = data.numberOfItems;


	Solution solution; // empty solution
	// --- Creation of the CPLEX environment ---
	IloEnv env;
	try {
		// --- Creation of the model associated with the environment ---
		IloModel model(env);
		int n = data.numberOfItems;
		int W = data.sizeBin;

		// --- Creation of the variables ---

		// primales variables
		IloArray<IloNumVarArray> x(env, n);
		for(int i = 0; i < n; i++){
			x[i] = IloNumVarArray(env, U, 0, 1);
		}
		IloNumVarArray y(env, U, 0, 1);

		// duales variables
		IloArray<IloNumVarArray> pi(env, n);
		for(int i = 0; i < n; i++){
			pi[i] = IloNumVarArray(env, U, 0, IloInfinity);
		}
		IloNumVarArray mu(env, U, 0, IloInfinity);

		// --- Creation of the objective function ---
		IloExpr obj(env);
		for(int i = 0; i < U; i++){
			obj += y[i];
		}
		model.add(IloMinimize(env, obj));
		obj.end();

		// --- Creation of the constraints ---
		IloExpr cstr(env);

		// assignement constraint
		for(int i = 0; i < n; i++){
			for(int j = 0; j < U; j++){
				cstr += x[i][j];
			}
			model.add(cstr == 1);
			cstr.clear();
		}

		// knapsack constraint
		for(int j = 0; j < U; j++){
			for(int i = 0; i < n; i++){
				cstr += data.sizeObj[i]*x[i][j] + pi[i][j];
			}
			cstr += mu[j]*maxDeviations;
			model.add(cstr <= W*y[j]);
			cstr.clear();
		}

		// dual constraint
		for(int j = 0; j < U; j++){
			for(int i = 0; i < n; i++){
				cstr = mu[j] + pi[i][j] - data.deviationObj[i]*x[i][j];
				model.add(cstr >= 0);
			}
		}
		cstr.end();

		// --- Solver configuration ---
		IloCplex solver(model); //< creates the solver for the model
		solver.setParam(IloCplex::Param::Threads, 1);  //< Limitation du nombre de threads à 1
		solver.setParam(IloCplex::Param::TimeLimit, 300); //< Temps limite (en secondes)

		//solver.exportModel("model.lp"); //< écrit le modèle au format LP (optionnel)

		// --- Solver launch ---
		solver.solve();

		// --- Solver results retrieval ---
		cout <<"\n----------------------------------"<<endl;
		bool existsSol = false;
		if (solver.getStatus() == IloAlgorithm::Optimal){
			cout<<"Status de la solution: OPTIMAL"<<endl;
			existsSol = true;
		}
		if (solver.getStatus() == IloAlgorithm::Feasible){
			cout<<"Status de la solution: TEMPS LIMITE"<<endl;
			existsSol = true;
		}
		if (solver.getStatus() == IloAlgorithm::Infeasible || solver.getStatus() == IloAlgorithm::InfeasibleOrUnbounded){
			cout<<"Status de la solution: IRREALISABLE"<<endl;
		}
		if (solver.getStatus() == IloAlgorithm::Unbounded){
			cout<<"Status de la solution: NON BORNE"<<endl;
		}
		if (existsSol)
			cout<<"Valeur de la fonction objectif : "<< (solver.getObjValue())<<endl; // ne pas oublier d'arrondir si le coût doit être entier
		cout<<"CPU Time (s) : "<<solver.getTime()<<endl;
		cout<<"----------------------------------\n"<<endl;

		if (existsSol) {

			// build solution


		}

	} catch (IloException &e) {
		cerr << "Concert Technology exception" << e << endl;
		exit(EXIT_FAILURE);
	}

	env.end(); //< frees the environment
	return solution;
}
