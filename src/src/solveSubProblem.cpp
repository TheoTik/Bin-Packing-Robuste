#include <iostream>
#include <ilcplex/ilocplex.h>
#include "SubProblem.hpp"



vector<int> SubProblem::solveSubProblemSameDev(Data const &data, vector<double> const &pi, int const &maxDeviations, int const &numeroSP, IloNumVarArray x, IloExpr obj, IloModel modele){

	vector<int> sol;

	try {
		// --- Creation of the model associated with the environment ---
		//int n = data.numberOfItems;
		int W = data.sizeBin;

		// --- Creation of the objective function ---
		for(int i = 0; i < pi.size(); i++){
			obj += x[i]*pi[i];
		}
		objective = IloMaximize(env, obj);
		modele.add(objective);


		// --- Creation of the constraints ---
		IloExpr cstr(env);

		// knapsack constraint
		if(numeroSP == 0){ // 1st subproblem
			for(int i = 0; i < pi.size(); i++){
				cstr += (data.sizeObj[i]+data.deviation)*x[i];
			}
			modele.add(cstr <= W);
			cstr.clear();
			for(int i = 0; i < pi.size(); i++){
				cstr += x[i];
			}
			modele.add(cstr <= maxDeviations);
		}
		else if(numeroSP == 1){ // 2nd subproblem
			for(int i = 0; i < pi.size(); i++){
				cstr += data.sizeObj[i]*x[i];
			}
			modele.add(cstr <= (W-maxDeviations*data.deviation));
		}
		cstr.end();

		// --- Solver configuration ---
		IloCplex solver(modele); //< creates the solver for the model
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
			cout<<"Valeur de la fonction objectif : "<< (solver.getObjValue())<<endl;
		cout<<"CPU Time (s) : "<<solver.getTime()<<endl;
		cout<<"----------------------------------\n"<<endl;

		if (existsSol) {

			// build solution
			for(int i = 0; i < pi.size(); i++){
				sol.push_back(solver.getValue(x[i]));
			}

			optSP = solver.getObjValue();
		}

	} catch (IloException &e) {
		cerr << "Concert Technology exception" << e << endl;
		exit(EXIT_FAILURE);
	}

	return sol;
}


vector<int> SubProblem::solveSubProblemNotSameDev(Data const &data, vector<double> const &pi, int const &maxDeviations, IloNumVarArray x, IloNumVarArray delta, IloNumVar mu, IloExpr obj, IloModel modelSP) {
	vector<int> sol;


	try {
		// --- Creation of the model associated with the environment ---
		//int n = data.numberOfItems;
		int W = data.sizeBin;

		// --- Creation of the objective function ---
		for(int i = 0; i < pi.size(); i++){
			obj += x[i]*pi[i];
		}
		objective = IloMaximize(env, obj);
		modelSP.add(objective);


		// --- Creation of the constraints ---
		IloExpr cstr(env);

		// knapsack constraint
		for(int i = 0; i < pi.size(); i++){
			cstr += data.sizeObj[i]*x[i] + delta[i];
		}
		cstr += mu*maxDeviations;
		modelSP.add(cstr <= W);
		cstr.clear();


		for(int i = 0; i < pi.size(); i++){
			cstr += delta[i] - data.deviationObj[i]*x[i];
		}
		modelSP.add(cstr + mu >= 0);


		// --- Solver configuration ---
		IloCplex solver(modelSP); //< creates the solver for the model
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
			cout<<"Valeur de la fonction objectif : "<< (solver.getObjValue())<<endl;
		cout<<"CPU Time (s) : "<<solver.getTime()<<endl;
		cout<<"----------------------------------\n"<<endl;

		if (existsSol) {

			// build solution
			for(int i = 0; i < pi.size(); i++){
				sol.push_back(solver.getValue(x[i]));
			}

			optSP = solver.getObjValue();

		}

	} catch (IloException &e) {
		cerr << "Concert Technology exception" << e << endl;
		exit(EXIT_FAILURE);
	}

	return sol;
}
