#include <iostream>
#include <ilcplex/ilocplex.h>
#include "Master.hpp"



Solution Master::solveMasterNotSameDev(Data const &data, vector<vector<int>> const &P, int const &maxDeviations, IloNumVarArray lambda, IloExpr obj, IloRangeArray constraint, IloModel modele) {

	Solution solution; // empty solution

	try {

		// --- Creation of the objective function ---
		for(int i = 0; i < P[0].size(); i++){
			obj += lambda[i];
		}
		objective = IloMinimize(env, obj);
		modele.add(objective);
		obj.end();


		// --- Creation of the constraints ---
		IloExpr c(env);
		IloRangeArray cstr1(env);

		// assignement constraint
		for(int i = 0; i < P.size(); i++){
			for(int j = 0; j < P[0].size(); j++){
				c += P[i][j]*lambda[j];
			}
			cstr1.add(c == 1);
			c.clear();
		}
		c.end();
		modele.add(cstr1);

		// --- Solver configuration ---
		IloCplex solver(modele); //< creates the solver for the model
		solver.setParam(IloCplex::Param::Threads, 1);  //< Limitation du nombre de threads à 1
		solver.setParam(IloCplex::Param::TimeLimit, 300); //< Temps limite (en secondes)

		//solver.exportModel("model.lp"); //< écrit le modèle au format LP (optionnel)

		// --- Solver launch ---
		solver.solve();

		// get dual values
		IloNumArray v(env);
		cout << "Duals : ";
		solver.getDuals(v, cstr1);
		cout << v << endl;
		for(int i = 0; i < v.getSize(); i++){
			pi.push_back(v[i]);
		}

		solMaster = solver.getObjValue();
		constraint.add(cstr1);


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

		if(existsSol){
			// build solution
			vector<double> a;
			for(int i = 0; i < P.size(); i++){
				solution.bins.push_back(a);
				for(int j = 0; j < P[0].size(); j++){
					solution.bins[i].push_back(P[i][j]* solver.getValue(lambda[j]));
				}
			}

		}


	} catch (IloException &e) {
		cerr << "Concert Technology exception" << e << endl;
		exit(EXIT_FAILURE);
	}

	return solution;

}



Solution Master::solveMasterSameDev(Data const &data, vector<vector<int>> const &P, vector<vector<int>> const &P1, int const &maxDeviation, IloNumVarArray lambdap, IloNumVarArray lambdaq, IloExpr obj, IloRangeArray constraint, IloModel modele){

	Solution solution;

	try {

		// --- Creation of the objective function ---
		for(int i = 0; i < P[0].size(); i++){
			obj += lambdap[i];
		}
		for(int i = 0; i < P1[0].size(); i++){
			obj += lambdaq[i];
		}
		objective = IloMinimize(env, obj);
		modele.add(objective);
		obj.end();


		// --- Creation of the constraints ---
		IloExpr c(env);
		IloRangeArray cstr1(env);

		// assignement constraint
		for(int i = 0; i < data.numberOfItems; i++){
			for(int j = 0; j < P[0].size(); j++){
				c += P[i][j]*lambdap[j];
			}
			for(int j = 0; j < P1[0].size(); j++){
				c += P1[i][j]*lambdaq[j];
			}
			cstr1.add(c == 1);
			c.clear();
		}
		c.end();
		modele.add(cstr1);

		// --- Solver configuration ---
		IloCplex solver(modele); //< creates the solver for the model
		solver.setParam(IloCplex::Param::Threads, 1);  //< Limitation du nombre de threads à 1
		solver.setParam(IloCplex::Param::TimeLimit, 300); //< Temps limite (en secondes)

		//solver.exportModel("model.lp"); //< écrit le modèle au format LP (optionnel)

		// --- Solver launch ---
		solver.solve();

		// get dual values
		IloNumArray v(env);
		cout << "Duals : ";
		solver.getDuals(v, cstr1);
		cout << v << endl;
		for(int i = 0; i < v.getSize(); i++){
			pi.push_back(v[i]);
		}

		solMaster = solver.getObjValue();
		constraint.add(cstr1);


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

		if(existsSol){
			// build solution
			vector<double> a;
			for(int i = 0; i < P.size(); i++){
				solution.bins.push_back(a);
				for(int j = 0; j < P[0].size(); j++){
					solution.bins[i].push_back(P[i][j]* solver.getValue(lambdap[j]));
				}
			}
			for(int i = 0; i < P1.size(); i++){
				solution.bins.push_back(a);
				for(int j = 0; j < P1[0].size(); j++){
					solution.bins[i].push_back(P1[i][j]* solver.getValue(lambdaq[j]));
				}
			}

		}


	} catch (IloException &e) {
		cerr << "Concert Technology exception" << e << endl;
		exit(EXIT_FAILURE);
	}

	return solution;
}







