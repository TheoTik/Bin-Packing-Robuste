#include <iostream>
#include <chrono>
#include <iostream>
#include <sys/time.h>
#include <ctime>
#include "Algorithm.hpp"
#include "First_Fit.hpp"
#include "Data.hpp"

using namespace std;
using std::cout; using std::endl;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::system_clock;

// ----------------------------------------------------------------------------
// Méthode principale (chemin vers le fichier de données en paramètre)
// ----------------------------------------------------------------------------
int main() {
	int method = 5;
	string file = "";
	cout << "Choose method : ";
	//cin >> method;
	cout << "Choose file (ex : r100_100_1) : ";
	//cin >> file;
	file = "r100_100_1";
	int maxDeviations = 0;
	cout << "Choose number of deviations (ex : 0,1, ..., 10) : ";
	//cin >> maxDeviations;

	// READ THE PROGRAM ARGUMENTS
	Data data(file);
	First_Fit Ffp;
	First_Fit Ffq;

	Algorithm* algo;

	auto ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
	cout<<"\n--------------------------------------"<<endl;
	switch (method) {
	case 1: {
		algo = new Rbpp();
		cout<<"Methode -> Dualité de programmation linéaire"<<endl;
		break;
	}
	case 2: {
		algo = new RelaxRbpp();
		cout<<"Methode -> Relaxation linéaire dualité"<<endl;
		break;
	}
	case 3: {
		Ffp.launch(data);
		algo = new GColNotSameDev();
		cout<<"Methode -> Generation de colonne à deviation d_i quelconques"<<endl;
		break;
	}
	case 4: {
		Ffp.launchSameDev(0, data, maxDeviations); // 0 : methode 0
		Ffq.launchSameDev(1, data, maxDeviations);
		algo = new GColSameDev();
		cout<<"Methode -> Generation de colonne à deviation d fixe"<<endl;
		break;
	}
	case 5:{
		Ffp.launch(data);
		algo = new GColNotSameDev();
		algo->numeroMethod = 1;
		cout<<"Methode -> Branch & Price à deviation d_i quelconques" << endl;
		break;
	}
	default: {
		cout << "Incorrect method value -> "<<method<<endl;
		exit(EXIT_FAILURE);
	}
	}
	cout<<"--------------------------------------"<<endl;

	Solution solution = algo->solve(data, maxDeviations, Ffp, Ffq);

	cout << "Number of Bins : " << solution.opt << endl;

	auto ms2 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

	double time = ms2 - ms;

	cout << "Temps d'execution : "<< time << " ms" << endl;
	delete algo;

	/*string cont = "n";
    cout << "continuer ? (o/n) : ";
    cin >> cont;

    if(cont == "o"){
    	main();
    }*/

	return 0;
}

