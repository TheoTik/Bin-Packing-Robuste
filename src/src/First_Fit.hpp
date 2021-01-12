/*
 * First_Fit.hpp
 *
 *  Created on: 23 oct. 2020
 *      Author: AMETANA & COMBALBERT
 */

#ifndef SRC_FIRST_FIT_HPP_
#define SRC_FIRST_FIT_HPP_

#include"Data.hpp"


/* la classe First_Fit permet de stocker les configurations realisables qui seront passées en parametre au probleme maitre
 *  en mettant les demandes a 1 et grace a un algorithme de First fit on trouve un sous ensemble de configuration realisable
 *  qui  une fois en parametre au maitre lui permet de trouver une solution realisable.
 *  La classe permet egalement d'ajouter la configuration realisable qui sera donnee par le sous probleme.
 */

class First_Fit{
public:
	int nb_bins;
	vector<vector<int>> sol_init;
public:
	First_Fit();
	void launch(const Data &instance);
	void launchSameDev(const int &a, const Data &instance, const int &maxDeviation);
	void affiche();
	void Ajout_colonne(const vector<int> &col);


};



#endif /* SRC_FIRST_FIT_HPP_ */
