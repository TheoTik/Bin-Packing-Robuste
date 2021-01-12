/*
 * First_Fit.cpp
 *
 *  Created on: 23 oct. 2020
 *      Author: AMETANA & COMBALBERT
 */

#include "First_Fit.hpp"
#include <iostream>


// constructeur par défaut
First_Fit::First_Fit(){
	nb_bins=0;
}

/**** aperçu des configurations***/
void First_Fit::affiche(){

	int nb_item=sol_init[0].size();

	cout <<"NOMBREE DE BINS : "<<nb_bins<<endl;

	for(int i=0;i<nb_bins;i++){
		cout<<"**************  Bin "<<i+1<<"  *************** "<<endl;
		for(int j=0;j<nb_item;j++){
			cout<<"Nombre item "<<j+1<<": ";
			cout<<sol_init[i][j]<<endl;
		}
	}
}


/*** ajout d'une colonne ou de configuration realisable***/
void First_Fit::Ajout_colonne(const vector<int> &col){

	int n=col.size();

	nb_bins +=1;
	sol_init.resize(nb_bins);
	sol_init[nb_bins-1].resize(n);

	for(int i=0;i<n;i++){
		sol_init[nb_bins-1][i]=col[i];
	}
}

void transpose(vector<vector<int> > &b)
{
	if (b.size() == 0)
		return;

	vector<vector<int> > trans_vec(b[0].size(), vector<int>());

	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < b[i].size(); j++)
		{
			trans_vec[j].push_back(b[i][j]);
		}
	}

	b = trans_vec;    // <--- reassign here
}

/*****création des configurations initiales***/

void First_Fit::launch(const Data &instance){

	int n= instance.numberOfItems;
	int W= instance.sizeBin;
	vector<int> A(n);

	for(int i=0;i<n;i++){
		A[i]=i;
	}

	// trie decroissant des items en fonction de la taille
	for(int i=0;i<n;i++){
		int max= instance.sizeObj[i];
		int pos=i;
		for(int j=i+1;j<n;j++){
			if(instance.sizeObj[j]>max){
				pos=j;
			}
		}
		int inter=A[pos];
		A[pos]=A[i];
		A[pos]=inter;
	}

	//construction de la solution
	sol_init.resize(n);
	vector<int> bins_open(1);
	bins_open[0]=W;
	sol_init.resize(1);
	sol_init[0].resize(n);
	for(int i=0;i<n;i++){
		sol_init[0][i]=0;
	}

	for(int i=0;i<n;i++){
		int item=A[i];
		int taille= instance.sizeObj[item] + instance.deviationObj[item];

		int k=0;

		while((taille>bins_open[k])){
			k +=1;
			int nb=bins_open.size();
			if(k>=nb){
				bins_open.push_back(W);
				sol_init.resize(nb+1);
				sol_init[nb].resize(n);
				for(int l=0;l<n;l++){
					sol_init[nb][l]=0;
				}
			}

		}
		bins_open[k] -= taille;
		sol_init[k][item] +=1;
	}
	nb_bins=bins_open.size();

	transpose(sol_init);
}


void First_Fit::launchSameDev(const int &numeroMethod, const Data &instance, const int &maxDeviation){


	if(numeroMethod == 0){
		int n= instance.numberOfItems;
		int W= instance.sizeBin;
		vector<int> A(n);

		for(int i=0;i<n;i++){
			A[i]=i;
		}

		int countDeviation = 0;
		int max = 0;
		int ifwecanaddaDeviation = 0;

		// trie decroissant des items en fonction de la taille
		for(int i=0;i<n;i++){

			if(countDeviation > maxDeviation){
				ifwecanaddaDeviation = 0;
			}
			else{
				ifwecanaddaDeviation = instance.deviation;
				countDeviation++;
			}

			max= instance.sizeObj[i] + ifwecanaddaDeviation;

			int pos=i;
			for(int j=i+1;j<n;j++){
				if(instance.sizeObj[j]+ifwecanaddaDeviation>max){
					max=instance.sizeObj[j]+ifwecanaddaDeviation;
					pos=j;
				}
			}
			int inter=A[pos];
			A[pos]=A[i];
			A[pos]=inter;
		}

		//construction de la solution
		sol_init.resize(n);
		vector<int> bins_open(1);
		bins_open[0]=W;
		sol_init.resize(1);
		sol_init[0].resize(n);
		for(int i=0;i<n;i++){
			sol_init[0][i]=0;
		}

		int taille = 0;
		countDeviation = 1;
		for(int i=0;i<n;i++){

			int item=A[i];

			if(countDeviation > maxDeviation){
				ifwecanaddaDeviation = 0;
			}
			else{
				ifwecanaddaDeviation = instance.deviation;
				countDeviation++;
			}
			taille= instance.sizeObj[item] + ifwecanaddaDeviation;

			int k=0;

			while((taille>bins_open[k])){
				k +=1;
				int nb=bins_open.size();
				if(k>=nb){
					bins_open.push_back(W);
					sol_init.resize(nb+1);
					sol_init[nb].resize(n);
					for(int l=0;l<n;l++){
						sol_init[nb][l]=0;
					}
				}

			}
			bins_open[k] -= taille;
			sol_init[k][item] +=1;
		}
		nb_bins=bins_open.size();

		transpose(sol_init);

	}

	else if(numeroMethod==1){
		int n= instance.numberOfItems;
		int W= instance.sizeBin - maxDeviation*instance.deviation;
		vector<int> A(n);

		for(int i=0;i<n;i++){
			A[i]=i;
		}

		// trie decroissant des items en fonction de la taille
		for(int i=0;i<n;i++){
			int max= instance.sizeObj[i];
			int pos=i;
			for(int j=i+1;j<n;j++){
				if(instance.sizeObj[j]>max){
					max=instance.sizeObj[j];
					pos=j;
				}
			}
			int inter=A[pos];
			A[pos]=A[i];
			A[pos]=inter;
		}

		//construction de la solution
		sol_init.resize(n);
		vector<int> bins_open(1);
		bins_open[0]=W;
		sol_init.resize(1);
		sol_init[0].resize(n);
		for(int i=0;i<n;i++){
			sol_init[0][i]=0;
		}

		for(int i=0;i<n;i++){
			int item=A[i];
			int taille= instance.sizeObj[item];

			int k=0;

			while((taille>bins_open[k])){
				k +=1;
				int nb=bins_open.size();
				if(k>=nb){
					bins_open.push_back(W);
					sol_init.resize(nb+1);
					sol_init[nb].resize(n);
					for(int l=0;l<n;l++){
						sol_init[nb][l]=0;
					}
				}

			}
			bins_open[k] -= taille;
			sol_init[k][item] +=1;
		}
		nb_bins=bins_open.size();

		transpose(sol_init);
	}

}















