//
// Created by afroger001 on 05/11/2020.
//

#include "Data.hpp"

#include <fstream>
#include <iostream>


Data::Data(string const & filePath){
	numberOfItems = 0;
	sizeBin = 0;
	cout << "Read the data from the file : " << filePath << endl;
	ifstream f("instances_bp_robuste/" + filePath + ".txt");
	if (!f.good()) {
		cerr << "Error: Could not open file " << filePath << endl;
		exit(EXIT_FAILURE);
	}
	string name;
	f >> sizeBin;
	f >> numberOfItems;
	char line[256];
	f.getline(line, 256);
	int s, d;
	for (int t = 0; t < numberOfItems; ++t) {
		f >> s >> d;
		sizeObj.push_back(s);
		deviationObj.push_back(d);
	}
	f.close();

	deviation = deviationObj[0];
}

Solution::Solution(): bins() { };



