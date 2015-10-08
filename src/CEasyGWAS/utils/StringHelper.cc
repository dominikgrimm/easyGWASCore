/*
*StringHelper class contains several methods to deal with string
*Max Planck Institute for Developmental Biology and MPI for Biological Cybernetics
*@author: Dominik Grimm
*@year: 2010
*/
#include "StringHelper.h"

using namespace std;

/*
*@Description: Split String by delimiter
*@param: String to split, Delimiter string
*@return: vector<string> with all substrings
*/
vector<string> StringHelper::split(string str,const string& delimiter) {
	vector<string> splitV;
	int l = 0;
	while(l!=-1) {
		l = str.find_first_of(delimiter,0);
		splitV.push_back(str.substr(0,l));
		str = str.substr(l+1,str.length());
	}
	return splitV;
}
