/*
*StringHelper class contains several methods to deal with string
*Max Planck Institute for Developmental Biology and MPI for Biological Cybernetics
*@author: Dominik Grimm
*@year: 2010
*/
#ifndef StringHelper_CLASS
#define StringHelper_CLASS

#include <string>
#include <vector>
#include <sstream>
#include "CEasyGWAS/types.h"

using namespace std;

class StringHelper {
	public:
		StringHelper() {};
		~StringHelper() {};

		static vector<string> split(string ,const string&);
		/*
		*@Desciption: string to type <float>, <double>, <int>, <bool>, <long> ...
		*@param: variable from typ T to store the transformed string
		*        the string to transform
		*@return: true if everything was ok, otherwise false
		*/
		template <typename T>
		static bool string_to(T& t, const string& s) {
			istringstream iss(s);
 			return !(iss >> t).fail();
		}
		
		/*
		*@Desciption: string to type <float>, <double>, <int>, <bool>, <long> ...
		*@param:  the string to transform
		*@return: the transformed type
		*/
		template <typename T>
		static T string_to( const string& s) {
			istringstream iss(s);
			T t;
			iss >> t;
 			return t;
		}

		/*
		*@Desciption: type <float>, <double>, <int>, <bool>, <long> to string ...
		*@param: variable from typ T to transform 
		*        the transformed string
		*@return: true if everything was ok, otherwise false
		*/
		template <typename T>
		static bool to_string(string& s, const T& t) {
			stringstream iss;
			bool flag =(iss << t);
			s = iss.str();
 			return flag;
		}
	
		/*
		*@Desciption: type <float>, <double>, <int>, <bool>, <long> to string ...
		*@param: variable from typ T to transform
		*@return: string
		*/
		template <typename T>
		static string to_string(const T& t) {
			stringstream iss;
			iss << t;
			return iss.str();
		}

		static string trim(const string& str,
				   const string& whitespace = " \t\n") {
			const size_t strBegin = str.find_first_not_of(whitespace);
		        if (strBegin == string::npos)
				return ""; // no content

			const size_t strEnd = str.find_last_not_of(whitespace);
			const size_t strRange = strEnd - strBegin + 1;
			return str.substr(strBegin, strRange);
		}

		/*
		*@Desciption: trim string
		*@param: string to trim
		*/
		static string& reduce(string& s) {
			size_t pos;
			while((pos = s.find(' ')) != string::npos) {
				s.erase(pos,1);
			}
			return s;	
		}
		
		/*
		*@Desciption: trim string with a certain pattern
		*@param: string to trim, trim pattern
		*/
		static string& reduce(string& s, const string& pattern) {
			size_t pos;
			while((pos = s.find(pattern)) != string::npos) {
				s.erase(pos,1);
			}
			return s;	
		}
		
		/*
		*@Desciption: count number of characters c in s
		*@param: number of characters
		*/
		static uint get_num_characters(string& s, const string& c) {
			uint num = 0;
			for(uint i=0; i<s.length(); i++) {
				if(s[i] == c[0]) num++;
			}
			return num;	
		}
	
		/*
		*@Desciption: string to upper case 
		*@param: string s
		*/
		static string& to_upper(string& s) {
			for(unsigned int i=0; i<s.length();i++) {
				s[i] = std::toupper(s[i]);
			}
			return s;
		}
		
		/*
		*@Desciption: string to lower case 
		*@param: string s
		*/
		static string& to_lower(string& s) {
			for(unsigned int i=0; i<s.length();i++) {
				s[i] = std::tolower(s[i]);
			}
			return s;
		}
};
#endif //StringHelper_CLASS
