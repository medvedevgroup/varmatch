#ifndef UTILITIES_H
#define UTILITIES_H

#ifdef DEBUG
#define dout cout
#else
#define dout 0 && cout
#endif

#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>

using namespace std;

/*split function*/
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (! item.empty()){
			elems.push_back(item);
		}
	}
	return elems;
}

/*This split function only support char as delim, string as delim please boost split function*/
std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

inline bool FileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

template <typename T>
inline
string ToString(const T & number){
	string String = static_cast<ostringstream*>( &(ostringstream() << number) )->str();
	return String;
}

void show_time()
{
 time_t nowtime;
 nowtime = time(NULL); //get int time number
 struct tm * ptm=localtime(&nowtime);  //get system time
 cout << ptm->tm_mon+1 << "/" << ptm->tm_mday << "/"<< ptm->tm_year+1900 << "," ;
 cout << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec <<" " << endl;
}

#endif