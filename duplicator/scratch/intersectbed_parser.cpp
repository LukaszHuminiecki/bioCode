#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <set>

using namespace std;


struct duplicationEvent {
int node;
string refseq1;
string refseq2;
};

enum taxon {HUMAN, Homo_Pan_Gorilla, Catarrhini, Eutheria , Theria, Amniota, Tetrapoda, Euteleostomi, Chordata, Deuterostomi, Bilateria, Eumetazoa, Metazoa, Fungi_Metazoa, Eukaryota} orderedTaxa;
enum color {red, orange, yellow, white, green, blue, brown, grey} orderedTaxaColors;
enum shape {circle, box} orderedTaxaShape;
taxon operator++(taxon& t) {t= (t==Eukaryota) ? HUMAN : taxon(t+1);return t;}
color operator++(color& c) {c= (c==grey) ? red : color(c+1);return c;}
shape operator++(shape& s) {s= (s==box) ? circle : shape(s+1);return s;}


struct duplicationWave {
vector<string> refseqs; // all refseqs mapping to given taxon
string t; // change to enum taxon
};


vector<duplicationEvent> duplicationWave2;


int main() {
  ifstream inputFile("intersectbed_refseq_Tfbs");
  string tempLine;
void printMap(map <string, string>);
void printMap(multimap <string, string>);
void printMap(map <string, set<string> >);
multimap <string, string> m;
map <string, set<string> > m2;
set<string> s;
  while ( getline(inputFile, tempLine, '\n') ) {
    vector <string> tempstr;
    stringstream ss(tempLine);
    string temp;
    while (getline(ss, temp, '\t')) {  
      tempstr.push_back(temp);
    }
m.insert(std::make_pair(tempstr[3], tempstr[15]));
m2.find(tempstr[3])->second.insert(tempstr[15]);

}  
//printMap(m);
double jaccard(string, string, multimap <string, string>);
cout << jaccard("NM_000815", "NM_000983", m) << endl;
double ji = 0;  // calculate jaccard correctly
return ji; 
}


double jaccard(string rf1, string rf2, multimap <string, string> m) {
double j;

//cout << m[rf1] << endl;
return j;
}


int jaccard(map <string, string> m, string refseq1, string refseq2) {
int jaccard;
// use std::set_intersection STL algorithm
// http://www.cplusplus.com/reference/algorithm/set_intersection/
// or do you need to iterate over both of them?
return jaccard;
}


void printMap(map <string, set<string> > m)
{
//for (map <string, set<string> >::const_iterator i = m.begin(); i != m.end(); i++)
//	{
//
//	for (set<string>::const_iterator ii = i->second.begin(); ii != i->second.end(); ii++) {
//		cout << "key:\t" << i->first << "\tvalue:\t" << *ii << endl;
//		}
//	}
}


void printMap(map <string, string> m)
{
for (map<string,string>::const_iterator i = m.begin(); i != m.end(); i++)
	{
	cout << "key:\t" << i->first << "\tvalue:\t" << i->second << endl;
	}
}


void printMap(multimap <string, string> m)
{
	for (multimap<string,string>::const_iterator i = m.begin(); i != m.end(); i++)
	{
	cout << "key:\t" << i->first << "\tvalue:\t" << i->second << endl;
	}
}
