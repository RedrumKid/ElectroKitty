#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void printV(vector<vector<int>>& v) {
  	for (auto i: v) {
        for (auto j: i) {
            cout << j << " ";
        }
        cout << endl;
    }
  	cout << endl;
}

int main() {
  
  	vector<vector<vector<int>>> v2(2, vector<vector<int>>(3, vector<int>(3, 1)));

    for(int i = 0; i<v2.size(); i++){
          	printV(v2[i]);
    }
    return 0;
}