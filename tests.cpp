#include <vector>
#include <iostream>

using namespace std;

int main(){

    vector<vector<vector<double>>> a(3);

    for(int i = 0; i < 3; i++){
        a[i][i].push_back(double(i));
    }

    for(int i = 0; i < 3; i++){
        cout<<a[i][i][0]<<endl;
    }

    return 0;
}