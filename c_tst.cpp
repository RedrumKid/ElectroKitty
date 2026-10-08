#include <vector>
#include <iostream>

using namespace std;

void update(vector<int>& pa){
    pa.push_back(1);
}

int main(){

    vector<int> vec;

    update(vec);

    for(int i = 0; i < int(vec.size()); i++){
        cout<<vec[i];
    }

    return 0;
}