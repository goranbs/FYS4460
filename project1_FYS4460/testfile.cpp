#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(){

  vector <int> particles;
  int N = 10;

  for(int i=0;i<N;++i){
    particles.push_back(i);
  }
  if (find(particles.begin(),particles.end(),3)==particles.end()){
    cout << particles.end() << " " << particles.begin() << endl;
  }
  return 0;
}
