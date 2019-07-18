#include <iostream>     // std::cout
#include <algorithm>    // std::rotate
#include <vector>       // std::vector

using namespace std;
int main () {
  std::vector<int> myvector;
  std::vector<int> Ciaone;
  // set some values:
  for (int i=1; i<10; ++i) myvector.push_back(i);
  std::vector<int> myvector2 (myvector.size());// 1 2 3 4 5 6 7 8 9
  copy(myvector.begin(),myvector.end(),myvector2.begin());
  random_shuffle(myvector2.begin(),myvector2.end());
  
  
  for (int i=0; i<5; ++i) Ciaone.push_back(myvector[i]);
	
	int j=0;
	 while(Ciaone.size()<myvector.size()){
  			if(std::find(Ciaone.begin(),Ciaone.end(),myvector2[j]) == Ciaone.end())  Ciaone.push_back(myvector2[j]);
				j++;
	}

	myvector.front()=2;
	myvector[5]=3;
  copy(myvector.begin(),myvector.end(),myvector2.begin());
  // print out content:
  std::cout << "myvector contains:";
  for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it)
    std::cout << ' ' << *it;
  std::cout << '\n';
  
   for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it)
    std::cout << ' ' << *it;
  std::cout << '\n';

for(int i=0;i<myvector.size()-1;i++){
	if (std::find(myvector.begin()+i+1, myvector.end(), myvector[i]) != myvector.end())
			std::cout << "Element found"<<endl;
			break;

}



   /*for (std::vector<int>::iterator it=Ciaone.begin(); it!=Ciaone.end(); ++it)
    std::cout << ' ' << *it;
  std::cout << '\n';*/
  return 0;
}
