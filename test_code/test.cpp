#include <map>
#include <iostream>
#include <chrono>

int main() {
   auto t1 = std::chrono::high_resolution_clock::now();
   int count = 1000000;
   int i;
   std::map<int,int> m;

   for (i=0;i<count;i++) {
      int a = i;
      int b = 1000;
      m.insert(std::pair<int,int>(a,b));
   }

   std::map<int,int>::iterator it;
   int sum=0;

   for (it=m.begin(); it!=m.end(); it++) {
      sum += it->second;
   }

   std::cout << "Sum: " << sum << std::endl;

   auto t2 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
   std::cout << "run time in milliseconds: " << fp_ms.count() << std::endl;
   return 0;
}
