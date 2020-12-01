#include <map>
#include <iostream>
#include <chrono>
#include <malloc.h>

int main() {
   auto t1 = std::chrono::high_resolution_clock::now();
   int count = 1000000;
   int i;
   std::map<int,int> m;
   /* int* m_array = (int*) malloc(count * sizeof(int)); */
   /* int* orig = m_array; */

   for (i=0;i<count;i++) {
      int a = i;
      int b = 1000;
      m.insert(std::pair<int,int>(a,b));
      /* *m_array = i; */
      /* m_array++; */
   }

   std::map<int,int>::iterator it;
   int sum=0;

   for (it=m.begin(); it!=m.end(); it++) {
      sum += it->second;
      /* sum += m_array[500]; */
   }

   /* free(orig); */

   std::cout << "Sum: " << sum << std::endl;

   auto t2 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
   std::cout << "run time in milliseconds: " << fp_ms.count() << std::endl;
   return 0;
}
