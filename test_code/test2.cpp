#include <map>
#include <iostream>
#include <chrono>
#include <fstream>

int main() {
   std::ofstream* outfile;
   outfile = new std::ofstream();
   outfile->open("test.dat", std::ios::out | std::ios::binary | std::ios::trunc);

   int a = 5;
   double b = 6.0;
   double c = 7.0;
   int d = 8;

   outfile->write(reinterpret_cast<const char *>(&a), sizeof(int));
   outfile->write(reinterpret_cast<const char *>(&b), sizeof(double));
   outfile->write(reinterpret_cast<const char *>(&c), sizeof(double));
   outfile->write(reinterpret_cast<const char *>(&d), sizeof(int));
   outfile->close();

   //reset all variable values
   a = 0;
   b = 0;
   c = 0;
   d = 0;

   int* iptr;
   double* dptr;
   unsigned char arr[100];

   std::ifstream *infile;
   infile = new std::ifstream();
   infile->open("test.dat", std::ios::in | std::ios::binary);

   infile->read(reinterpret_cast<char *>(arr), sizeof(int)*2 + sizeof(double) * 2);  
   iptr = (int*) arr;

   a = *iptr;
   iptr++;

   dptr = (double*) iptr;
   b = *dptr;
   dptr++;

   c = *dptr;
   dptr++;

   iptr = (int*) dptr;
   d = *iptr;

   std::cout << a << std::endl;
   std::cout << b << std::endl;
   std::cout << c << std::endl;
   std::cout << d << std::endl;

   infile->close();

   return 0;
}
