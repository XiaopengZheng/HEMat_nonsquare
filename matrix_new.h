#include <helib/helib.h>
#include <vector>
#include <fstream>

class Ptxt_new
{
private:
   std::vector<long> lookup_m;
   std::vector<long> power_new_m;
public:
   helib::Ptxt<helib::BGV> ptxt_m;
   helib::Ptxt<helib::BGV> getptxt(){
      return ptxt_m;
   }
   Ptxt_new(const helib::PubKey &public_key, std::vector<long> lookup,std::vector<long> power_new){
      helib::Ptxt<helib::BGV> ptxt(public_key);
      ptxt_m = ptxt;
      power_new_m = power_new;
      lookup_m = lookup;
   }
   Ptxt_new(helib::Ptxt<helib::BGV> ptxt, const helib::PubKey &public_key, std::vector<long> lookup,std::vector<long> power_new){
      ptxt_m = ptxt;
      power_new_m = power_new;
      lookup_m = lookup;
   }
   helib::BGV::SlotType& operator[](long i){
      return ptxt_m[lookup_m[power_new_m[i]]];
   }
};

void setp(long &p, long m0);

void setm(int &m1, int &m2, int &m3, int m, int l, int n);

std::vector<long> gen_lookup(helib::Context& context, long m0);

std::vector<long> gen_powernew(long g1, long g2, long g3, long m0, long slot, int p1, int p2, int p3);

std::vector<long> Findgen(int m1, int m2,int m3);

bool isPrime(long num);

int Modpower(int a, int b, int m);

int FindPrime(long m);

int Modpower(int a,int b,int m);

void plain_text_multiplication(std::vector<std::vector<long>> M1, std::vector<std::vector<long>> M2, int m, int l, int n, long p);

std::vector<std::vector<long>> plain_text_multiplication_output(std::vector<std::vector<long>> M1, std::vector<std::vector<long>> M2, int m, int l, int n, long p);

void plain_text_multiplication(std::vector<std::vector<long>> M1, std::vector<std::vector<long>> M2, int m, int l, int n, long p, int s);

std::vector<std::vector<long>> randmat(int m, int n, int q, long p, std::string s);

std::vector<std::vector<long>> randmat_sym(int n, int q, long p, std::string s);

void printf_matrix(const std::vector<std::vector<long>> &M, int m, int n, int w, long &p);

void printf_matrix(const std::vector<std::vector<long>> &M, int m, int n, int w, long &p, std::ofstream &fout);

void encode(Ptxt_new& Ptxt, std::vector<std::vector<long>> M, int t, int m, int l, int p1, int p2, int p3);

helib::Ctxt FastTrace(helib::Ctxt ctxt, long r, std::vector<long> tr);

void GenTrace(helib::SecKey& secret_key, long r, std::vector<long> tr);