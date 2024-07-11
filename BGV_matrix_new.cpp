#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <helib/helib.h>
#include "matrix_new.h"
#include <vector>
#include <omp.h>

using namespace std;
using namespace NTL;

void matrix_square(int n)
{
  int m1 = 19, m2 = 17, m3 = 23;
  int p1 = m1 - 1, p2 = m2 - 1, p3 = m3 - 1;
  int m0 = m1 * m2 * m3;
  long p;
  int b_size = 16;
  int b = n / b_size;
  // if (m == 128 && l == 128 && n == 1) p = 3666356737;
  setp(p, m0);
  unsigned long r = 1;
  unsigned long bits = 80;
  unsigned long c = 2;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  long slot = context.getNSlots();
  std::cout << " m =  " << m1 << "*" << m2 << "*" << m3 << "=" << m0
            << "\n"
            << " n = "
            << slot
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  secret_key.GenSecKey();
  // context.printout();
  std::vector<long> lookup = gen_lookup(context, m0);
  std::vector<long> tr1 = Findgen(m1, m2, m3);
  std::vector<long> tr2 = Findgen(m2, m1, m3);
  std::vector<long> tr3 = Findgen(m3, m1, m2);
  long g1 = tr1[1], g2 = tr2[1], g3 = tr3[1];
  HELIB_NTIMER_START(KeyGen);
  GenTrace(secret_key, p2, tr2);
  HELIB_NTIMER_STOP(KeyGen);
  helib::printNamedTimer(std::cout, "KeyGen");
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;
  std::vector<long> power_new = gen_powernew(g1, g2, g3, m0, slot, p1, p2, p3);
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(n, n, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(n, n, 10, p, "M2.txt");
  std::cout << std::endl;
  cout << "Plaintext multiplication: M1*M2 mod " << p << " is equal to: " << endl;
  plain_text_multiplication(M1, M2, n, n, n, p);
  // Block matrices
  vector<vector<vector<vector<long>>>> b_M1(b, vector<vector<vector<long>>>(b, vector<vector<long>>(b_size, vector<long>(b_size, 0))));
  vector<vector<vector<vector<long>>>> b_M2(b, vector<vector<vector<long>>>(b, vector<vector<long>>(b_size, vector<long>(b_size, 0))));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int i0 = i / b_size;
      int j0 = j / b_size;
      b_M1[i0][j0][i % b_size][j % b_size] = M1[i][j];
      b_M2[i0][j0][i % b_size][j % b_size] = M2[i][j];
    }
  }
  HELIB_NTIMER_START(Encode);
  std::vector<std::vector<Ptxt_new>> Ptxt1 = std::vector<std::vector<Ptxt_new>>(b, std::vector<Ptxt_new>(b, Ptxt_new(public_key, lookup, power_new)));
  std::vector<std::vector<Ptxt_new>> Ptxt2 = std::vector<std::vector<Ptxt_new>>(b, std::vector<Ptxt_new>(b, Ptxt_new(public_key, lookup, power_new)));

#pragma omp parallel for
  for (int L = 0; L < b * b; L++)
  {
    int i = L / b;
    int j = L % b;
    encode(Ptxt1[i][j], b_M1[i][j], 0, b_size, b_size, p1, p2, p3);
    encode(Ptxt2[i][j], b_M2[i][j], 1, b_size, b_size, p1, p2, p3);
  }
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  HELIB_NTIMER_START(Encrypt);
  vector<vector<helib::Ctxt>> ctxt1(b, vector<helib::Ctxt>(b, helib::Ctxt(public_key)));
  vector<vector<helib::Ctxt>> ctxt2(b, vector<helib::Ctxt>(b, helib::Ctxt(public_key)));
#pragma omp parallel for
  for (int L = 0; L < b * b; L++)
  {
    int i = L / b;
    int j = L % b;
    public_key.Encrypt(ctxt1[i][j], Ptxt1[i][j].getptxt());
    public_key.Encrypt(ctxt2[i][j], Ptxt2[i][j].getptxt());
  }
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");
  std::cout << "Begin" << std::endl;
  vector<vector<helib::Ctxt>> ctxt_result(b, vector<helib::Ctxt>(b, helib::Ctxt(public_key)));
  HELIB_NTIMER_START(HomMatMult);
  for (int i = 0; i < b; i++)
  {
    for (int j = 0; j < b; j++)
    {
      for (int k = 0; k < b; k++)
      {
        helib::Ctxt temp = ctxt1[i][k];
        temp.multLowLvl(ctxt2[k][j]);
        ctxt_result[i][j].addCtxt(temp);
      }
      ctxt_result[i][j].reLinearize();
      ctxt_result[i][j] = FastTrace(ctxt_result[i][j], p2, tr2);
    }
  }
  HELIB_NTIMER_STOP(HomMatMult);
  helib::printNamedTimer(std::cout, "HomMatMult");
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result(b, vector<helib::Ptxt<helib::BGV>>(b, helib::Ptxt<helib::BGV>(public_key)));
  for (int i = 0; i < b; i++)
  {
    for (int j = 0; j < b; j++)
    {
      secret_key.Decrypt(plaintext_result[i][j], ctxt_result[i][j]);
    }
  }

  vector<vector<Ptxt_new>> Ptxt_result(b, vector<Ptxt_new>(b, Ptxt_new(plaintext_result[0][0], public_key, lookup, power_new)));
  for (int i = 0; i < b; i++)
  {
    for (int j = 0; j < b; j++)
    {
      Ptxt_result[i][j].ptxt_m = plaintext_result[i][j];
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      long a = long(Ptxt_result[i / b_size][j / b_size][p1 * (i % b_size) + (j % b_size)]);
      if (a > p / 2)
      {
        a = a - p;
      }
      std::cout << setw(log10(p)) << a << " ";
    }
    std::cout << std::endl;
  }
}
void matrix_nonsquare(int m, int l, int n)
{
  int m1, m2, m3;
  setm(m1, m2, m3, n, l, m);
  int p1 = m1 - 1, p2 = m2 - 1, p3 = m3 - 1;
  int m0 = m1 * m2 * m3;
  long p;
  // if (m == 128 && l == 128 && n == 1) p = 3666356737;
  setp(p, m0);
  unsigned long r = 1;
  unsigned long bits = 110;
  unsigned long c = 2;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  long slot = context.getNSlots();
  std::cout << " m =  " << m1 << "*" << m2 << "*" << m3 << "=" << m0
            << "\n"
            << " n = "
            << slot
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  secret_key.GenSecKey();
  // context.printout();
  std::vector<long> lookup = gen_lookup(context, m0);
  std::vector<long> tr1 = Findgen(m1, m2, m3);
  std::vector<long> tr2 = Findgen(m2, m1, m3);
  std::vector<long> tr3 = Findgen(m3, m1, m2);
  long g1 = tr1[1], g2 = tr2[1], g3 = tr3[1];
  // for (int i = 0;i<tr2.size();i++){
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  HELIB_NTIMER_START(KeyGen);
  // for (int i = 0;i<tr2.size();i++){
  //   //cout << tr2[i] << std::endl;
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  GenTrace(secret_key, p2, tr2);
  HELIB_NTIMER_STOP(KeyGen);
  helib::printNamedTimer(std::cout, "KeyGen");
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;
  std::vector<long> power_new = gen_powernew(g1, g2, g3, m0, slot, p1, p2, p3);
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(m, l, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(l, n, 10, p, "M2.txt");
  std::cout << std::endl;
  cout << "Plaintext multiplication: M1*M2 mod " << p << " is equal to: " << endl;
  plain_text_multiplication(M1, M2, m, l, n, p);
  Ptxt_new Ptxt1 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt2 = Ptxt_new(public_key, lookup, power_new);
  encode(Ptxt1, M1, 0, m, l, p1, p2, p3);
  encode(Ptxt2, M2, 1, l, n, p1, p2, p3);
  helib::Ctxt ctxt1(public_key), ctxt2(public_key);
  public_key.Encrypt(ctxt1, Ptxt1.getptxt());
  public_key.Encrypt(ctxt2, Ptxt2.getptxt());
  std::cout << "Begin" << std::endl;
  helib::Ctxt ctxt_result(public_key);
  HELIB_NTIMER_START(HomMatMult);
  ctxt2.multiplyBy(ctxt1);
  ctxt_result = FastTrace(ctxt2, p2, tr2);
  HELIB_NTIMER_STOP(HomMatMult);
  helib::printNamedTimer(std::cout, "HomMatMult");
  helib::Ptxt<helib::BGV> plaintext_result = helib::Ptxt<helib::BGV>(context);
  //std::cout << ctxt_result.capacity() << std::endl;
  secret_key.Decrypt(plaintext_result, ctxt_result);
  Ptxt_new Ptxt_result = Ptxt_new(plaintext_result, public_key, lookup, power_new);
  for (int k = 0; k < 1; k++)
  {
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        long long a = long(Ptxt_result[p1 * p3 * k + p1 * i + j]);
        if (a>p/2) {a = a-p;}
        std::cout << setw(log10(p)) << a << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void example1()
{
  int m = 4;
  int l = 100;
  int n = 50;
  int m1, m2, m3;
  setm(m1, m2, m3, n, l, m);
  int p1 = m1 - 1, p2 = m2 - 1, p3 = m3 - 1;
  int m0 = m1 * m2 * m3;
  long p;
  // if (m == 128 && l == 128 && n == 1) p = 3666356737;
  setp(p, m0);
  unsigned long r = 1;
  unsigned long bits = 300;
  unsigned long c = 2;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  long slot = context.getNSlots();
  std::cout << " m =  " << m1 << "*" << m2 << "*" << m3 << "=" << m0
            << "\n"
            << " n = "
            << slot
            << "\n"
            << "bit of plaintext"
            << log2(p)
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  secret_key.GenSecKey();
  // context.printout();
  std::vector<long> lookup = gen_lookup(context, m0);
  std::vector<long> tr1 = Findgen(m1, m2, m3);
  std::vector<long> tr2 = Findgen(m2, m1, m3);
  std::vector<long> tr3 = Findgen(m3, m1, m2);
  long g1 = tr1[1], g2 = tr2[1], g3 = tr3[1];
  // for (int i = 0;i<tr2.size();i++){
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  HELIB_NTIMER_START(KeyGen);
  // for (int i = 0;i<tr2.size();i++){
  //   //cout << tr2[i] << std::endl;
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  GenTrace(secret_key, p2, tr2);
  GenTrace(secret_key, p1, tr1);
  HELIB_NTIMER_STOP(KeyGen);
  helib::printNamedTimer(std::cout, "KeyGen");
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;
  std::vector<long> power_new = gen_powernew(g1, g2, g3, m0, slot, p1, p2, p3);
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(4, 100, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(100, 50, 10, p, "M2.txt");
  std::cout << std::endl;
  std::cout << "The matrix M3 is equal to: " << std::endl;
  std::vector<std::vector<long>> M3 = randmat(50, 25, 10, p, "M3.txt");
  std::cout << std::endl;
  std::cout << "The matrix M4 is equal to: " << std::endl;
  std::vector<std::vector<long>> M4 = randmat(25, 5, 10, p, "M4.txt");
  std::cout << std::endl;
  std::cout << "The matrix M5 is equal to: " << std::endl;
  std::vector<std::vector<long>> M5 = randmat(5, 1, 10, p, "M5.txt");
  std::cout << std::endl;

  std::vector<std::vector<long>> M;
  M = plain_text_multiplication_output(M1, M2, 4, 100, 50, p);
  M = plain_text_multiplication_output(M, M3, 4, 50, 25, p);
  M = plain_text_multiplication_output(M, M4, 4, 25, 5, p);
  M = plain_text_multiplication_output(M, M5, 4, 5, 1, p);
  cout << "Plaintext multiplication: M1*M2*M3*M4*M5 mod " << p << " is equal to: " << endl;
  printf_matrix(M, 4, 1, log10(p) + 2, p);
  Ptxt_new Ptxt1 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt2 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt3 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt4 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt5 = Ptxt_new(public_key, lookup, power_new);
  encode(Ptxt1, M1, 0, 4, 100, p1, p2, p3);
  encode(Ptxt2, M2, 1, 100, 50, p1, p2, p3);
  encode(Ptxt3, M3, 2, 50, 25, p1, p2, p3);
  encode(Ptxt4, M4, 1, 25, 5, p1, p2, p3);
  encode(Ptxt5, M5, 2, 5, 1, p1, p2, p3);

  helib::Ctxt ctxt1(public_key), ctxt2(public_key), ctxt3(public_key), ctxt4(public_key), ctxt5(public_key);
  public_key.Encrypt(ctxt1, Ptxt1.getptxt());
  public_key.Encrypt(ctxt2, Ptxt2.getptxt());
  public_key.Encrypt(ctxt3, Ptxt3.getptxt());
  public_key.Encrypt(ctxt4, Ptxt4.getptxt());
  public_key.Encrypt(ctxt5, Ptxt5.getptxt());
  std::cout << "Begin" << std::endl;
  helib::Ctxt ctxt_result(public_key);
  HELIB_NTIMER_START(HomMatMult);
  ctxt2.multiplyBy(ctxt1);
  ctxt_result = FastTrace(ctxt2, p2, tr2);
  ctxt_result.multiplyBy(ctxt3);
  ctxt_result = FastTrace(ctxt_result, p1, tr1);
  ctxt_result.multiplyBy(ctxt4);
  ctxt_result = FastTrace(ctxt_result, p2, tr2);
  ctxt_result.multiplyBy(ctxt5);
  ctxt_result = FastTrace(ctxt_result, p1, tr1);
  HELIB_NTIMER_STOP(HomMatMult);
  helib::printNamedTimer(std::cout, "HomMatMult");
  helib::Ptxt<helib::BGV> plaintext_result = helib::Ptxt<helib::BGV>(context);
  secret_key.Decrypt(plaintext_result, ctxt_result);
  Ptxt_new Ptxt_result = Ptxt_new(plaintext_result, public_key, lookup, power_new);
   cout << "Cyphertext multiplication: M1*M2*M3*M4*M5 is equal to: " << endl;
  for (int k = 0; k < 1; k++)
  {
    for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 1; j++)
      {
        int a = long(Ptxt_result[p1 * p3 * k + p1 * i + j]) > p / 2 ? long(Ptxt_result[p1 * p3 * k + p1 * i + j]) - p : long(Ptxt_result[p1 * p3 * k + p1 * i + j]);
        std::cout << setw(log10(p)) << a << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void example2()
{
  int m = 5;
  int l = 50;
  int n = 50;
  int m1, m2, m3;
  setm(m1, m2, m3, n, l, m);
  int p1 = m1 - 1, p2 = m2 - 1, p3 = m3 - 1;
  int m0 = m1 * m2 * m3;
  long p = 21889001;
  // if (m == 128 && l == 128 && n == 1) p = 3666356737;
  unsigned long r = 1;
  unsigned long bits = 250;
  unsigned long c = 2;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  long slot = context.getNSlots();
  std::cout << " m =  " << m1 << "*" << m2 << "*" << m3 << "=" << m0
            << "\n"
            << " n = "
            << slot
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  secret_key.GenSecKey();
  // context.printout();
  std::vector<long> lookup = gen_lookup(context, m0);
  std::vector<long> tr1 = Findgen(m1, m2, m3);
  std::vector<long> tr2 = Findgen(m2, m1, m3);
  std::vector<long> tr3 = Findgen(m3, m1, m2);
  long g1 = tr1[1], g2 = tr2[1], g3 = tr3[1];
  // for (int i = 0;i<tr2.size();i++){
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  HELIB_NTIMER_START(KeyGen);
  // for (int i = 0;i<tr2.size();i++){
  //   //cout << tr2[i] << std::endl;
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  GenTrace(secret_key, p2, tr2);
  GenTrace(secret_key, p1, tr1);
  HELIB_NTIMER_STOP(KeyGen);
  helib::printNamedTimer(std::cout, "KeyGen");
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;
  std::vector<long> power_new = gen_powernew(g1, g2, g3, m0, slot, p1, p2, p3);
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(5, 50, 50, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(50, 50, 10, p, "M2.txt");
  std::cout << std::endl;
  std::cout << "The matrix M3 is equal to: " << std::endl;
  std::vector<std::vector<long>> M3 = randmat(50, 50, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M4 is equal to: " << std::endl;
  std::vector<std::vector<long>> M4 = randmat(50, 50, 10, p, "M2.txt");
  std::cout << std::endl;
  std::cout << "The matrix M5 is equal to: " << std::endl;
  std::vector<std::vector<long>> M5 = randmat(50, 50, 10, p, "M2.txt");
  std::cout << std::endl;

  std::vector<std::vector<long>> M;
  M = plain_text_multiplication_output(M1, M2, 5, 50, 50, p);
  M = plain_text_multiplication_output(M, M3, 5, 50, 50, p);
  M = plain_text_multiplication_output(M, M4, 5, 50, 50, p);
  M = plain_text_multiplication_output(M, M5, 5, 50, 50, p);
  cout << "Plaintext multiplication: M1*M2*M3*M4*M5 is equal to: " << endl;
  printf_matrix(M, 5, 50, log10(p) + 2, p);
  Ptxt_new Ptxt1 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt2 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt3 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt4 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt5 = Ptxt_new(public_key, lookup, power_new);
  encode(Ptxt1, M1, 0, 5, 50, p1, p2, p3);
  encode(Ptxt2, M2, 1, 50, 50, p1, p2, p3);
  encode(Ptxt3, M3, 2, 50, 50, p1, p2, p3);
  encode(Ptxt4, M4, 1, 50, 50, p1, p2, p3);
  encode(Ptxt5, M5, 2, 50, 50, p1, p2, p3);

  helib::Ctxt ctxt1(public_key), ctxt2(public_key), ctxt3(public_key), ctxt4(public_key), ctxt5(public_key);
  public_key.Encrypt(ctxt1, Ptxt1.getptxt());
  public_key.Encrypt(ctxt2, Ptxt2.getptxt());
  public_key.Encrypt(ctxt3, Ptxt3.getptxt());
  public_key.Encrypt(ctxt4, Ptxt4.getptxt());
  public_key.Encrypt(ctxt5, Ptxt5.getptxt());
  std::cout << "Begin" << std::endl;
  helib::Ctxt ctxt_result(public_key);
  HELIB_NTIMER_START(HomMatMult);
  ctxt2.multiplyBy(ctxt1);
  ctxt_result = FastTrace(ctxt2, p2, tr2);
  ctxt_result.multiplyBy(ctxt3);
  ctxt_result = FastTrace(ctxt_result, p1, tr1);
  ctxt_result.multiplyBy(ctxt4);
  ctxt_result = FastTrace(ctxt_result, p2, tr2);
  ctxt_result.multiplyBy(ctxt5);
  ctxt_result = FastTrace(ctxt_result, p1, tr1);
  HELIB_NTIMER_STOP(HomMatMult);
  helib::printNamedTimer(std::cout, "HomMatMult");
  helib::Ptxt<helib::BGV> plaintext_result = helib::Ptxt<helib::BGV>(context);
  secret_key.Decrypt(plaintext_result, ctxt_result);
  Ptxt_new Ptxt_result = Ptxt_new(plaintext_result, public_key, lookup, power_new);
  cout << "Cyphertext multiplication: M1*M2*M3*M4*M5 is equal to: " << endl;
  for (int i = 0; i < 5; i++)
  {
    for (int k = 0; k < 50; k++)
    {
      for (int j = 0; j < 1; j++)
      {
        int a = long(Ptxt_result[p1 * p3 * k + p1 * i + j]) > p / 2 ? long(Ptxt_result[p1 * p3 * k + p1 * i + j]) - p : long(Ptxt_result[p1 * p3 * k + p1 * i + j]);
        std::cout << setw(log10(p)) << a << " ";
      }
    }
    std::cout << std::endl;
  }
}

void example3()
{
  int m = 4;
  int l = 50;
  int n = 50;
  int m1, m2, m3;
  setm(m1, m2, m3, n, l, m);
  int p1 = m1 - 1, p2 = m2 - 1, p3 = m3 - 1;
  int m0 = m1 * m2 * m3;
  long p;
  setp(p,m0);
  // if (m == 128 && l == 128 && n == 1) p = 3666356737;
  unsigned long r = 1;
  unsigned long bits = 240;
  unsigned long c = 2;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  long slot = context.getNSlots();
  std::cout << " m =  " << m1 << "*" << m2 << "*" << m3 << "=" << m0
            << "\n"
            << " n = "
            << slot
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  secret_key.GenSecKey();
  // context.printout();
  std::vector<long> lookup = gen_lookup(context, m0);
  std::vector<long> tr1 = Findgen(m1, m2, m3);
  std::vector<long> tr2 = Findgen(m2, m1, m3);
  std::vector<long> tr3 = Findgen(m3, m1, m2);
  long g1 = tr1[1], g2 = tr2[1], g3 = tr3[1];
  // for (int i = 0;i<tr2.size();i++){
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  HELIB_NTIMER_START(KeyGen);
  // for (int i = 0;i<tr2.size();i++){
  //   //cout << tr2[i] << std::endl;
  //   secret_key.GenKeySWmatrix(1, tr2[i], 0, 0);
  // }
  GenTrace(secret_key, p2, tr2);
  GenTrace(secret_key, p1, tr1);
  HELIB_NTIMER_STOP(KeyGen);
  helib::printNamedTimer(std::cout, "KeyGen");
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;
  std::vector<long> power_new = gen_powernew(g1, g2, g3, m0, slot, p1, p2, p3);
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(4, 50, 10, p, "M1.txt");
  std::vector<std::vector<long>> M2 = randmat_sym(50, 10,p, "M2.txt");
  std::cout << std::endl;

  std::vector<std::vector<long>> M;
  M = plain_text_multiplication_output(M1, M2, 4, 50, 50, p);
  for (int i = 0;i<3;i++){
     M = plain_text_multiplication_output(M, M2, 4, 50, 50, p); 
  }
  cout << "Plaintext multiplication: M1*M2*M2*M2*M2 is equal to: " << endl;
  printf_matrix(M, 4, 50, log10(p) + 2, p);
  Ptxt_new Ptxt1 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt2 = Ptxt_new(public_key, lookup, power_new);
  Ptxt_new Ptxt3 = Ptxt_new(public_key, lookup, power_new);
  encode(Ptxt1, M1, 0, 4, 50, p1, p2, p3);
  encode(Ptxt2, M2, 1, 50, 50, p1, p2, p3);
  encode(Ptxt3, M2, 2, 50, 50, p1, p2, p3);

  helib::Ctxt ctxt1(public_key), ctxt2(public_key), ctxt3(public_key), ctxt4(public_key), ctxt5(public_key);
  public_key.Encrypt(ctxt1, Ptxt1.getptxt());
  public_key.Encrypt(ctxt2, Ptxt2.getptxt());
  public_key.Encrypt(ctxt3, Ptxt3.getptxt());
  std::cout << "Begin" << std::endl;
  helib::Ctxt ctxt_result(public_key);
  HELIB_NTIMER_START(HomMatMult);
  ctxt1.multiplyBy(ctxt2);
  ctxt_result = FastTrace(ctxt1, p2, tr2);
  ctxt_result.multiplyBy(ctxt3);
  ctxt_result = FastTrace(ctxt_result, p1, tr1);
  ctxt_result.multiplyBy(ctxt2);
  ctxt_result = FastTrace(ctxt_result, p2, tr2);
  ctxt_result.multiplyBy(ctxt3);
  ctxt_result = FastTrace(ctxt_result, p1, tr1);
  HELIB_NTIMER_STOP(HomMatMult);
  helib::printNamedTimer(std::cout, "HomMatMult");
  helib::Ptxt<helib::BGV> plaintext_result = helib::Ptxt<helib::BGV>(context);
  secret_key.Decrypt(plaintext_result, ctxt_result);
  Ptxt_new Ptxt_result = Ptxt_new(plaintext_result, public_key, lookup, power_new);
  cout << "Cyphertext multiplication: M1*M2*M2*M2*M2 mod " << p << " is equal to: " << endl;
  for (int i = 0; i < 4; i++)
  {
    for (int k = 0; k < 50; k++)
    {
      for (int j = 0; j < 1; j++)
      {
        int a = long(Ptxt_result[p1 * p3 * k + p1 * i + j]) > p / 2 ? long(Ptxt_result[p1 * p3 * k + p1 * i + j]) - p : long(Ptxt_result[p1 * p3 * k + p1 * i + j]);
        std::cout << setw(log10(p)) << a << " ";
      }
    }
    std::cout << std::endl;
  }
}

int main()
{
  int m, l, n, thread, t;
  cout << "            ------------------------------- " << endl;
  cout << "             Momomorpic Matrix Operations" << endl;
  cout << "            ------------------------------- " << endl;

  cout << "Please enter 0 or 1 (0 for square, 1 for nonsquare)" << endl;
  cin >> t;
  if (t == 1)
  {
    cout << "Please input the size of matrix (m*l*n):" << endl;
    cin >> m;
    cin >> l;
    cin >> n;
    matrix_nonsquare(m, l, n);
  }
  else
  {
    cout << "Please input the size of matrix (n >=16):" << endl;
    cin >> n;
    matrix_square(n);
  }
  // example3();
  // example1();
  // eigvector(2);
  return 0;
}
