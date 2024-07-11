
#include "matrix_new.h"
#include <helib/helib.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <helib/helib.h>
#include <vector>
#include <fstream>
#include <omp.h>

using namespace std;
using namespace NTL;

std::vector<long> gen_powernew(long g1, long g2, long g3, long m0, long slot, int p1, int p2, int p3)
{
  std::vector<long> power_new;
  for (int L = 0; L < slot; L++)
  {
    int i = L / (p1 * p3);
    int j = L % (p1 * p3) / p1;
    int k = L % p1;
    long a = Modpower(g2, i, m0);
    long b = Modpower(g3, j, m0);
    long c = Modpower(g1, k, m0);
    long pro = ((a * b * c) % m0);
    if (pro < 0)
    {
      pro = pro + m0;
    }
    power_new.push_back(pro);
  }
  return power_new;
}

std::vector<long> Findgen(int m1, int m2, int m3)
{
  std::vector<long> tr;
  long m = m1 * m2 * m3;
  for (int i = 1; i <= m1; i++)
  {
    if ((i * m2 * m3 + 1) % m1 != 0)
    {
      tr.push_back((i * m2 * m3 + 1) % m);
    }
  }
  for (int i = 0; i < tr.size(); i++)
  {
    std::vector<long> tr0;
    tr0.push_back(1);
    int j;
    for (j = 1; j < m1 - 1; j++)
    {
      // std::cout << Modpower(tr[i],j,m) <<  " ";
      int a = Modpower(tr[i], j, m);
      tr0.push_back(Modpower(tr[i], j, m));
      if (a == 1 || (Modpower(tr[i], j, m) == Modpower(tr[i], j - 1, m)))
      {
        break;
      }
    }
    // std::cout << std::endl;
    if (j == m1 - 1)
    {
      // std::cout << tr[i] << " ";
      return tr0;
      break;
    }
  }
  return tr;
}

bool isPrime(long num)
{
  if (num == 1 || num == 4)
    return 0;
  if (num == 2 || num == 3)
    return 1;
  if (num % 6 != 1 && num % 6 != 5)
    return 0;
  int tmp = sqrt(num);
  for (int i = 5; i <= tmp; i += 6)
  {
    if (num % i == 0 || num % (i + 2) == 0)
      return 0;
  }
  return 1;
}

void setp(long &p, long m0)
{
  long k = 65537/m0;
  while (!isPrime(m0 * k + 1) || log2(m0 * k + 1) < 16)
  {
    k++;
  }
  p = m0 * k + 1;
  std::cout << "number of bit of plaintext: " << log2(p) << std::endl;
}

void setm(int &m1, int &m2, int &m3, int m, int l, int n)
{
  m1 = FindPrime(m);
  m2 = FindPrime(l);
  while (m2 == m1)
  {
    m2 = FindPrime(m2);
  }
  m3 = FindPrime(n);
  while (m3 == m2 || m3 == m1)
  {
    m3 = FindPrime(m3);
  }
}

std::vector<long> gen_lookup(helib::Context &context, long m0)
{
  int num = context.getEA().getPAlgebra().numOfGens();
  long slot = context.getNSlots();
  std::vector<long> gen, ord;
  for (int i = 0; i < num; i++)
  {
    gen.push_back(context.getEA().getPAlgebra().genToPow(i, 1));
    ord.push_back(context.getZMStar().OrderOf(i));
  }
  std::vector<long> power_old = std::vector<long>(slot, 1);
  int P = 1;
  for (int M = 0; M < num; M++)
  {
    P = P * ord[M];
    for (int L = 0; L < slot; L++)
    {
      power_old[L] = (power_old[L] * Modpower(gen[M], L / (slot / P), m0)) % m0;
    }
  }
  std::vector<long> lookup = std::vector<long>(m0 + 1, 0);
  for (int L = 0; L < slot; L++)
  {
    lookup[power_old[L]] = L;
  }
  return lookup;
}

int Modpower(int a, int b, int m)
{
  if (b == 0)
  {
    return 1;
  }
  if (b == 1)
  {
    return a;
  }
  else
  {
    if (b != 1)
    {
      return (a * Modpower(a, b - 1, m)) % m;
    }
    else
    {
      int t = Modpower(a, b - 1, m);
      return (t * t) % m;
    }
  }
}

int FindPrime(long m)
{
  int k = m + 1;
  while (!isPrime(k))
  {
    k++;
  }
  return (k);
}

vector<vector<long>> plain_text_multiplication_output(vector<vector<long>> M1, vector<vector<long>> M2, int m, int l, int n, long p)
{
  vector<vector<long>> M(m, vector<long>(n));
  // long p = INT32_MAX;
  std::ofstream fout;
  fout.open("M1*M2_plaintext.txt", ios::out);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      M[i][j] = 0;
      for (int k = 0; k < l; k++)
      {
        M[i][j] = M[i][j] + M1[i][k] * M2[k][j];
      }
      if (M[i][j]>p/2){
        M[i][j] = M[i][j] - p;
      }
    }
  }
  return M;
}

void plain_text_multiplication(vector<vector<long>> M1, vector<vector<long>> M2, int m, int l, int n, long p)
{
  vector<vector<long>> M(m, vector<long>(n));
  // long p = INT32_MAX;
  std::ofstream fout;
  fout.open("M1*M2_plaintext.txt", ios::out);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      M[i][j] = 0;
      for (int k = 0; k < l; k++)
      {
        M[i][j] = M[i][j] + M1[i][k] * M2[k][j];
      }
      long a = M[i][j] % p;
      if (abs(a) > p / 2)
      {
        a = (a > 0 ? a - p : a + p);
      }
      std::cout << setw(log10(p) + 2) << a << " ";
      fout << setw(log10(p) + 2) << a << " ";
    }
    cout << endl;
    fout << endl;
  }
  cout << "See M1*M2_plaintext.txt" << endl;
  fout.close();
  cout << endl;
}

void plain_text_multiplication(vector<vector<long>> M1, vector<vector<long>> M2, int m, int l, int n, long p, int s)
{
  vector<vector<long>> M(m, vector<long>(n));
  vector<vector<long>> M_copy(m, vector<long>(n));
  // long p = INT32_MAX;
  std::ofstream fout;
  fout.open("M1*M2_plaintext.txt", ios::out);
  M = M1;
  for (int k = 0; k < s; k++)
  {
    M_copy = M;
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        M[i][j] = 0;
        for (int k = 0; k < l; k++)
        {
          M[i][j] = M[i][j] + M_copy[i][k] * M2[k][j];
        }
      }
    }
  }
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      long a = M[i][j] % p;
      if (abs(a) > p / 2)
      {
        a = (a > 0 ? a - p : a + p);
      }
      std::cout << setw(log10(p) + 2) << a << " ";
      fout << setw(log10(p) + 2) << a << " ";
    }
    cout << endl;
    fout << endl;
  }
  cout << "See M1*M2_plaintext.txt" << endl;
  fout.close();
  cout << endl;
}

void printf_matrix(const std::vector<std::vector<long>> &M, int m, int n, int w, long &p)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      std::cout << setw(w) << M[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

void printf_matrix(const vector<vector<long>> &M, int m, int n, int w, long &p, ofstream &fout)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int a;
      if (M[i][j] % p < p / 2)
        a = M[i][j] % p;
      else
        a = M[i][j] % p - p;
      fout << setw(w) << a << " ";
    }
    fout << std::endl;
  }
}

void encode(Ptxt_new &Ptxt, std::vector<std::vector<long>> M, int t, int m, int l, int p1, int p2, int p3)
{
  if (t == 0)
  {
    for (int i = 0; i < m; i++)
    {
      for (int k = 0; k < l; k++)
      {
        for (int j = 0; j < p1; j++)
        {
          Ptxt[k * p1 * p3 + p1 * i + j] = M[i][k];
        }
      }
    }
  }
  if (t==1)
  {
    for (int i = 0; i < m; i++)
    {
      for (int k = 0; k < l; k++)
      {
        for (int j = 0; j < p3; j++)
        {
          Ptxt[i * p1 * p3 + p1 * j + k] = M[i][k];
        }
      }
    }
  }
  if (t == 2) 
  {
    for (int i = 0; i < l; i++)
    {
      for (int k = 0; k < m; k++)
      {
        for (int j = 0; j < p3; j++)
        {
          Ptxt[i * p1 * p3 + p1 * j + k] = M[k][i];
        }
      }
    }
  }
}

std::vector<std::vector<long>> randmat(int m, int n, int q, long p, std::string s)
{
  std::vector<std::vector<long>> M(m, std::vector<long>(n, 0));
  int d = int(log10(q) + 2);
  srand((unsigned)time(NULL) + long(&M));
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      long a = rand() % q;
      if (a > q / 2)
      {
        a = a - q;
      }
      M[i][j] = a;
    }
  }
  if (n <= 64)
  {
    printf_matrix(M, m, n, d, p);
  }
  std::ofstream fout;
  fout.open(s, ios::out);
  printf_matrix(M, m, n, d, p, fout);
  cout << "See " << s << endl;
  fout.close();
  return M;
}

helib::Ctxt FastTrace(helib::Ctxt ctxt, long r, std::vector<long> tr)
{
  if (r == 1)
  {
    return ctxt;
  }
  if (r % 2 == 0)
  {
    helib::Ctxt temp = FastTrace(ctxt, r / 2, tr);
    helib::Ctxt temp2 = temp;
    temp2.smartAutomorph(tr[r / 2]);
    temp.addCtxt(temp2);
    return temp;
  }
  else
  {
    helib::Ctxt temp = FastTrace(ctxt, (r - 1) / 2, tr);
    helib::Ctxt temp2 = temp;
    helib::Ctxt temp3 = ctxt;
    temp2.smartAutomorph(tr[(r - 1) / 2]);
    temp3.smartAutomorph(tr[r - 1]);
    temp.addCtxt(temp2);
    temp.addCtxt(temp3);
    return temp;
  }
}

void GenTrace(helib::SecKey &secret_key, long r, std::vector<long> tr)
{
  if (r == 1)
  {
    return;
  }
  if (r % 2 == 0)
  {

    GenTrace(secret_key, r / 2, tr);
    secret_key.GenKeySWmatrix(1, tr[r / 2], 0, 0);
    secret_key.setKeySwitchMap();
    return;
  }
  else
  {
    GenTrace(secret_key, (r - 1) / 2, tr);
    secret_key.GenKeySWmatrix(1, tr[(r - 1) / 2], 0, 0);
    secret_key.GenKeySWmatrix(1, tr[r - 1], 0, 0);
    secret_key.setKeySwitchMap();
    return;
  }
}

std::vector<std::vector<long>> randmat_sym(int n, int q, long p, std::string s)
{
  std::vector<std::vector<long>> M(n, std::vector<long>(n, 0));
  int d = int(log10(q) + 2);
  srand((unsigned)time(NULL) + long(&M));
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      long a = rand() % q;
      if (a > q / 2)
      {
        a = a - q;
      }
      M[i][j] = a;
      M[j][i] = a;
    }
  }
  if (n <= 64)
  {
    printf_matrix(M, n, n, d, p);
  }
  std::ofstream fout;
  fout.open(s, ios::out);
  printf_matrix(M, n, n, d, p, fout);
  cout << "See " << s << endl;
  fout.close();
  return M;
}
// for (int i = 0; i < tr2.size(); i++)
// {
//   helib::Ctxt temp = ctxt1;
//   temp.smartAutomorph(tr2[i]);
//   ctxt_result.addCtxt(temp);
// }

// for (int k = 0; k < p2; k++)
// {
//   for (int i = 0; i < p1 * p3; i++)
//   {
//     Ptxt[p1 * p3 * k + i] = i + k + 1;
//   }
// }