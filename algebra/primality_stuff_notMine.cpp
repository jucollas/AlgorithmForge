/*
Accepted version of the solution to the Factoring Large Integers problem
Use Factorization with Pollards Rho + Primality de Miller-Rabin
Note: The IO part was cut down by edition

Author: Juan Diego Collazos
Date : 20/5/2025
Problem : https://onlinejudge.org/external/114/11476.pdf

Reference:
Pollards Rho Algorithm for Prime Factorization: https://www.geeksforgeeks.org/pollards-rho-algorithm-prime-factorization/
Primality Test | Set 3 (Miller-Rabin) : https://www.geeksforgeeks.org/primality-test-set-3-miller-rabin/
*/

#include <bits/stdc++.h>

using namespace std;
using ll = long long;

mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

ll gcd(ll a, ll b) {
  while (b != 0) {
    ll t = b;
    b = a % b;
    a = t;
  }
  return a;
}

ll mulmod(ll a, ll b, ll mod) {
  ll ans = 0;
  a %= mod;
  while (b > 0) {
    if (b & 1)
      ans = (ans + a) % mod;
    a = (a * 2) % mod;
    b >>= 1;
  }
  return ans;
}

ll powmod(ll base, ll exp, ll mod) {
  ll ans = 1;
  base %= mod;
  while (exp > 0) {
    if (exp & 1)
      ans = mulmod(ans, base, mod);
    base = mulmod(base, base, mod);
    exp >>= 1;
  }
  return ans;
}

bool is_prime(ll n, int k = 5){
  bool ans;
  if (n < 2)
    ans = false;
  else if (n == 2 || n == 3)
    ans = true;
  else if (n % 2  == 0)
    ans = false;
  else{
    ll r = 0, d = n - 1;
    while (d % 2 == 0){
      d >>= 1;
      r++;
    }
    ans = true;
    int i = 0;
    while(i < k && ans){
      ll a = rng() % (n - 3) + 2;
      ll x = powmod(a, d, n);
      if (x != 1 && x != n -1){
        ll j = 0;
        bool composite = true;
        while(j < r - 1 && composite){
          x = powmod(x, 2, n);
          if (x == n - 1)
            composite = false;
          ++j;
        }
        if (composite)
          ans = false;
      }
      ++i;
    }
  }
  return ans;
}

ll f(ll x, ll n, ll c){
  return (mulmod(x, x, n) + c) % n;
}

ll pollards_rho(ll n) {
  ll ans;
  if (n % 2 == 0)
    ans = 2;
  else if(is_prime(n))
    ans = n;
  else{
    ans = -1;
    while(ans == -1){
      ll c = rng() % n;
      ll x = 2, y = 2, d = 1;
      while(d == 1) {
        x = f(x, n, c);
        y = f(f(y, n, c), n, c);
        d = gcd(abs(x - y), n);
      }
      if (d != n)
        ans = d;
    }
  }
  return ans;
}

void factor(ll n, vector<ll>& factors){
  if (n != 1){
    if (is_prime(n)){
      factors.push_back(n);
    }else{
      ll k = pollards_rho(n);
      factor(k, factors);
      factor(n / k, factors);
    }
  }
}