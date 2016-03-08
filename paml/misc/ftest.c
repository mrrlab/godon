#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#define N 100

typedef struct Btest {
  int n;
  double a;
  double b;
  bool median; 
} Btest;

Btest btests[] = {
  // n, a, b, median
  {4, 0.5, 10, false},
  {4, 0.5, 10, true},
  {8, 2, .1, false},
  {7, 15, 1, true},
  {4, 1.16, 3.54, false},
  {4, 1.16, 3.54, true},
};

void parr(double *a, int n) {
  printf("{");
  for (int i=0; i < n; i++) {
    printf("%f", a[i]);
    if (i != n-1) {
      printf(",");
    }
  }
  printf("}\n");
}

int main() {
  double eps = 0, alneps = 0, sml = 0, alnsml = 0;
  eps = pow((double)FLT_RADIX, -(double)DBL_MANT_DIG);
  alneps = log(eps);
  sml = DBL_MIN;
  alnsml = log(sml);

  printf("%g,%g,%g,%g\n", eps, alneps, sml, alnsml);
  double rates[N], x[N];
  printf("...beta...\n");
  for (int i=0; i < sizeof(btests)/sizeof(Btest); i++) {
    Btest t = btests[i];
    DiscreteBeta(rates, x, t.a, t.b, t.n, t.median);
    parr(rates, t.n);
    parr(x, t.n);
  }
  printf("...gamma...\n");
  for (int i=0; i < sizeof(btests)/sizeof(Btest); i++) {
    Btest t = btests[i];
    DiscreteGamma(rates, x, t.a, t.b, t.n, t.median);
    parr(rates, t.n);
    parr(x, t.n);
  }


  return 0;
  
}
