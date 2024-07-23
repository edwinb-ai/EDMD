#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#define MAXN 5000
#define MAXBATCH 10
#define MAXFILES 600

#define tprintf \
  if (test)     \
  printf
#include "calcdiff.h"
// Event driven MD.

// particle particles[MAXBATCH][MAXN];
// particle* celllist[CEL][CEL][CEL];
int npart;
double xsize = 0.0;
double ysize = 0.0;
double zsize = 0.0;

double delta;
char inputfile[1024];
char outputfile[1024];
char buffer[1024];

double r2s[MAXBATCH][MAXFILES];
double times[MAXFILES];
double starttime = -1;

double ts[MAXFILES];
int bs[MAXFILES];
int ss[MAXFILES];
int maxb;
int maxs;
double maxtime;
double tbatch = 0;

double x[MAXN], y[MAXN], z[MAXN];

int totbonds = 0;

double startconf[MAXBATCH][MAXN][3];

char *getlasttoken(char *str, char *del)
{
  char *point1;
  char *point2;
  point1 = strtok(str, del);
  point2 = point1;
  while (point2)
  {
    point1 = point2;
    point2 = strtok(NULL, del);
  }
  return point1;
}

char *movepastnextdot(char *p)
{
  while (p[0] != '\0' && p[0] != '.')
  {
    p++;
  }
  if (p[0] == '\0')
    return NULL;
  p++;
  return p;
}

void mygetline(char *str, FILE *f)
{
  int comment = 1;
  while (comment)
  {
    if (!fgets(str, 1024, f))
      return;
    if (str[0] != '#')
      comment = 0;
  }
}

double gettime(char* filename)
{
  FILE* f = fopen(filename, "r");
  int dummy; 
  double dummyd; 
  double t;
  int r = fscanf(f, "%d\n%lf %lf %lf Time: %lf", &dummy, &dummyd, &dummyd, &dummyd, &t);
  if (r != 5) 
  {
    printf ("Failed to read time from snapshot!\n");
    exit(3);
  }
  return t;  
}

int main(int argc, char *argv[])
{
  FILE *f;
  if (argc < 2)
  {
    printf("No filename given!\n");
    exit(3);
  }
  int i, j, k, b, s, tmp;
  double t;
  int numfiles = argc - 1;

  for (i = 0; i < MAXBATCH; i++)
  {
    for (j = 0; j < MAXFILES; j++)
    {
      r2s[i][j] = 0;
    }
  }

  for (j = 0; j < MAXFILES; j++) {
    times[j] = 0.0;
  }

  for (i = 0; i < numfiles; i++)
  {
    strcpy(inputfile, argv[i + 1]);
    char *point = getlasttoken(argv[i + 1], "/");
    point = movepastnextdot(point);
    while (point)
    {
      switch (point[0])
      {
      case 'b':
        point++;
        sscanf(point, "%d", &b);
        break;
      case 's':
        point++;
        sscanf(point, "%d", &s);
        break;
      }
      point = movepastnextdot(point);
    }
    t = gettime(inputfile);
    ts[i] = t;
    bs[i] = b;
    ss[i] = s;
    // Set starttime once
    if (starttime < 0) {
        starttime = t;
    }
    t -= starttime;
    printf("Starttime: %lf\n", starttime);
    printf("t: %.16lf, b: %d, s: %d\n", t, b, s);

    if (s == 0) //First snapshot in batch: store the starting configuration
    {
      if (b == 1){
        tbatch = t;
      }
      f = fopen(inputfile, "r");
      tmp = fscanf(f, "%d\n", &npart);
      if (tmp != 1)
      {
        printf("error reading npart %d\n", npart);
        exit(3);
      }
      mygetline(buffer, f);

      tmp = sscanf(buffer, "%lf %lf %lf\n", &xsize, &ysize, &zsize);
      if (tmp != 3)
      {
        printf("error reading box\n");
        exit(3);
      }

      for (j = 0; j < npart; j++)
      {
        mygetline(buffer, f);
        char letter;
        int nsides;
        double rad;
        tmp = sscanf(buffer, "%c %lf %lf %lf %lf",
          &letter, &(startconf[b][j][0]), &(startconf[b][j][1]),
          &(startconf[b][j][2]), &rad);
        if (tmp < 5)
        {
          printf("can't read from string: %s\n", buffer);
          exit(3);
        }
      }
      fclose(f);
    }
    if (b > maxb) {
      maxb = b;
    }

    if (s > maxs) {
      maxs = s;
    }

    double r2 = 0.0;

    f = fopen(inputfile, "r");
    tmp = fscanf(f, "%d\n", &npart);
    if (tmp != 1)
    {
      printf("error reading npart %d\n", npart);
      exit(3);
    }
    
    mygetline(buffer, f);
    tmp = sscanf(buffer, "%lf %lf %lf", &xsize, &ysize, &zsize);
    
    if (t > maxtime) {
      maxtime = t;
    }

    if (b == 0) {
      times[s] = t;
      // printf("t: %.16lf, b: %d, s: %d\n", t, b, s);
      // printf("%lf %d\n", times[s], s);
    // printf("%lf %lf\n", times[s], t);
    }
    if (tmp != 3)
    {
      printf("error reading box \n");
      exit(3);
    }

    for (j = 0; j < npart; j++)
    {
      mygetline(buffer, f);
      char letter;
      double rad;
      tmp = sscanf(buffer, "%c %lf %lf %lf %lf",
        &letter, &(x[j]), &(y[j]), &(z[j]), &rad);
      if (tmp < 5)
      {
        printf("can't read from string: %s\n", buffer);
        exit(3);
      }
    }
    fclose(f);

    for (j = 0; j < npart; j++)
    {
      x[j] = x[j] - startconf[b][j][0];
      y[j] = y[j] - startconf[b][j][1];
      z[j] = z[j] - startconf[b][j][2];
      r2 = x[j] * x[j] + y[j] * y[j] + z[j] * z[j];
      r2s[b][s] += r2;
    }
  }
  for (k = 0; k < b; k++)
  {
    times[maxs + k] = k * tbatch;
  }
  printf("%d batches, maxtime = %lf, maxs: %d, tbatch: %lf\n",
    maxb, maxtime, maxs, tbatch);

  strcpy(outputfile, "r2_clang.dat");

  printf("Output: %s\n", outputfile);

  double totalr2 = 0.0;
  f = fopen(outputfile, "w");
  for (i = 0; i < maxs; i++)
  {
    totalr2 = 0.0;
    fprintf(f, "%.12lf ", times[i]);
    for (j = 0; j < maxb; j++)
    {
      r2s[j][i] /= npart;
      totalr2 += r2s[j][i];
      // if (j == 0 || r2s[j][i] > 0)
      //   fprintf(f, "%.12lf ", r2s[j][i]);
    }
    fprintf(f, "%.12lf\n", totalr2 / maxb);
  }
  fclose(f);

  return 0;
}
