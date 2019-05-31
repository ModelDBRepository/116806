#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gd.h>

#define N 1000
#define T 1000
#define T2 3000

double **data1, **data2, **c;
double avg[T], std[T];

void data_initialize(void)
{
  int i;
  data1 = (double **)malloc(T2*sizeof(double *));
  for(i = 0; i < T2; i++){
    data1[i] = (double *)malloc(N*sizeof(double));
  }
  data2 = (double **)malloc(T2*sizeof(double *));
  for(i = 0; i < T2; i++){
    data2[i] = (double *)malloc(N*sizeof(double));
  }
}
void data_finalize(void)
{
  int i;
  for(i = 0; i < T2; i++){
    free(data1[i]);
  }
  free(data1);
  for(i = 0; i < T2; i++){
    free(data2[i]);
  }
  free(data2);
}
void c_initialize(void)
{
  int i;
  c = (double **)malloc(T2*sizeof(double *));
  for(i = 0; i < T2; i++){
    c[i] = (double *)malloc(T2*sizeof(double));
  }
}
void c_finalize(void)
{
  int i;
  for(i = 0; i < T2; i++){
    free(c[i]);
  }
  free(c);
}
void input(char *infile1, char *infile2)
{
  FILE *file;
  int t, i;
  char filename[1024], buf[1024];

  sprintf(filename, "%s", infile1);
  file = fopen(filename, "r");
  for(t = 0; t < T2; t++){
    for(i = 0; i < N; i++){
      fgets(buf, 1024, file);
      data1[t][i] = atof(buf);
    }
  }
  fclose(file);

  sprintf(filename, "%s", infile2);
  file = fopen(filename, "r");
  for(t = 0; t < T2; t++){
    for(i = 0; i < N; i++){
      fgets(buf, 1024, file);
      data2[t][i] = atof(buf);
    }
  }
  fclose(file);
}
void similarity_index(void)
{
  int t1, t2, t, i;
  double norm1[T2], norm2[T2], r;

  for(t = 0; t < T2; t++){
    r = 0;
    for(i = 0; i < N; i++){
      r += data1[t][i]*data1[t][i];
    }
    norm1[t] = sqrt(r);
  }
  for(t = 0; t < T2; t++){
    r = 0;
    for(i = 0; i < N; i++){
      r += data2[t][i]*data2[t][i];
    }
    norm2[t] = sqrt(r);
  }

  for(t1 = 0; t1 < T2; t1++){
    for(t2 = 0; t2 < T2; t2++){
      r = 0;
      for(i = 0; i < N; i++){
	r += data1[t1][i]*data2[t2][i];
      }
      if (norm1[t1]*norm2[t2] == 0){
	c[t1][t2] = 0;
      }else{
	c[t1][t2] = r/(norm1[t1]*norm2[t2]);
      }
    }
  }
}
void output(char *outprefix)
{
  FILE *file;
  int t1, t2, i;
  char filename[1024];
  gdImagePtr im;
  int gray[256], idx;
  int dt;
  double avg, std;

  im= gdImageCreate(T, T);
  for(i = 0; i < 256; i++){
    gray[i] = gdImageColorAllocate(im, i, i, i);
  }
  for(t1 = 0; t1 < T; t1++){
    for(t2 = 0; t2 < T; t2++){
      idx = floor(255*c[t1][t2]);
      gdImageSetPixel(im, t1, t2, gray[idx]);
    }
  }

  // PNG file of Similarity matrix
  sprintf(filename, "%s.png", outprefix);
  file = fopen(filename, "wb");
  gdImagePng(im, file);
  fclose(file);
  gdImageDestroy(im);

  /*
  // Rows of Similarity matrix
  for(t1 = 0; t1 < T2; t1+=100){
    sprintf(filename, "%s.%03d", outprefix, t1);
    file = fopen(filename, "w");
    for(t2 = 0; t2 < T2; t2++){
      fprintf(file, "%d %f\n", t2, c[t1][t2]);
    }
    fclose(file);
  }
  */

  sprintf(filename, "%s.si", outprefix);
  file = fopen(filename, "w");
  for(dt = -T/2; dt < T/2; dt++){
    avg = 0;
    for(t1 = T+T/2; t1 < T2-T/2; t1++){
      t2 = t1 + dt;
      avg += c[t1][t2]/T;
    }
    std = 0;
    for(t1 = T+T/2; t1 < T2-T/2; t1++){
      t2 = t1 + dt;
      std += (c[t1][t2]-avg)*(c[t1][t2]-avg)/T;
    }
    std = sqrt(std);
    fprintf(file, "%d %f %f %f %f\n", dt, avg, std, avg+std, avg-std);
  }
  fclose(file);


  sprintf(filename, "%s.d", outprefix);
  file = fopen(filename, "w");
  for(dt = 0; dt < T2; dt++){
    fprintf(file, "%d %f\n", dt, c[dt][dt]);
  }
  fclose(file);

}
int main(int argc, char *argv[])
{
  if (argc != 4){
    fprintf(stderr, "usage: %s <infile1> <infile2> <outprefix>\n", argv[0]);
    exit(0);
  }

  data_initialize();
  c_initialize();

  input(argv[1], argv[2]);

  similarity_index();

  output(argv[3]);

  data_finalize();
  c_finalize();

  return 0;
}
