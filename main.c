#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 1000			// # of excitatory and inhibitory neurons
#define NN ((N)*(N))
#define N2 ((N)*(2))
#define T 3000			// # of steps
#define TAU_EX	50.0		// time constant for excitatory neurons
#define TAU_INH	70.0		// time constant for inhibitory neurons
#define THETA_EX	0.1	// threshold for excitatory neurons
#define THETA_INH	0.1	// threshold for inhibitory neurons
// UNUSED #define C_EXEX		2.0	// ex <- ex connection weights 
#define C_INHEX		4.0	// inh <- ex connection weights 
#define C_EXINH		16.0	// ex <- inh connection weights 
#define C_INHINH	6.0	// inh <- inh connection weights 
#define ETA		0.5	// noise level
#define ALPHA		0.005	// learning coefficient of ex <- ex weights
#define TAU		10.0	// decay of ex <- ex weights
#define I_CS	1.0		// intensity of input signals
#define T_CS	50.0		// duration of input signals
#define N_CS	200		// # of neurons receiving input signals

#define EX(x)	(x)
#define INH(x)	((N)+(x))
#define id(i,j) (((j)+(N)*(i)))
#define did(t,i) (((i)+(N2)*(t)))

double *w_exex, *w_inhex, *w_exinh, *w_inhinh; // connection weights
double *data;			// all neurons' activities

extern void init_genrand(unsigned long);
extern double genrand_real2(void);

void initialize(void)
{
  w_exex = (double *)malloc(NN * sizeof(double));
  w_inhex = (double *)malloc(NN * sizeof(double));
  w_exinh = (double *)malloc(NN * sizeof(double));
  w_inhinh = (double *)malloc(NN * sizeof(double));
  data = (double *)malloc(T * N2 * sizeof(double));
}

void finalize(void)
{
  free(w_exex);
  free(w_inhex);
  free(w_exinh);
  free(w_inhinh);
  free(data);
}

void set_connections(const unsigned long seed, const char *infile)
{
  int i, j;
  FILE *file;
  char buf[1024];

  // set ex <- ex connections by reading infile
  file = fopen(infile, "r");
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      fgets(buf, 1024, file);
      w_exex[id(i,j)] = atof(buf);
    }
  }
  fclose(file);

  // set inh <- ex connections
  init_genrand(seed);
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      if (genrand_real2() < 0.5){
	w_inhex[id(i,j)] = 2*C_INHEX/(double)N;
      }else{
	w_inhex[id(i,j)] = 0;
      }
    }
  }

  // set ex <- inh connections
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      w_exinh[id(i,j)] = 0;
    }
  }
  for(i = 0; i < N; i++){
    w_exinh[id(i,i)] = C_EXINH;
  }

  // set inh <- inh connections
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      w_inhinh[id(i,j)] = C_INHINH/(double)N;
    }
  }
}

void dudt(double du[], const double u[], const double z[],
	  const int t)
{
  int i, j;
  double ex[N2], inh[N2], I[N];

  // compute PSPs
  for(i = 0; i < N2; i++){
    ex[i] = inh[i] = 0;
  }
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      ex[EX(i)] += w_exex[id(i,j)]*z[EX(j)];
    }
  }
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      inh[EX(i)] += w_exinh[id(i,j)]*z[INH(j)];
    }
  }
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      ex[INH(i)] += w_inhex[id(i,j)]*z[EX(j)];
    }
  }
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      inh[INH(i)] += w_inhinh[id(i,j)]*z[INH(j)];
    }
  }

  // set input signals
  if (t < T_CS){
    for(i = 0; i < N_CS; i++){
      I[i] = I_CS*(1+ETA*2*(1-genrand_real2()));
    }
    for(i = N_CS; i < N; i++){
      I[i] = I_CS*ETA*genrand_real2();
    }
  }else{
    for(i = 0; i < N; i++){
      I[i] = I_CS*ETA*genrand_real2();
    }
  }

  // compute du
  for(i = 0; i < N; i++){
    du[EX(i)] = (1.0/TAU_EX)*(-u[EX(i)] + I[EX(i)] + ex[EX(i)] - inh[EX(i)]);
  }
  for(i = 0; i < N; i++){
    du[INH(i)] = (1.0/TAU_INH)*(-u[INH(i)] + ex[INH(i)] - inh[INH(i)]);
  }
  
}

void exec(const unsigned long seed)
{
  int t, i;
  double u[N2], du[N2], z[N2];

  for(i = 0; i < N2; i++){
    u[i] = z[i] = 0;
  }

  init_genrand(seed);
  for(t = 0; t < T; t++){

    dudt(du, u, z, t);

    for(i = 0; i < N2; i++){
      u[i] += du[i];
    }

    for(i = 0; i < N; i++){
      if (u[EX(i)] > THETA_EX){
	z[EX(i)] = u[EX(i)];
      }else{
	z[EX(i)] = 0;
      }
      data[did(t,EX(i))] = z[EX(i)];
    }
    for(i = 0; i < N; i++){
      if (u[INH(i)] > THETA_INH){
	z[INH(i)] = u[INH(i)];
      }else{
	z[INH(i)] = 0;
      }
      data[did(t,INH(i))] = z[INH(i)];
    }
  }
}

void update_connections(const char *outfile)
{
  int t, i, j;
  double *w_new, norm[N], r;
  FILE *file;

  w_new = (double *)malloc(NN * sizeof(double));

  for(i = 0; i < N; i++){
    norm[i] = 0;
    for(t = 0; t < T; t++){
      norm[i] += data[did(t,EX(i))]*data[did(t,EX(i))];
    }
    norm[i] = sqrt(norm[i]);
  }

  for(i = 0; i < N; i++){
    for(j = i+1; j < N; j++){
      r = 0;
      for(t = 0; t < T; t++){
	r += data[did(t,EX(i))]*data[did(t,EX(j))];
      }
      if (norm[i] != 0 && norm[j] != 0){
	w_new[id(i,j)] = r/(norm[i]*norm[j]);
      }else{
	w_new[id(i,j)] = 0;
      }
      w_new[id(j,i)] = w_new[id(i,j)];
    }
  }

  for(i = 0; i < N; i++){
    w_new[id(i,i)] = 0;
  }

  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      w_new[id(i,j)] = ((1-1.0/TAU)*w_exex[id(i,j)]
			+(ALPHA/TAU)*w_new[id(i,j)]);
    }
  }

  file = fopen(outfile, "w");
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      fprintf(file, "%f\n", w_new[id(i,j)]);
    }
  }
  fclose(file);

  free(w_new);
}

void output_data(const char *outprefix)
{
  FILE *file;
  char buf[1024];
  int t, i;

  sprintf(buf, "%s.r", outprefix);
  file = fopen(buf, "w");
  for(t = 0; t < T; t++){
    for(i = 0; i < N; i++){
      if (data[did(t,EX(i))] >= THETA_EX){
	fprintf(file, "%d %d\n", t, i);
      }
    }
  }
  fclose(file);

  sprintf(buf, "%s.a", outprefix);
  file = fopen(buf, "w");

  for(t = 0; t < T; t++){
    for(i = 0; i < N; i++){
      fprintf(file, "%f\n", data[did(t,EX(i))]);
    }
  }

  fclose(file);
}

int main(int argc, char *argv[])
{
  if (argc < 6){
    fprintf(stderr, "usage: %s <seed> <output> <w.in> <w.out> <noiseseed>\n",
	    argv[0]);
    exit(1);
  }

  initialize();

  set_connections(atol(argv[1]), argv[3]);

  exec(atol(argv[5]));

  update_connections(argv[4]);

  output_data(argv[2]);

  finalize();

  return 0;
}
