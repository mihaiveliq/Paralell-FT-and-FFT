#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	double *vec_in;
	int n;
	int num_thrds;
} Var;

//Useful global var
double *re_out, *im_out;
Var var;

int min(int a, int b) {
	return a > b ? b : a;
}

//FT par
void *par_FT(void *aux_param) {

	int thrd_id = *(int*) aux_param;

    //Extract values
	int num_thrds = var.num_thrds;
	int n = var.n;

	//Preparing the calculation
	double angle;
	int t, k;
	const double pi = 3.14159265359;
	int start = thrd_id * (n / num_thrds);
	int end = (thrd_id == num_thrds - 1) ? n : ((thrd_id + 1) * (n / num_thrds));

	//The calculation
    for (k = start; k < end; ++k) {
	    for (t = 0; t < n; ++t) {
			angle = 2 * pi * t * k / n;
			re_out[k] += var.vec_in[t] * cos(angle);
			im_out[k] += var.vec_in[t] * sin(angle) * (-1);
	    }
    }

	return NULL;
}

//FT seq
void FT(double *vec_in, int n, char *name_out) {

	//Opening output file
    FILE *out_file;
	out_file = fopen(name_out, "wt");
  	if (out_file == NULL) {
    	fprintf(stdout, "Failed to open the output file\n.");
    	exit(1);
  	}

	//Printing the number of values in freq domain
    fprintf(out_file, "%d\n", n);

	//Preparing the computation
	double re_out, im_out, angle;
	int t, k;
	const double pi = 3.14159265359;

	//The computation
	for (k = 0; k < n; ++k) {
		re_out = 0;
		im_out = 0;
		for (t = 0; t < n; ++t) {
			angle = 2 * pi * t * k / n;
			re_out += vec_in[t] * cos(angle);
			im_out += vec_in[t] * sin(angle) * (-1);
		}
		//Printing the values in freq domain
		fprintf(out_file, "%f %f\n", re_out, im_out);
	}

	//Closing the output file
    fclose(out_file);
}

int main(int argc, char * argv[]) {

	if (argc < 3) {
    	fprintf(stdout, "Usage: %s <in_file> <out_file> <numThreads>\n", argv[0]);
    	exit(1);
  	}

	//Opening input file
	FILE *in_file;
	in_file = fopen(argv[1], "rt");
  	if (in_file == NULL) {
    	fprintf(stdout, "Failed to open the input file\n.");
    	exit(1);
  	}

	int n, i;
	int num_thrds = atoi(argv[3]);

	//Reading the number of values in the time domain
	fscanf(in_file, "%d", &n);

	//Building thread function's vars
	var.vec_in = malloc(n * sizeof(double));
	var.n = n;
	var.num_thrds = num_thrds;
	pthread_t tid[num_thrds];

	//Reading the values in the time domain
	for(i = 0; i < n; ++i) {
		fscanf(in_file, "%lf", &var.vec_in[i]);
	}

	//Closing input file
	fclose(in_file);

	//Usefull tools
	re_out = malloc(n * sizeof(double));
	im_out = malloc(n * sizeof(double));
	int ids[num_thrds];

	for (i = 0; i < num_thrds; ++i) {
		ids[i] = i;
	}

	//Creating threads   
	for (i = 0; i < num_thrds; ++i) {
		pthread_create(&(tid[i]), NULL, par_FT, &(ids[i]));
	}

	//Joining them
	for (i = 0; i < num_thrds; ++i) {
		pthread_join(tid[i], NULL);
	}

	//Opening output file
    FILE *out_file;
	out_file = fopen(argv[2], "wt");
  	if (out_file == NULL) {
    	fprintf(stdout, "Failed to open the output file\n.");
    	exit(1);
  	}

	//Printing the number of values in freq domain
    fprintf(out_file, "%d\n", n);

	for (i = 0; i < n; ++i) {
		//Printing the values in freq domain
	    fprintf(out_file, "%f %f\n", re_out[i], im_out[i]);
	}

	//Closing the output file
    fclose(out_file);

	free(re_out);
	free(im_out);

	/*Seq FT*/
	//Print values in freq domain
	//FT(var.vec_in, n, argv[2]);

    free(var.vec_in);

	return 0;
}
