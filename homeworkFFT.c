#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
 
double PI;
typedef double complex cplx;

int n;
cplx *buf_t;
cplx *out_t;
int num_thrds;
pthread_barrier_t bar;
pthread_mutex_t mut;

void _fft(cplx *buf, cplx *out, int n, int step) {
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2] = out[i] + t;
			buf[(i + n) / 2] = out[i] - t;
		}
	}
}

void *thrd_fft(void *var) {
	int thrd_id = *(int*) var;
	int step;

	switch(num_thrds) {
		case 1:
			/*For root
				R
			*/
			_fft(buf_t, out_t, n, 1);
			
			break;

		case 2:
			/*For first layer
				 R
			 	/ \
			  C'1 C'2
			
			From _fft(buf_t, out_t, n, step = 1) =>
			=> 2 rec: _fft(out_t + 0, buf_t + 0, n, step = 1 * 2);
					  _fft(out_t + 1, buf_t + 1, n, step = 1 * 2);
			*/
			step = 1;
			_fft(out_t + thrd_id, buf_t + thrd_id, n, step * 2);

			//For R
			pthread_barrier_wait(&bar);
			if (thrd_id == 0) {
				for (int i = 0; i < n; i += 2 * step) {
					cplx t = cexp(-I * PI * i / n) * out_t[i + step];
					buf_t[i / 2] = out_t[i] + t;
					buf_t[(i + n) / 2] = out_t[i] - t;
				}
			}

			break;

		case 4:
		/*For second layer
				  _R_
			 	/     \
			  C'1     C'2
			  / \    /  \
		   C"1 C"2  C"3 C"4
			
			From _fft(buf_t, out_t, n, step = 1) =>
								Step = 1
			=> 2 rec: _fft(out_t + 0, buf_t + 0, n, step = 1 * 2)| =>
					  _fft(out_t + 1, buf_t + 1, n, step = 1 * 2)|
								Step = 2
			=> 4 rec: _fft(buf_t + 0 + 0, out_t + 0 + 0, n, step = 2 + 2 * 2)
					  _fft(buf_t + 0 + 2, out_t + 0 + 2, n, step = 2 + 2 * 2) 
					  _fft(buf_t + 1 + 0, out_t + 1 + 0, n, step = 2 + 2 * 2)
					  _fft(buf_t + 1 + 2, out_t + 1 + 2, n, step = 2 + 2 * 2) 
			*/
			step = 2;
			_fft(buf_t + thrd_id, out_t + thrd_id, n, step * 2);

			pthread_barrier_wait(&bar);

			//For C'1
			if (thrd_id == 0) {
				for (int i = 0; i < n; i += 2 * step) {
						cplx t = cexp(-I * PI * i / n) * buf_t[i + step];
						out_t[i / 2] = buf_t[i] + t;
						out_t[(i + n) / 2] = buf_t[i] - t;
				}
			}

			//For C'2
			if (thrd_id == 1) {
				for (int i = 0; i < n; i += 2 * step) {
						cplx t = cexp(-I * PI * i / n) * (buf_t + 1)[i + step];
						(out_t + 1)[i / 2] = (buf_t + 1)[i] + t;
						(out_t + 1)[(i + n) / 2] = (buf_t + 1)[i] - t;
				}
			}

			pthread_barrier_wait(&bar);

			//For R
			if (thrd_id == 2) {
				step = 1;

				for (int i = 0; i < n; i += 2 * step) {
						cplx t = cexp(-I * PI * i / n) * out_t[i + step];
						buf_t[i / 2] = out_t[i] + t;
						buf_t[(i + n) / 2] = out_t[i] - t;
				}
			}
	}
	return NULL;
}
 
void fft()
{
	out_t = malloc(n * sizeof(cplx));

	for (int i = 0; i < n; i++) {
		out_t[i] = buf_t[i];
	}

	pthread_t tid[num_thrds];
	int ids[num_thrds];
	pthread_barrier_init(&bar, NULL, num_thrds);
	pthread_mutex_init(&mut, NULL);

	for (int i = 0; i < num_thrds; ++i) {
		ids[i] = i;
	}

	//Creating threads   
	for (int i = 0; i < num_thrds; ++i) {
		pthread_create(&(tid[i]), NULL, thrd_fft, &(ids[i]));
	}

	//Joining them
	for (int i = 0; i < num_thrds; ++i) {
		pthread_join(tid[i], NULL);
	}

	pthread_barrier_destroy(&bar);
	pthread_mutex_destroy(&mut);

	//_fft(buf_t, out_t, n, 1);
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

	num_thrds = atoi(argv[3]);

	//Reading the number of values in the time domain
	fscanf(in_file, "%d", &n);
	
	PI = atan2(1, 1) * 4;
	double a;
	buf_t = malloc(n * sizeof(cplx));

	//Reading the values in the time domain
	for(int i = 0; i < n; ++i) {
		fscanf(in_file, "%lf", &a);
		buf_t[i] = CMPLX(a, 0.0);
	}

	//Closing input file
	fclose(in_file);
	
	//Compute
	fft();

	//Opening output file
    FILE *out_file;
	out_file = fopen(argv[2], "wt");
  	if (out_file == NULL) {
    	fprintf(stdout, "Failed to open the output file\n.");
    	exit(1);
  	}

	//Printing the number of values in freq domain
    fprintf(out_file, "%d\n", n);

	for (int i = 0; i < n; ++i) {
		//Printing the values in freq domain
	    fprintf(out_file, "%f %f\n", creal(buf_t[i]), cimag(buf_t[i]));
	}

	//Closing the output file
    fclose(out_file);
	free(buf_t);
	free(out_t);
 
	return 0;
}
