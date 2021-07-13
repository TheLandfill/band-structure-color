#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "xorshift.h"

SHUFFLE_FUNC_DECL(double)

double * gen_band(double low, double high, size_t num);

int main() {
	struct xorshift32_state state = { time(NULL) };
	const int NUM_ELEMENTS = 1e6;
	FILE * data = fopen("results", "w");
	if (!data) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		goto clean_file;
	}
	double * valence_band = gen_band(2.0, 5.0, NUM_ELEMENTS);
	if (!valence_band) {
		goto clean_valence_band;
	}
	double * conduction_band = gen_band(7.2, 10, NUM_ELEMENTS);
	if (!conduction_band) {
		goto clean_conduction_band;
	}
	for (size_t j = 0; j < 50; j++) {
		shuffle_array(valence_band, NUM_ELEMENTS, &state);
		shuffle_array(conduction_band, NUM_ELEMENTS, &state);
		for (size_t i = 0; i < NUM_ELEMENTS; i++) {
			fprintf(data, "%f\n", conduction_band[i] - valence_band[i]);
		}
	}
	free(conduction_band);
	clean_conduction_band:
	free(valence_band);
	clean_valence_band:
	fclose(data);
	clean_file:
	return 0;
}

double * gen_band(double low, double high, size_t num) {
	double * band = malloc(num * sizeof(double));
	if (!band) {
		return band;
	}
	double diff = high - low;
	diff /= (double)num;
	band[0] = low;
	for (size_t i = 1; i < num; i++) {
		band[i] = band[i - 1] + diff;
	}
	return band;
}

SHUFFLE_FUNC_DEFN(double)
