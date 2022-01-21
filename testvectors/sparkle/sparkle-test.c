#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "sparkle_vectors.h"
#include "sparkle.c"


static bool compare_and_print (size_t len, uint32_t* comp, uint32_t* exp) {
	bool ok = true;
	for (size_t i = 0; i < len; i++)
    	ok = ok & (exp[i] == comp[i]);
  	if (!ok) {
      printf("***FAILURE***\n");
    	printf("computed:");
  		for (size_t i = 0; i < len; i++)
    		printf("0x%0x, " , comp[i]);
      printf("\n");

  		printf("expected:");
  		for (size_t i = 0; i < len; i++)
    		printf("0x%0x, " , exp[i]);
  		printf("\n");
    }
    return ok;
}


static bool sparkle_wrapper(sparkle_tv tv) {
	uint32_t* temp = malloc (sizeof (uint32_t) * tv.input_len);
	for (size_t i = 0; i < tv.input_len; i++)
	{
		temp[i] = tv.input[i];
	}

	sparkle(temp, tv.brans, tv.steps);
	return compare_and_print(tv.input_len, temp, tv.output);
}


int main()
{
	bool ok = true;
  	for (int i = 0; i < sizeof(vectors)/sizeof(sparkle_tv); ++i) {
    	ok &= sparkle_wrapper(vectors[i]);
	}
  if (ok)
    printf("The testing executed successfully! \n");
	return 0;

}