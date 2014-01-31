#include <stdio.h>
#include <stdlib.h>
#include "rbuf.h"

void debug_print(rbuf_t *rbuf, int *dat)
{
    int i;
    for (i=-1; rbuf_next(rbuf, &i); ) printf(" %d", i); printf("\n");
    for (i=-1; rbuf_next(rbuf, &i); ) printf(" %d", dat[i]); printf("\n");
}

int main(int argc, char **argv)
{
    int i, j, *dat = (int*)calloc(10,sizeof(int));
    rbuf_t rbuf;
    rbuf_init(&rbuf,10);

    rbuf.f = 5; // force wrapping
    for (i=0; i<9; i++)
    {
        j = rbuf_add(&rbuf);
        dat[j] = i+1;
    }
    printf("Inserted 1-9 starting at offset 5:\n");
    debug_print(&rbuf, dat);

    i = rbuf_kth(&rbuf, 3);
    printf("4th is %d\n", dat[i]);

    printf("Deleting 1-2:\n");
    rbuf_shift_n(&rbuf, 2);
    debug_print(&rbuf, dat);

    free(dat);
    return 0;
}

