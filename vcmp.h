#ifndef __VCMP_H__
#define __VCMP_H__

typedef struct _vcmp_t vcmp_t;

vcmp_t *vcmp_init();
void vcmp_destroy(vcmp_t *vcmp);

/*
 *  vcmp_set_ref() - sets and compares reference alleles
 *  Returns 0 on success or -1 if alleles not compatible
 */
int vcmp_set_ref(vcmp_t *vcmp, char *ref1, char *ref2);
int vcmp_find_allele(vcmp_t *vcmp, char **als1, int nals1, char *al2);

#endif
