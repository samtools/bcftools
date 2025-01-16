/* The MIT License

   Copyright (c) 2019-2025 Genome Research Ltd.

   Author: Pierre Lindenbaum PhD Institut-du-Thorax. U1087. Nantes. France

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <locale.h>
#include <wchar.h>
#include <unistd.h>     // for isatty
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>

#define ASSERT_NOT_NULL(a) do {if(a==NULL) {fprintf(stderr,"[%s:%d]NULL Ptr exception\n",__FILE__,__LINE__);abort();}} while(0)
#define WHERE fprintf(stderr,"[%s:%s:%d]\n",__FUNCTION__,__FILE__,__LINE__)
typedef unsigned char color_t;
typedef struct RGB {
    color_t r,g,b,a;
} RGB;
typedef struct Cell {
    kstring_t text;
    kstring_t url;
    RGB color;
} Cell,*CellPtr;

typedef struct Row {
    unsigned int size;
    CellPtr* cells;
} Row,*RowPtr;

typedef struct Table {
    RowPtr header;
    unsigned int size;
    RowPtr* rows;
} Table,*TablePtr;


typedef struct KStringArray {
    unsigned int size;
    kstring_t* strings;
} KStringArray,*KStringArrayPtr;


typedef struct
{
bcf_hdr_t* header;
FILE* out;
int ascii;
unsigned long n_variants;
}  args_t;

static args_t args;

/** build a new Cell */
static CellPtr CellNew() {
    CellPtr ptr = (CellPtr)calloc(1UL,sizeof(Cell));
    ASSERT_NOT_NULL(ptr);
    ks_initialize(&(ptr->text));
    ks_initialize(&(ptr->url));
    return ptr;
    }
static CellPtr CellClear(CellPtr ptr) {
    ASSERT_NOT_NULL(ptr);
    ks_clear(&(ptr->text));
    return ptr;
    }
 static CellPtr CellAppendText(CellPtr ptr, const char* s) {
    ASSERT_NOT_NULL(ptr);
    if(s!=NULL) kputs(s,&(ptr->text));
    return ptr;
    }
 
  static CellPtr CellAppendTextN(CellPtr ptr, const char* s,unsigned int n) {
    ASSERT_NOT_NULL(ptr);
    if(s!=NULL) kputsn(s,n,&(ptr->text));
    return ptr;
    }
 
 
/** build a new Cell */
static CellPtr CellSetText(CellPtr ptr, const char* s) {
    CellClear(ptr);
    CellAppendText(ptr,s);
    return ptr;
    }
static CellPtr CellSetLL(CellPtr ptr, long long v) {
    CellClear(ptr);
    kputll(v,&(ptr->text));
    return ptr;
    }
static CellPtr CellSetD(CellPtr ptr, double v) {
    CellClear(ptr);
    kputd(v,&(ptr->text));
    return ptr;
    }


/** build a new Cell with string content */
static CellPtr CellNewStr(const char* s) {
    CellPtr ptr = CellNew();
    ASSERT_NOT_NULL(ptr);
    return CellSetText(ptr,s);
    }
static unsigned int CellWidth(CellPtr ptr) {
    ASSERT_NOT_NULL(ptr);
    return ks_len(&(ptr->text));
    }
 
static void CellPrint(CellPtr ptr,args_t* args) {
    ASSERT_NOT_NULL(ptr);
    if(args->ascii) {
        fwrite((void*)ks_c_str(&(ptr->text)), sizeof(char), ks_len(&(ptr->text)), args->out);
        }
    else
        {
        fwprintf(args->out,L"%s",ks_c_str(&(ptr->text)));
        }
    }
static void CellFree(CellPtr ptr) {
    if(ptr==NULL) return;
    ks_free(&(ptr->text));
    ks_free(&(ptr->url));
    free(ptr);
    }
    
static unsigned int RowSize(RowPtr row) {
    return row->size;
    }

static void RowFree(RowPtr ptr) {
    if(ptr==NULL) return;
    if(ptr->cells!=NULL)  {
        unsigned int i;
        for(i=0;i< ptr->size;++i) {
            CellFree(ptr->cells[i]);
            }
        free(ptr->cells);
        }
    free(ptr);
    }

 
static RowPtr RowNew(unsigned int size) {
    unsigned int i;
    RowPtr ptr = (RowPtr)calloc(1UL,sizeof(Row));
    ASSERT_NOT_NULL(ptr);
    ptr->cells = (CellPtr*)calloc(size,sizeof(CellPtr));
    ASSERT_NOT_NULL(ptr->cells);
    ptr->size = size;
    for(i=0;i< size;++i) {
        ptr->cells[i] = CellNew();
        ASSERT_NOT_NULL(ptr->cells[i]);
        }
    return ptr;
    }

static RowPtr RowAppend(RowPtr ptr,CellPtr cell) {
    ASSERT_NOT_NULL(ptr);
    ASSERT_NOT_NULL(cell);
    ptr->cells = (CellPtr*)realloc(ptr->cells,(ptr->size+1)*sizeof(CellPtr));
    ASSERT_NOT_NULL(ptr->cells);
    ptr->cells[ptr->size] = cell;
    ptr->size++;
    return ptr;
    }

static RowPtr RowAppendStr(RowPtr row,const char* s) {
    return RowAppend(row,CellNewStr(s));
    }


static CellPtr RowGet(RowPtr row,unsigned int idx) {
    assert(idx < RowSize(row));
    return row->cells[idx];
    }

static unsigned int TableNCols(TablePtr t) {
return RowSize(t->header);
}

static unsigned int TableNRows(TablePtr t) {
ASSERT_NOT_NULL(t);
return t->size;
}

static void TableFree(TablePtr ptr) {
    if(ptr==NULL) return;
    RowFree(ptr->header);
    if(ptr->rows!=NULL)  {
        unsigned int i;
        for(i=0;i< ptr->size;++i) {
            RowFree(ptr->rows[i]);
            }
        free(ptr->rows);
        }
    free(ptr);
    }    

static TablePtr TableNew(unsigned int ncols) {
    TablePtr ptr = (TablePtr)(calloc(1UL,sizeof(Table)));
    ASSERT_NOT_NULL(ptr);
    ptr->size = 0UL;
    ptr->rows = NULL;
    ptr->header= RowNew(ncols);
    ASSERT_NOT_NULL(ptr->header);
    return ptr;
    }
 
static RowPtr TableRowAt(TablePtr ptr,unsigned int y) {
ASSERT_NOT_NULL(ptr);
assert(y < TableNRows(ptr));
ASSERT_NOT_NULL(ptr->rows);
ASSERT_NOT_NULL(ptr->rows[y]);
return ptr->rows[y];
}

static TablePtr TableAppendColumn(TablePtr ptr,const char* title) {
    unsigned int y;
    ASSERT_NOT_NULL(ptr);
    RowAppendStr(ptr->header,title);
    for(y=0;y< TableNRows(ptr);++y) {
        RowAppendStr(TableRowAt(ptr,y),"");
        }
    return ptr;
    }
   
static TablePtr TableNewStr(const char* str,...) {
    va_list arg;
    TablePtr ptr = TableNew(0UL);
    ASSERT_NOT_NULL(ptr);
    va_start(arg, str);
    while (str) {
        TableAppendColumn(ptr,str);
        str = va_arg(arg, const char *);
        }
    va_end(arg);
    return ptr;
    }





static CellPtr TableAt(TablePtr ptr,unsigned int x,unsigned int y) {
RowPtr row = TableRowAt(ptr,y);
return RowGet(row,x);
}

static RowPtr TableNewRow(TablePtr ptr) {
    RowPtr row = RowNew(TableNCols(ptr));
    ptr->rows = (RowPtr*)realloc(ptr->rows,(ptr->size+1)*sizeof(RowPtr));
    ASSERT_NOT_NULL(ptr->rows);
    ptr->rows[ptr->size] = row;
    ptr->size++;
    return row;
    }

static void printSymbol(args_t* args,unsigned int repeat, wchar_t wc, char c) {
unsigned int i;
if(args->ascii==1) {
     for(i=0;i< repeat;i++) {
        fputc(c,args->out);
        }
     }
else
    {
    for(i=0;i< repeat;i++) {
        fputwc(wc,args->out);
        }
    }
}
#define FPUTC(C) do { if(args->ascii) fputc(C,args->out); else fputwc(C,args->out);} while(0)

static void TablePrint(TablePtr ptr,args_t* args) {
unsigned int y,x;
unsigned int* widths = calloc(TableNCols(ptr),sizeof(unsigned int));
ASSERT_NOT_NULL(ptr);

for(x=0; x<TableNCols(ptr); ++x) {
    unsigned int width =  CellWidth(RowGet(ptr->header,x));
    if(width>widths[x]) widths[x] = width;
    }

for(y=0;y< TableNRows(ptr);++y) {
    for(x=0; x<TableNCols(ptr); ++x) {
        ASSERT_NOT_NULL(TableAt(ptr,x,y));
        unsigned int width =  CellWidth(TableAt(ptr,x,y));
        if(width>widths[x]) widths[x] = width;
        }
    }
	
		//print header
		
		// line 1 of header
		for(x=0;x<TableNCols(ptr);++x) {
			printSymbol(args,1,(x==0?L'\u250C':L'\u252C'),'+');
			printSymbol(args,2+widths[x],L'\u2500','-');
			}
		printSymbol(args,1,L'\u2510','+');
		FPUTC('\n');
		
		//line 2 of header
		for(int x=0;x<TableNCols(ptr);++x) {
			printSymbol(args,1,L'\u2502','|');
			FPUTC(' ');
			CellPrint(RowGet(ptr->header,x),args);
			printSymbol(args,widths[x]-CellWidth(RowGet(ptr->header,x)),' ',' ');
			FPUTC(' ');
			}
		printSymbol(args,1,L'\u2502','|');
		FPUTC('\n');
		
		//line 3 of header
		for(int x=0;x<TableNCols(ptr);++x) {
		    if(x==0 && TableNRows(ptr)==0)
				{	
				printSymbol(args,1,L'\u2514','+');
				}
			else if(x==0)
				{
				printSymbol(args,1,L'\u251C','+');
				}
			else
				{
				printSymbol(args,1,(TableNRows(ptr)==0?L'\u2534':L'\u253C'),'+');
				}
			printSymbol(args,2+widths[x],L'\u2500','-');
			}
		
		
	    printSymbol(args,1,(TableNRows(ptr)==0 ? L'\u2518':L'\u2524'),'+');
		FPUTC('\n');
		
		//print body
		for(y=0;y< TableNRows(ptr);++y) {
			RowPtr row = TableRowAt(ptr,y);
			//line  of data
			for(x=0;x<TableNCols(ptr);++x) {
				CellPtr cell = RowGet(row,x);
				printSymbol(args,1,L'\u2502','|');
				FPUTC(' ');
				CellPrint(cell,args);
				printSymbol(args,widths[x]-CellWidth(cell),' ',' ');
				FPUTC(' ');
				}
			printSymbol(args,1,L'\u2502','|');
		FPUTC('\n');
			}
		//last line
		if(TableNRows(ptr)>0)
			{
			for(x=0;x<TableNCols(ptr);++x) {
				printSymbol(args,1,(x==0?L'\u2514':L'\u2534'),'+');
				printSymbol(args,2+widths[x],L'\u2500','-');
				}
			printSymbol(args,1,L'\u2518','+');
		FPUTC('\n');
			}

free(widths);
}

const char *about(void)
{
    return "Convert VCF to table.\n";
}

/*
    Called once at startup, it initializes local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args.header = in;
    args.ascii = 0;
    args.out = stdout;
    args.n_variants=0L;
    
    if(args.ascii==0) {
        if (setlocale(LC_CTYPE, "") == NULL)
	        {
	        fprintf(stderr,"setlocale failed. Switching to ascii\n");
		    args.ascii=1;
	        }
        }
    
    return 1;//suppress VCF/BCF header
}


static KStringArrayPtr KStringArrayNew(const char* str,char delim) {
KStringArrayPtr ptr = calloc(1UL,sizeof(KStringArray));
ASSERT_NOT_NULL(ptr);
char* prev=(char*)str;
char* p =(char*)str;
for(;;) {
    if(*p==delim || *p==0) {
        ptr->strings=(kstring_t*)realloc(ptr->strings,sizeof(kstring_t)*(ptr->size+1));
        ASSERT_NOT_NULL(ptr->strings);
        ks_initialize(&ptr->strings[ptr->size]);
        kputsn(prev,p-prev,&ptr->strings[ptr->size]);
        ptr->size++;
        if(*p==0) break;
        prev=p+1;
        }
    p++;
    }
return ptr;
}
void KStringArrayFree(KStringArray* ptr) {
unsigned int i;
if(ptr==NULL) return;
for(i=0;i< ptr->size;++i) {
    ks_free(&ptr->strings[i]);
    }
free(ptr);
}


/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t *process(bcf1_t *v) {
unsigned int i;
args.n_variants++;
bcf_hdr_t *h = args.header;
kstring_t vcf_line = KS_INITIALIZE;
vcf_format(args.header, v, &vcf_line); 
//remove last CR/LF
if(vcf_line.s[vcf_line.l-1]=='\n') {
    vcf_line.s[vcf_line.l-1]=0;
    vcf_line.l--;
    }
KStringArrayPtr tokens = KStringArrayNew(vcf_line.s,'\t');
 
TablePtr p = TableNewStr("KEY","VALUE",NULL);

RowPtr row = TableNewRow(p);
CellSetText(RowGet(row,0), "CHROM");
CellSetText(RowGet(row,1), bcf_seqname(h, v));

row = TableNewRow(p);
CellSetText(RowGet(row,0), "POS");
CellSetLL(RowGet(row,1), v->pos + 1);

row = TableNewRow(p);
CellSetText(RowGet(row,0), "ID");
CellSetText(RowGet(row,1),tokens->strings[2].s);

row = TableNewRow(p);
CellSetText(RowGet(row,0), "REF");
CellSetText(RowGet(row,1),tokens->strings[3].s);

row = TableNewRow(p);
CellSetText(RowGet(row,0), "ALT");
CellSetText(RowGet(row,1),tokens->strings[4].s);

row = TableNewRow(p);
CellSetText(RowGet(row,0), "QUAL");
CellSetText(RowGet(row,1),tokens->strings[5].s);

row = TableNewRow(p);
CellSetText(RowGet(row,0), "FILTER");
CellSetText(RowGet(row,1),tokens->strings[6].s);


TablePrint(p,&args);
TableFree(p);

if(tokens->size>7 && strcmp(tokens->strings[7].s,".")!=0) {
    KStringArrayPtr infos = KStringArrayNew(tokens->strings[7].s,';');
    TablePtr p = TableNewStr("KEY","IDX","VALUE",NULL);
    for(i=0;i< infos->size;i++) {
        unsigned int j;
        char* info = infos->strings[i].s;
        char* eq  = strchr(info,'=');
        if(eq==NULL || eq==info) continue;
        KStringArrayPtr values = KStringArrayNew(eq+1,',');
        for(j=0;j< values->size;j++) {
            row = TableNewRow(p);
            CellAppendTextN(RowGet(row,0),info,eq-info);
            if(values->size>1) CellSetD(RowGet(row,1),(int)(j+1));
            CellSetText(RowGet(row,2),values->strings[j].s);
            }
        KStringArrayFree(values);
        }
    TablePrint(p,&args);
    TableFree(p);
    KStringArrayFree(infos);
    }


//fputc('\n',args.out);

if(tokens->size>9) {
    
    KStringArrayPtr formats = KStringArrayNew(tokens->strings[8].s,':');
    TablePtr p = TableNewStr("SAMPLE",NULL);
    TableAppendColumn(p, "GTYPE");
    int gt_col = -1;
    for(i=0; i<formats->size;i++) { 
      TableAppendColumn(p, formats->strings[i].s);
      if(strcmp("GT",formats->strings[i].s)==0) gt_col=(int)i;
      }
    
    for(i=9;i< tokens->size;i++) {
        kstring_t gtype_name = KS_INITIALIZE;
        int count_allele_0=0;
        int count_allele_1=0;
        int count_allele_missing=0;
        int count_allele_other=0;
        int print_it = 1;
        KStringArrayPtr values = KStringArrayNew(tokens->strings[i].s,':');
        unsigned int j;
        if(gt_col!=-1) {
            char* gt = strdup(values->strings[gt_col].s);
            for(j=0;gt[j]!=0;j++) {
                if(gt[j]=='|') gt[j]='/';
                }
            KStringArrayPtr alleles = KStringArrayNew(gt,'/');
            for(j=0;j< alleles->size;++j) {
                char* allele = alleles->strings[j].s;
                if(strcmp(allele,"0")==0) count_allele_0++;
                else if(strcmp(allele,"1")==0) count_allele_1++;
                else if(strcmp(allele,".")==0) count_allele_missing++;
                else count_allele_other++;
                }
            if(alleles->size==2) {
                if(count_allele_0==0 && count_allele_1==0 && count_allele_other==0) kputs("NO_CALL",&gtype_name);
                else if(count_allele_0==2) kputs("HOM_REF",&gtype_name);
                else if(count_allele_1==2) kputs("HOM_VAR",&gtype_name);
                else if(count_allele_0==1 && count_allele_1==1) kputs("HET",&gtype_name);
                }
            else if(alleles->size==1) {
                if(count_allele_0==1) kputs("REF",&gtype_name);
                else if(count_allele_1==1) kputs("VAR",&gtype_name);
                else if(count_allele_missing==1) kputs("NO_CALL",&gtype_name);
                }
            KStringArrayFree(alleles);
            free(gt);
            }
       if(print_it) {
           row = TableNewRow(p);
           CellSetText(RowGet(row,0),args.header->samples[i-9]);
           CellSetText(RowGet(row,1),gtype_name.s);
           
           for(j=0; j< values->size;j++) {
              CellSetText(RowGet(row,j+2), values->strings[j].s);
              }
           }
       KStringArrayFree(values);
       ks_free(&gtype_name);
       }
    TablePrint(p,&args);
    TableFree(p);
    KStringArrayFree(formats);
    }


fflush(args.out);
fputc('\n',args.out);





for(i=0;i< tokens->size;i++) {
    fprintf(stderr,"[%d] = %s\n",i,tokens->strings[i].s);
    }
fflush(args.out);
fprintf(args.out,"<<< %s:%s:%s\n",tokens->strings[0].s,tokens->strings[1].s,tokens->strings[3].s);
fflush(args.out);

ks_free(&vcf_line);
KStringArrayFree(tokens);
return NULL;
}

void destroy(void)
{
}


