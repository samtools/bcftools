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
#include <unistd.h>
#include <locale.h>
#include <wchar.h>
#include <unistd.h>     // for isatty
#include "hts_internal.h"
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <regex.h>
#include "../bcftools.h"

#define ASSERT_NOT_NULL(a) do {if(a==NULL) {fprintf(stderr,"[%s:%d]NULL Ptr exception\n",__FILE__,__LINE__);abort();}} while(0)
#define WHERE fprintf(stderr,"[%s:%s:%d]\n",__FUNCTION__,__FILE__,__LINE__)

#define DEFINE_ANSI_IOMANIP(NAME,OPCODE) const char* COLOR_##NAME="\033[" #OPCODE  "m";

DEFINE_ANSI_IOMANIP(RESET,0)
DEFINE_ANSI_IOMANIP(BLACK,30)
DEFINE_ANSI_IOMANIP(RED,31)
DEFINE_ANSI_IOMANIP(GREEN,32)
DEFINE_ANSI_IOMANIP(YELLOW,33)
DEFINE_ANSI_IOMANIP(BLUE,34)
DEFINE_ANSI_IOMANIP(MAGENTA,35)
DEFINE_ANSI_IOMANIP(CYAN,36)
DEFINE_ANSI_IOMANIP(WHITE,37)



KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

KHASH_MAP_INIT_STR(hdict, bcf_hrec_t*)
typedef khash_t(hdict) hdict_t;



typedef unsigned char color_t;

typedef struct Cell {
    kstring_t text;
    kstring_t url;
    const char* color;
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


enum build_t {
	undefined,
	human_hg19,
	human_hg38
	};

typedef struct StringList_t {
    unsigned int size;
    char** strings;
} StringList,*StringListPtr;




typedef struct  {
    bcf_hdr_t* header;
    FILE* out;
    int ascii;
    /** columns for VEP predictions */
    StringListPtr vepTokens;
    /** columns for bcftools csq predictions */
    StringListPtr bcsqTokens;
    /** columns for SNPEFF ANN predictions */
    StringListPtr annTokens;
    regex_t regex_rsid;
    unsigned long n_variants;
    enum build_t build;
    int hide_HOM_REF;
    int hide_NO_CALL;
    }  args_t;

static args_t args;

/** build a new Cell */
static CellPtr CellNew() {
    CellPtr ptr = (CellPtr)calloc(1UL,sizeof(Cell));
    ASSERT_NOT_NULL(ptr);
    ptr->color = COLOR_BLACK;
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
/*
static CellPtr CellSetLL(CellPtr ptr, long long v) {
    CellClear(ptr);
    kputll(v,&(ptr->text));
    return ptr;
    }*/
    
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
 
 /*
 static const char* CellCStr(CellPtr ptr) {
 	ASSERT_NOT_NULL(ptr);
 	return ks_c_str(&(ptr->text));
 	}
 */
 
static void CellPrint(CellPtr ptr,args_t* args) {
    ASSERT_NOT_NULL(ptr);
    if(args->ascii==0 && ptr->color!=NULL && ptr->color!=COLOR_BLACK) {
    	fputs(ptr->color,args->out);
    	}
    fwrite((void*)ks_c_str(&(ptr->text)), sizeof(char), ks_len(&(ptr->text)), args->out);
     if(args->ascii==0 && ptr->color!=NULL && ptr->color!=COLOR_BLACK) {
    	fputs(COLOR_RESET,args->out);
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

static RowPtr RowRemoveAt(RowPtr ptr,unsigned int idx) {
	 ASSERT_NOT_NULL(ptr);
	 assert(idx < ptr->size);
	 CellFree(ptr->cells[idx]);
	 memmove(&ptr->cells[idx], &ptr->cells[idx+1], sizeof(CellPtr)*((ptr->size-1)-idx));
	 ptr->size--;
	 return ptr;
	}

static RowPtr RowAppendStr(RowPtr row,const char* s) {
    return RowAppend(row,CellNewStr(s));
    }


static CellPtr RowAt(RowPtr row,unsigned int idx) {
    assert(idx < RowSize(row));
    return row->cells[idx];
    }

/** set content of idx-th column */
static CellPtr RowSetText(RowPtr row,unsigned int idx,const char* value) {
    return CellSetText(RowAt(row,idx),value);
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
	return RowAt(row,x);
	}
static int TableIsColumnEmpty(TablePtr ptr,unsigned int x) {
	unsigned int y;
	ASSERT_NOT_NULL(ptr);
	assert(x < TableNCols(ptr));
	for(y=0;y < TableNRows(ptr);++y) {
		CellPtr cell = TableAt(ptr,x,y);
		if(CellWidth(cell)!=0) return 0;
		}
	return 1;
	}

static TablePtr TableRemoveColumn(TablePtr ptr,unsigned int x) {
	RowRemoveAt(ptr->header,x);
	for(int i=0;i< ptr->size;i++) {
		RowRemoveAt(ptr->rows[i],x);
		}
	return ptr;
	}

static TablePtr TableRemoveEmptyColumns(TablePtr ptr) {
	unsigned int x=0;
	ASSERT_NOT_NULL(ptr);
	while(x <TableNCols(ptr)) {
		if(TableIsColumnEmpty(ptr,x)) {
			TableRemoveColumn(ptr,x);
			}
		else
			{
			x++;
			}
		}
	return ptr;
	}



static RowPtr TableNewRow(TablePtr ptr) {
    RowPtr row = RowNew(TableNCols(ptr));
    ptr->rows = (RowPtr*)realloc(ptr->rows,(ptr->size+1)*sizeof(RowPtr));
    ASSERT_NOT_NULL(ptr->rows);
    ptr->rows[ptr->size] = row;
    ptr->size++;
    return row;
    }


static StringListPtr StringListNew(const char* str,char delim) {
StringListPtr ptr = calloc(1UL,sizeof(StringList));
ASSERT_NOT_NULL(ptr);
char* prev=(char*)str;
char* p =(char*)str;
for(;;) {
    if(*p==delim || *p==0) {
        ptr->strings=(char**)realloc(ptr->strings,sizeof(char*)*(ptr->size+1));
        ASSERT_NOT_NULL(ptr->strings);
        ptr->strings[ptr->size] = strndup(prev,p-prev);
     	ASSERT_NOT_NULL( ptr->strings[ptr->size]);
        ptr->size++;
        if(*p==0) break;
        prev=p+1;
        }
    p++;
    }
return ptr;
}
/**
Dispose list of String
*/


void StringListFree(StringList* ptr) {
unsigned int i;
if(ptr==NULL) return;
for(i=0;i< ptr->size;++i) {
    free(ptr->strings[i]);
    }
free(ptr);
}

const char* StringListAt(StringList* ptr,unsigned int idx) {
ASSERT_NOT_NULL(ptr);
assert(idx < ptr->size);
return ptr->strings[idx];
}


/**
print symbol
*/


static void printSymbol(args_t* args,unsigned int repeat, const char* wc, char c) {
unsigned int i;
if(args->ascii==1) {
     for(i=0;i< repeat;i++) {
        fputc(c,args->out);
        }
     }
else
    {
    for(i=0;i< repeat;i++) {
        fputs(wc,args->out);
        }
    }
}


static void TablePrint(TablePtr ptr,args_t* args) {
unsigned int y,x;
unsigned int* widths = calloc(TableNCols(ptr),sizeof(unsigned int));
ASSERT_NOT_NULL(ptr);

for(x=0; x<TableNCols(ptr); ++x) {
    unsigned int width =  CellWidth(RowAt(ptr->header,x));
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
			printSymbol(args,1,(x==0?"\u250C":"\u252C"),'+');
			printSymbol(args,2+widths[x],"\u2500",'-');
			}
		printSymbol(args,1,"\u2510",'+');
		fputc('\n',args->out);
		
		//line 2 of header
		for(int x=0;x<TableNCols(ptr);++x) {
			printSymbol(args,1,"\u2502",'|');
			fputc(' ',args->out);
			CellPrint(RowAt(ptr->header,x),args);
			printSymbol(args,widths[x]-CellWidth(RowAt(ptr->header,x))," ",' ');
			fputc(' ',args->out);
			}
		printSymbol(args,1,"\u2502",'|');
		fputc('\n',args->out);
		
		//line 3 of header
		for(int x=0;x<TableNCols(ptr);++x) {
		    if(x==0 && TableNRows(ptr)==0)
				{	
				printSymbol(args,1,"\u2514",'+');
				}
			else if(x==0)
				{
				printSymbol(args,1,"\u251C",'+');
				}
			else
				{
				printSymbol(args,1,(TableNRows(ptr)==0?"\u2534":"\u253C"),'+');
				}
			printSymbol(args,2+widths[x],"\u2500",'-');
			}
		
		
	    printSymbol(args,1,(TableNRows(ptr)==0 ? "\u2518":"\u2524"),'+');
		fputc('\n',args->out);
		
		//print body
		for(y=0;y< TableNRows(ptr);++y) {
			RowPtr row = TableRowAt(ptr,y);
			//line  of data
			for(x=0;x<TableNCols(ptr);++x) {
				CellPtr cell = RowAt(row,x);
				printSymbol(args,1,"\u2502",'|');
				fputc(' ',args->out);
				CellPrint(cell,args);
				printSymbol(args,widths[x]-CellWidth(cell)," ",' ');
				fputc(' ',args->out);
				}
			printSymbol(args,1,"\u2502",'|');
			fputc('\n',args->out);
			}
		//last line
		if(TableNRows(ptr)>0)
			{
			for(x=0;x<TableNCols(ptr);++x) {
				printSymbol(args,1,(x==0?"\u2514":"\u2534"),'+');
				printSymbol(args,2+widths[x],"\u2500",'-');
				}
			printSymbol(args,1,"\u2518",'+');
			fputc('\n',args->out);
			}

free(widths);
}

static void HyperLinkTableAdd(TablePtr table, const char* allele, const char* key, const char* value) {
if(value==NULL || strcmp(value,"")==0) return;
RowPtr row = TableNewRow(table);
CellSetText(RowAt(row,0),value);
}

static int  findContigs(bcf_hdr_t *hdr_in, const char* ctg1a, uint64_t len1, const char* ctg2a, uint64_t len2) {
	char ctg1b[10];
	char ctg2b[10];
	sprintf(ctg1b, "chr%s", ctg1a);
	sprintf(ctg2b, "chr%s", ctg2a);
	int found=0;
    int i,n_contigs= hdr_in->n[BCF_DT_CTG];
    for(i=0;i< n_contigs && found<2 ;i++) {
    	uint64_t len;
    	bcf_idpair_t  c = hdr_in->id[BCF_DT_CTG][i];
    	len = c.val->info[0];
    	const char* contig_name = c.key;
    	if(len == len1 && (strcmp(ctg1a,contig_name)==0 || strcmp(ctg1b,contig_name)==0)) {
    		found++;
    		}
    	else if(len == len2 && (strcmp(ctg2a,contig_name)==0 || strcmp(ctg2b,contig_name)==0)) {
    		found++;
    		}
    	}
    return found==2;
    }

const char *about(void)
{
    return "Convert VCF to table.\n";
}

/*
    Called once at startup, it initializes local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char **argv, bcf_hdr_t *hdr_in, bcf_hdr_t *out)
{
	 int c;
    args.header = hdr_in;
    args.ascii = 0;
    args.out = stdout;
    args.n_variants=0L;
    args.vepTokens = NULL;
    args.bcsqTokens = NULL;
    args.hide_HOM_REF = 0;
    args.hide_NO_CALL = 0;
    args.annTokens = StringListNew(
        "Allele,Annotation,Annotation_Impact,Gene_Name,Gene_ID,"
        "Feature_Type,Feature_ID,Transcript_BioType,Rank,"
        "HGVS.c,HGVS.p,cDNA.pos/length,CDS.pos/length,AA.pos/length,Distance,Message",
        ',');
    
    
	c = regcomp(&args.regex_rsid,"rs[0-9]+",REG_EXTENDED|REG_ICASE|REG_NOSUB);
    assert(c==0);
    
    static struct option loptions[] =
    {
        {"hide",required_argument,NULL,'x'},
        {0,0,0,0}
    };
   
    while ((c = getopt_long(argc, argv, "hx:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'x':
            	{
            	int i;
            	StringListPtr hide = StringListNew(optarg,',');
            	for(i=0; i< hide->size;++i) {
            		const char* hidden =  StringListAt(hide,i);
            		if(strcasecmp(hidden, "HOM_REF")==0) {
            			args.hide_HOM_REF = 1;
            			}
            		else if(strcasecmp(hidden, "NO_CALL")==0 || strcasecmp(hidden, "MISSING")==0) {
            			args.hide_NO_CALL = 1;
            			}
            		}
            	StringListFree(hide);
                break;
                }
            case 'h':
            case '?':
            default: error("wrong arguments"); break;
        }
    }

    if ( !isatty(fileno((FILE *)stdout)) ) {
    	args.ascii=1;
    	}
     
    if(args.ascii==0) {
        if (setlocale(LC_CTYPE, "") == NULL)
	        {
	        fprintf(stderr,"setlocale failed. Switching to ascii\n");
		    args.ascii=1;
	        }
        }
    /* guess the build ?*/
    {
    if( findContigs(hdr_in,"1",249250621,"2",243199373)) {
    	args.build = human_hg19;
    	}
   	else if( findContigs(hdr_in,"1",248956422,"2",242193529)) {
    	args.build = human_hg38;
    	}
   else {
    	args.build = undefined;
    	}
	}
    
    /** find INFO/CSQ and decode it */
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr_in, BCF_HL_INFO, NULL, "CSQ", NULL);
    if(hrec!=NULL) {
    	int ret = bcf_hrec_find_key(hrec, "Description");
		char *format = ret < 0 ? NULL: strstr(hrec->vals[ret], "Format: ");
		if(format!=NULL) {
				format += 8;
				char* vep_format = strdup(format);
				//remove trailing quote
				if(vep_format[strlen(vep_format)-1]=='"') {
					vep_format[strlen(vep_format)-1] = 0;
					}
				args.vepTokens = StringListNew(vep_format,'|');
				free(vep_format);
			}
	   }
	/** find INFO/BCSQ */
    hrec = bcf_hdr_get_hrec(hdr_in, BCF_HL_INFO, NULL, "BCSQ", NULL);
    if(hrec!=NULL) {
    	int ret = bcf_hrec_find_key(hrec, "Description");
		char *format = ret < 0 ? NULL: strstr(hrec->vals[ret], "Format: ");
		if(format!=NULL) {
				format += 8;
				args.bcsqTokens = StringListNew(format,'|');
			}
	   }
            

    return 1;//suppress VCF/BCF header
}



static void escapeHttp(kstring_t* k,const char* s) {
char* p=(char*)s;
while(*p!=0) {
	if(isalnum(*p)) {
		kputc(*p,k);
		}
	else
		{
		kputc(*p,k);
		//ksprintf(k,"%%%02x", (unsigned char)*p);
		}
	p++	;
	}
}

#define PRINT_HEADER \
	switch(args.build) {\
	case human_hg19 : fputs(" GRCh37 : ",args.out); break;\
	case human_hg38 : fputs(" GRCh38 : ",args.out); break;\
	default:break;\
	}\
	fprintf(args.out,"%s:%s:%s (%ld)\n",tokens->strings[0],tokens->strings[1],tokens->strings[3], args.n_variants)


/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t *process(bcf1_t *v) {
    TablePtr vepTable = NULL;
    TablePtr bcsqTable = NULL;
    TablePtr annTable = NULL;
    TablePtr hyperlinksTable = TableNewStr("DB","Allele","URL",NULL);
    unsigned int i;
    args.n_variants++;
    kstring_t vcf_line = KS_INITIALIZE;
    vcf_format(args.header, v, &vcf_line); 
    //remove last CR/LF
    if(vcf_line.s[vcf_line.l-1]=='\n') {
        vcf_line.s[vcf_line.l-1]=0;
        vcf_line.l--;
        }



    StringListPtr tokens = StringListNew(vcf_line.s,'\t');
    StringListPtr alt_alleles = StringListNew(StringListAt(tokens,4),',');


    fputs("<<<",args.out);
    PRINT_HEADER;
     
    TablePtr p = TableNewStr("KEY","VALUE",NULL);

    RowPtr row = TableNewRow(p);
    CellSetText(RowAt(row,0), "CHROM");
    CellSetText(RowAt(row,1),StringListAt(tokens,0));

    row = TableNewRow(p);
    CellSetText(RowAt(row,0), "POS");
    CellSetText(RowAt(row,1),StringListAt(tokens,1));

    row = TableNewRow(p);
    CellSetText(RowAt(row,0), "ID");
    CellSetText(RowAt(row,1),StringListAt(tokens,2));
    if(regexec(&args.regex_rsid, StringListAt(tokens,2),0,NULL,0)==0) {
	    HyperLinkTableAdd(hyperlinksTable,NULL, "RSID", StringListAt(tokens,2));
	    }

    row = TableNewRow(p);
    CellSetText(RowAt(row,0), "REF");
    CellSetText(RowAt(row,1),StringListAt(tokens,3));

    row = TableNewRow(p);
    CellSetText(RowAt(row,0), "ALT");
    CellSetText(RowAt(row,1),StringListAt(tokens,4));

    row = TableNewRow(p);
    CellSetText(RowAt(row,0), "QUAL");
    CellSetText(RowAt(row,1),tokens->strings[5]);

    row = TableNewRow(p);
    CellSetText(RowAt(row,0), "FILTER");
    CellSetText(RowAt(row,1),StringListAt(tokens,6));
    if(strcmp(StringListAt(tokens,6),".")!=0 && strcmp(StringListAt(tokens,6),"PASS")!=0) {
	    RowAt(row,1)->color = COLOR_RED;
	    }

     fprintf(args.out, "# Variant\n");
    TablePrint(p,&args);
    TableFree(p);
    fputc('\n',args.out);


    /* ADD HYPERLINKS */
    if(args.build == human_hg19 || args.build==human_hg38) {
	    kstring_t url = KS_INITIALIZE;
	    for(i=0;i< alt_alleles->size;++i) {
		    ks_clear(&url);
		    const char* alt_allele= StringListAt(alt_alleles,i);
		    
		    
		    RowPtr annot = TableNewRow(hyperlinksTable);
		    CellSetText(RowAt(annot,0), "GNOMAD");
		    CellSetText(RowAt(annot,1), alt_allele);
		    kputs("https://gnomad.broadinstitute.org/variant/",&url);
		    escapeHttp(&url,StringListAt(tokens,0));
		    escapeHttp(&url,"-");
		    kputs(StringListAt(tokens,1),&url);
		    escapeHttp(&url,"-");
		    escapeHttp(&url,StringListAt(tokens,3));
		    escapeHttp(&url,"-");
		    escapeHttp(&url,alt_allele);
		    
		    
		    CellSetText(RowAt(annot,2), url.s);
		    }
	    //			StringUtils.escapeHttp(ensemblCtg) + "-" + ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()+"?dataset=gnomad_r2_1"
	    //HyperLinkTableAdd(hyperlinksTable,NULL, "RSID",url.s);
	    
	    
	    
	    ks_free(&url);
	    }


    if(tokens->size>7 && strcmp(tokens->strings[7],".")!=0) {
        StringListPtr infos = StringListNew(tokens->strings[7],';');
        TablePtr p = TableNewStr("KEY","IDX","VALUE",NULL);
        for(i=0;i< infos->size;i++) {
            unsigned int j;
            const char* info = StringListAt(infos,i);
            char* eq  = strchr(info,'=');
            if(eq==NULL || eq==info) continue;
            
            
           	
            StringListPtr values = StringListNew(eq+1,',');
            for(j=0;j< values->size;j++) {
                //skip CSQ
            	if(args.vepTokens!=NULL && strncmp(info,"CSQ=",4)==0) {
            		unsigned int k;
            		//build VEP table if needed
		        	if(vepTable==NULL) {
		        		vepTable = TableNew(0);
		        		for(k=0;k< args.vepTokens->size;++k) {
		        			TableAppendColumn(vepTable,StringListAt( args.vepTokens,k));
		        			}	
		        		}
		        	// fill VEP table
		        	row = TableNewRow(vepTable);
		        	StringListPtr veps = StringListNew( StringListAt(values,j),'|');
		        	for(k=0;k< args.vepTokens->size && k < veps->size;++k) {
		        			CellSetText(RowAt(row,k),StringListAt( veps,k));
		        			
		        			if(strcmp(StringListAt(args.vepTokens,k), "SYMBOL")==0) {
		        				HyperLinkTableAdd(hyperlinksTable, NULL,StringListAt(args.vepTokens,k), StringListAt( veps,k));
		        				}
		        			
		        			}
		        	StringListFree(veps);
		        	continue;
		        	}
		    
                 //skip BCSQ
            	if(args.bcsqTokens!=NULL && strncmp(info,"BCSQ=",5)==0) {
            		unsigned int k;
            		//build BCSQ table if needed
		        	if(bcsqTable==NULL) {
		        		bcsqTable = TableNew(0);
		        		for(k=0;k< args.bcsqTokens->size;++k) {
		        			TableAppendColumn(bcsqTable,StringListAt( args.bcsqTokens,k));
		        			}	
		        		}
		        
		        	// fill BCSQ table
		        	row = TableNewRow(bcsqTable);
		        	StringListPtr bcsq = StringListNew( StringListAt(values,j),'|');
		        	for(k=0;k< args.bcsqTokens->size && k < bcsq->size;++k) {
		        			CellSetText(RowAt(row,k),StringListAt( bcsq,k));
		        			
		        			if(strcmp(StringListAt(args.bcsqTokens,k), "SYMBOL")==0) {
		        				HyperLinkTableAdd(hyperlinksTable, NULL,StringListAt(args.bcsqTokens,k), StringListAt( bcsq,k));
		        				}
		        			
		        			}
		        	StringListFree(bcsq);
		        	continue;
		        	}
                
                //skip SNPEFF/ANN
                if(args.annTokens!=NULL && strncmp(info,"ANN=",4)==0) {
                    unsigned int k;
                    //build BCSQ table if needed
		        	if(annTable==NULL) {
		        		annTable = TableNew(0);
		        		for(k=0;k< args.annTokens->size;++k) {
		        			TableAppendColumn(annTable,StringListAt( args.annTokens,k));
		        			}	
		        		}
                    // fill ANN table
		        	row = TableNewRow(annTable);
		        	StringListPtr ann = StringListNew( StringListAt(values,j),'|');
		        	for(k=0;k< args.annTokens->size && k < ann->size;++k) {
		        			RowSetText(row,k,StringListAt(ann,k));
		        			}
		            StringListFree(ann);
                    continue;
                    }
                
                row = TableNewRow(p);
                CellAppendTextN(RowAt(row,0),info,eq-info);
                if(values->size>1) CellSetD(RowAt(row,1),(int)(j+1));
                RowSetText(row,2,values->strings[j]);
                }
            StringListFree(values);
            }
        fprintf(args.out, "# INFO\n");
        TablePrint(p,&args);
	    TableFree(p);
        StringListFree(infos);
        fputc('\n',args.out);
        }


    if(TableNRows(hyperlinksTable)>0) {
	    fprintf(args.out, "# HYPERLINKS\n");
	    TablePrint(hyperlinksTable,&args);
	    fputc('\n',args.out);
	    }



    if(vepTable!=NULL && TableNRows(vepTable)>0) {
	    fprintf(args.out, "# VEP/CSQ\n");
	    TableRemoveEmptyColumns(vepTable);
	    TablePrint(vepTable,&args);
	    fputc('\n',args.out);
	    }

    if(bcsqTable!=NULL && TableNRows(bcsqTable)>0) {
	    fprintf(args.out, "# BCSQ\n");
	    TableRemoveEmptyColumns(bcsqTable);
	    TablePrint(bcsqTable,&args);
	    fputc('\n',args.out);
	    }
	
	 if(annTable!=NULL && TableNRows(annTable)>0) {
	    fprintf(args.out, "# ANN/SNPEFF\n");
	    TableRemoveEmptyColumns(annTable);
	    TablePrint(annTable,&args);
	    fputc('\n',args.out);
	    }
	
    if(tokens->size>9) {

        StringListPtr formats = StringListNew(tokens->strings[8],':');
        TablePtr p = TableNewStr("SAMPLE",NULL);
        TableAppendColumn(p, "GTYPE");
        int gt_col = -1;
        for(i=0; i<formats->size;i++) { 
          TableAppendColumn(p, formats->strings[i]);
          if(strcmp("GT",formats->strings[i])==0) gt_col=(int)i;
          }
        
        for(i=9;i< tokens->size;i++) {
            kstring_t gtype_name = KS_INITIALIZE;
            int count_allele_0=0;
            int count_allele_1=0;
            int count_allele_missing=0;
            int count_allele_other=0;
            int print_it = 1;
            StringListPtr values = StringListNew(tokens->strings[i],':');
            const char* color = COLOR_BLACK;
            unsigned int j;
            if(gt_col!=-1) {
                char* gt = strdup(values->strings[gt_col]);
                for(j=0;gt[j]!=0;j++) {
                    if(gt[j]=='|') gt[j]='/';
                    }
                StringListPtr alleles = StringListNew(gt,'/');
                for(j=0;j< alleles->size;++j) {
                    char* allele = alleles->strings[j];
                    if(strcmp(allele,"0")==0) count_allele_0++;
                    else if(strcmp(allele,"1")==0) count_allele_1++;
                    else if(strcmp(allele,".")==0) count_allele_missing++;
                    else count_allele_other++;
                    }
                if(alleles->size==2) {
                    if(count_allele_0==0 && count_allele_1==0 && count_allele_other==0) {kputs("NO_CALL",&gtype_name);if(args.hide_NO_CALL) print_it=0;}
                    else if(count_allele_0==2) { kputs("HOM_REF",&gtype_name); color=COLOR_GREEN;if(args.hide_HOM_REF) print_it=0;}
                    else if(count_allele_missing==0 && strcmp(StringListAt(alleles,0),StringListAt(alleles,1))==0) {kputs("HOM_VAR",&gtype_name); color=COLOR_RED;}
                    else if(count_allele_missing==0 && strcmp(StringListAt(alleles,0),StringListAt(alleles,1))!=0) {kputs("HET",&gtype_name); color=COLOR_CYAN;}
                    }
                else if(alleles->size==1) {
                    if(count_allele_0==1) {kputs("REF",&gtype_name);color=COLOR_GREEN;if(args.hide_HOM_REF) print_it=0;}
                    else if(count_allele_1==1) {kputs("ALT",&gtype_name);color=COLOR_RED;}
                    else if(count_allele_missing==1)  {kputs("NO_CALL",&gtype_name);if(args.hide_NO_CALL) print_it=0;}
                    }
                else
                	{
                	if(count_allele_0==alleles->size) {kputs("HOM_REF",&gtype_name);; color=COLOR_GREEN;if(args.hide_HOM_REF) print_it=0;}
                	else if(count_allele_missing==alleles->size) {kputs("NO_CALL",&gtype_name);if(args.hide_NO_CALL) print_it=0;}
                	}
                StringListFree(alleles);
                free(gt);
                }
           if(print_it) {
               row = TableNewRow(p);
               CellSetText(RowAt(row,0),args.header->samples[i-9]);
               CellSetText(RowAt(row,1),gtype_name.s);
               RowAt(row,1)->color = color;
               
               for(j=0; j< values->size;j++) {
                  CellSetText(RowAt(row,j+2), values->strings[j]);
                  }
               }
           StringListFree(values);
           ks_free(&gtype_name);
           }
        fprintf(args.out, "# GENOTYPES\n");
        TablePrint(p,&args);
        TableFree(p);
        StringListFree(formats);
        }




    fputs(">>>",args.out);
    PRINT_HEADER;

    fputc('\n',args.out);

    ks_free(&vcf_line);
    StringListFree(tokens);
    StringListFree(alt_alleles);
    TableFree(hyperlinksTable);
    TableFree(bcsqTable);
    TableFree(vepTable);
    TableFree(annTable);
    return NULL;/* suppress bcf output */
    }

void destroy(void) {
    regfree(&args.regex_rsid);
    StringListFree(args.vepTokens);
    StringListFree(args.bcsqTokens);
    StringListFree(args.annTokens);
    }


