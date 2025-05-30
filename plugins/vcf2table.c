/* The MIT License

   Copyright (c) 2019-2025  Pierre Lindenbaum

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
#include <getopt.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>  // for isatty

#include "../bcftools.h"

#define ASSERT_NOT_NULL(a)                                                \
  do {                                                                    \
    if (a == NULL) {                                                      \
      fprintf(stderr, "[%s:%d]NULL Ptr exception\n", __FILE__, __LINE__); \
      abort();                                                            \
    }                                                                     \
  } while (0)
#define WHERE fprintf(stderr, "[%s:%s:%d]\n", __FUNCTION__, __FILE__, __LINE__)

#define DEFINE_ANSI_IOMANIP(NAME, OPCODE) \
  const char* COLOR_##NAME = "\033[" #OPCODE "m";

/** colors for ANSI output */
DEFINE_ANSI_IOMANIP(RESET, 0)
DEFINE_ANSI_IOMANIP(BLACK, 30)
DEFINE_ANSI_IOMANIP(RED, 31)
DEFINE_ANSI_IOMANIP(GREEN, 32)
DEFINE_ANSI_IOMANIP(YELLOW, 33)
DEFINE_ANSI_IOMANIP(BLUE, 34)
DEFINE_ANSI_IOMANIP(MAGENTA, 35)
DEFINE_ANSI_IOMANIP(CYAN, 36)
DEFINE_ANSI_IOMANIP(WHITE, 37)

/** A Cell in a Row in a table */
typedef struct Cell {
  kstring_t text;  // text content
  // kstring_t url; for future use
  const char* color;
} Cell, *CellPtr;

/** a Row in a table */
typedef struct Row {
  // number of cells
  unsigned int size;
  // malloc-ed cells
  CellPtr* cells;
} Row, *RowPtr;

/** a table of data */
typedef struct Table {
  // special row containing the table header = number of column
  RowPtr header;
  // number of rows
  unsigned int size;
  // malloc-ed rows
  RowPtr* rows;
} Table, *TablePtr;

/** a list of C-strings */
typedef struct StringList_t {
  // number of strings
  unsigned int size;
  // malloc-ed strings
  char** strings;
} StringList, *StringListPtr;

/** 'common' build */
enum build_t {
  undefined,
  human_hg19,
  human_hg38,
  rotavirus_rf  // for tests...
};

/** global arguments */
typedef struct {
  /** vcf header in */
  bcf_hdr_t* header;
  /** output stream (stdout) */
  FILE* out;
  /** force ascii only */
  int ascii;
  /** columns for VEP predictions */
  StringListPtr vepTokens;
  /** columns for bcftools csq predictions */
  StringListPtr bcsqTokens;
  /** table for SNPEFF ANN predictions */
  TablePtr annTable;
  /** table for SNPEFF LOF predictions */
  TablePtr lofTable;
  /** general info about the variant */
  TablePtr vcTable;
  /** table for hyperlinks */
  TablePtr hyperlinksTable;
  /** table for spliceai */
  TablePtr spliceaiTable;
  /** table for INFO col */
  TablePtr infoTable;
  /** table for genotype types */
  TablePtr gtypeTable;
  /** number of variant seen so far */
  unsigned long n_variants;
  /** genome build */
  enum build_t build;

  /* show/hide flags */
  int hide_HOM_REF;
  int hide_NO_CALL;
  int hide_HET;
  int hide_HOM_VAR;
  int hide_OTHER;
  int hide_VC_table;
  int hide_INFO_table;
  int hide_GT_table;
  int hide_GTTYPE_table;
  int hide_VEP_table;
  int hide_BCSQ_table;
  int hide_ANN_table;
  int hide_LOF_table;
  int hide_SPLICEAI_table;
  int hide_colors;
  int hide_links;
} args_t;

/** global arguments for this plugin */
static args_t args;

/** build a new Cell */
static CellPtr CellNew() {
  CellPtr ptr = (CellPtr)calloc(1UL, sizeof(Cell));
  ASSERT_NOT_NULL(ptr);
  ptr->color = COLOR_BLACK;
  ks_initialize(&(ptr->text));
  // ks_initialize(&(ptr->url)); future use
  return ptr;
}

/** clear the content of a cell */
static CellPtr CellClear(CellPtr ptr) {
  ASSERT_NOT_NULL(ptr);
  ks_clear(&(ptr->text));
  return ptr;
}

/** append the content of a cell */
static CellPtr CellAppendText(CellPtr ptr, const char* s) {
  ASSERT_NOT_NULL(ptr);
  if (s != NULL) kputs(s, &(ptr->text));
  return ptr;
}

/** append n bytes to the content of a cell */
static CellPtr CellAppendTextN(CellPtr ptr, const char* s, unsigned int n) {
  ASSERT_NOT_NULL(ptr);
  if (s != NULL) kputsn(s, n, &(ptr->text));
  return ptr;
}

/** set content from a C string */
static CellPtr CellSetText(CellPtr ptr, const char* s) {
  CellClear(ptr);
  CellAppendText(ptr, s);
  return ptr;
}

/** set content from an integer */
static CellPtr CellSetLL(CellPtr ptr, long long v) {
  CellClear(ptr);
  kputll(v, &(ptr->text));
  return ptr;
}

/** set content from an floating number */
static CellPtr CellSetD(CellPtr ptr, double v) {
  CellClear(ptr);
  kputd(v, &(ptr->text));
  return ptr;
}

/** build a new Cell with string content */
static CellPtr CellNewStr(const char* s) {
  CellPtr ptr = CellNew();
  ASSERT_NOT_NULL(ptr);
  return CellSetText(ptr, s);
}

/** return the length of the content for this cell */
static unsigned int CellWidth(CellPtr ptr) {
  ASSERT_NOT_NULL(ptr);
  return ks_len(&(ptr->text));
}

/** print the content of a cell */
static void CellPrint(CellPtr ptr, args_t* args) {
  int color_flag = 0;
  ASSERT_NOT_NULL(ptr);
  // begin color
  if (args->ascii == 0 && ptr->color != NULL && ptr->color != COLOR_BLACK) {
    color_flag = 1;
    fputs(ptr->color, args->out);
  }
  fwrite((void*)ks_c_str(&(ptr->text)), sizeof(char), ks_len(&(ptr->text)),
         args->out);

  // end color
  if (color_flag) {
    fputs(COLOR_RESET, args->out);
  }
}

/** return the content of a Cell as a C-string */
static const char* CellCStr(CellPtr ptr) {
  ASSERT_NOT_NULL(ptr);
  return ks_c_str(&(ptr->text));
}

/** destroy a cell */
static void CellFree(CellPtr ptr) {
  if (ptr == NULL) return;
  ks_free(&(ptr->text));
  // ks_free(&(ptr->url));
  free(ptr);
}

/** return the number of cell in a row */
static unsigned int RowSize(RowPtr row) { return row->size; }

/** destroy a RowPtr */
static void RowFree(RowPtr ptr) {
  if (ptr == NULL) return;
  if (ptr->cells != NULL) {
    unsigned int i;
    for (i = 0; i < ptr->size; ++i) {
      CellFree(ptr->cells[i]);
    }
    free(ptr->cells);
  }
  free(ptr);
}

/** create a RowPtr with 'size' empty cells */
static RowPtr RowNew(unsigned int size) {
  unsigned int i;
  RowPtr ptr = (RowPtr)calloc(1UL, sizeof(Row));
  ASSERT_NOT_NULL(ptr);
  ptr->cells = (CellPtr*)calloc(size, sizeof(CellPtr));
  ASSERT_NOT_NULL(ptr->cells);
  ptr->size = size;
  for (i = 0; i < size; ++i) {
    ptr->cells[i] = CellNew();
    ASSERT_NOT_NULL(ptr->cells[i]);
  }
  return ptr;
}

/** append the cell to this Row */
static RowPtr RowAppend(RowPtr ptr, CellPtr cell) {
  ASSERT_NOT_NULL(ptr);
  ASSERT_NOT_NULL(cell);
  ptr->cells = (CellPtr*)realloc(ptr->cells, (ptr->size + 1) * sizeof(CellPtr));
  ASSERT_NOT_NULL(ptr->cells);
  ptr->cells[ptr->size] = cell;
  ptr->size++;
  return ptr;
}

/** remove the idx-th cell of a row */
static RowPtr RowRemoveAt(RowPtr ptr, unsigned int idx) {
  ASSERT_NOT_NULL(ptr);
  assert(idx < ptr->size);
  CellFree(ptr->cells[idx]);
  memmove((void*)&ptr->cells[idx], (void*)&ptr->cells[idx + 1],
          sizeof(CellPtr) * ((ptr->size - 1) - idx));
  ptr->size--;
  return ptr;
}

/** append a new Cell with the content 's' */
static RowPtr RowAppendStr(RowPtr row, const char* s) {
  return RowAppend(row, CellNewStr(s));
}

/** return the idx-th row */
static CellPtr RowAt(RowPtr row, unsigned int idx) {
  assert(idx < RowSize(row));
  return row->cells[idx];
}

/** set content of idx-th column */
static CellPtr RowSetText(RowPtr row, unsigned int idx, const char* value) {
  return CellSetText(RowAt(row, idx), value);
}

/** return the number of columns in a table */
static unsigned int TableNCols(TablePtr t) { return RowSize(t->header); }

/** return the number of rows in a table */
static unsigned int TableNRows(TablePtr t) {
  ASSERT_NOT_NULL(t);
  return t->size;
}

/** remove all lines of a TablePtr */
static TablePtr TableClear(TablePtr ptr) {
  unsigned int i;
  ASSERT_NOT_NULL(ptr);
  for (i = 0; i < ptr->size; ++i) {
    RowFree(ptr->rows[i]);
    ptr->rows[i] = NULL;
  }
  ptr->size = 0UL;
  return ptr;
}

/** dispose a TablePtr */
static void TableFree(TablePtr ptr) {
  if (ptr == NULL) return;
  TableClear(ptr);
  free(ptr->rows);
  RowFree(ptr->header);
  free(ptr);
}

/** create a new table with 'ncrols' columns */
static TablePtr TableNew(unsigned int ncols) {
  TablePtr ptr = (TablePtr)(calloc(1UL, sizeof(Table)));
  ASSERT_NOT_NULL(ptr);
  ptr->size = 0UL;
  ptr->rows = NULL;
  ptr->header = RowNew(ncols);
  ASSERT_NOT_NULL(ptr->header);
  return ptr;
}

/** return the y-th row in the table */
static RowPtr TableRowAt(TablePtr ptr, unsigned int y) {
  ASSERT_NOT_NULL(ptr);
  assert(y < TableNRows(ptr));
  ASSERT_NOT_NULL(ptr->rows);
  ASSERT_NOT_NULL(ptr->rows[y]);
  return ptr->rows[y];
}

/** append a new column named 'title' in the table */
static TablePtr TableAppendColumn(TablePtr ptr, const char* title) {
  unsigned int y;
  ASSERT_NOT_NULL(ptr);
  RowAppendStr(ptr->header, title);
  for (y = 0; y < TableNRows(ptr); ++y) {
    RowAppendStr(TableRowAt(ptr, y), "");
  }
  return ptr;
}

/** create a new empty table with the header 'str' until NULL */
static TablePtr TableNewStr(const char* str, ...) {
  va_list arg;
  TablePtr ptr = TableNew(0UL);
  ASSERT_NOT_NULL(ptr);
  va_start(arg, str);
  while (str) {
    TableAppendColumn(ptr, str);
    str = va_arg(arg, const char*);
  }
  va_end(arg);
  return ptr;
}

/** return the x-th Cell in the y-th row */
static CellPtr TableAt(TablePtr ptr, unsigned int x, unsigned int y) {
  return RowAt(TableRowAt(ptr, y), x);
}

/** test wether the content of a column is empty */
static int TableIsColumnEmpty(TablePtr ptr, unsigned int x) {
  unsigned int y;
  ASSERT_NOT_NULL(ptr);
  assert(x < TableNCols(ptr));
  for (y = 0; y < TableNRows(ptr); ++y) {
    CellPtr cell = TableAt(ptr, x, y);
    if (CellWidth(cell) != 0) return 0;
  }
  return 1;
}

/** remove the x-th column in the table */
static TablePtr TableRemoveColumn(TablePtr ptr, unsigned int x) {
  RowRemoveAt(ptr->header, x);
  for (int i = 0; i < ptr->size; i++) {
    RowRemoveAt(ptr->rows[i], x);
  }
  return ptr;
}

/** remove any empty column in the table */
static TablePtr TableRemoveEmptyColumns(TablePtr ptr) {
  unsigned int x = 0;
  ASSERT_NOT_NULL(ptr);
  while (x < TableNCols(ptr)) {
    if (TableIsColumnEmpty(ptr, x)) {
      TableRemoveColumn(ptr, x);
    } else {
      x++;
    }
  }
  return ptr;
}

/** create and insert a new row in the table, return this new row */
static RowPtr TableNewRow(TablePtr ptr) {
  RowPtr row = RowNew(TableNCols(ptr));
  ptr->rows = (RowPtr*)realloc(ptr->rows, (ptr->size + 1) * sizeof(RowPtr));
  ASSERT_NOT_NULL(ptr->rows);
  ptr->rows[ptr->size] = row;
  ptr->size++;
  return row;
}

/**
 * Create a new empty StringList
 */
static StringListPtr StringListNew() {
  StringListPtr ptr = (StringListPtr)calloc(1UL, sizeof(StringList));
  ASSERT_NOT_NULL(ptr);
  return ptr;
}

/**
 * Create a new  StringList by splitting 'str' with 'delim'
 */
static StringListPtr StringListMake(const char* str, char delim) {
  StringListPtr ptr = StringListNew();
  char* prev = (char*)str;
  char* p = (char*)str;
  for (;;) {
    if (*p == delim || *p == 0) {
      size_t len = p - prev;
      ptr->strings = (char**)realloc(ptr->strings, sizeof(char*) * (ptr->size + 1));
      ASSERT_NOT_NULL(ptr->strings);
      /* strndup not implemented on windows , use malloc */
      ptr->strings[ptr->size] = malloc((1 + len) * sizeof(char));
      ASSERT_NOT_NULL(ptr->strings[ptr->size]);
      memcpy(ptr->strings[ptr->size], prev, len);
      ptr->strings[ptr->size][len] = 0;
      ptr->size++;
      if (*p == 0) break;
      prev = p + 1;
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
  if (ptr == NULL) return;
  for (i = 0; i < ptr->size; ++i) {
    free(ptr->strings[i]);
  }
  free(ptr->strings);
  free(ptr);
}

/** return content of idx-th item as a const char* */
const char* StringListAt(StringList* ptr, unsigned int idx) {
  ASSERT_NOT_NULL(ptr);
  assert(idx < ptr->size);
  return ptr->strings[idx];
}

/**
print symbol used by TablePrint to print multiple unicode/plain characters
*/
static void printSymbol(args_t* args, unsigned int repeat, const char* wc,
                        char c) {
  unsigned int i;
  if (args->ascii == 1) {
    for (i = 0; i < repeat; i++) {
      fputc(c, args->out);
    }
  } else {
    for (i = 0; i < repeat; i++) {
      fputs(wc, args->out);
    }
  }
}

/** print the content of a table */
static void TablePrint(TablePtr ptr, args_t* args) {
  unsigned int y, x;
  unsigned int* widths = calloc(TableNCols(ptr), sizeof(unsigned int));
  ASSERT_NOT_NULL(ptr);

  for (x = 0; x < TableNCols(ptr); ++x) {
    unsigned int width = CellWidth(RowAt(ptr->header, x));
    if (width > widths[x]) widths[x] = width;
  }

  for (y = 0; y < TableNRows(ptr); ++y) {
    for (x = 0; x < TableNCols(ptr); ++x) {
      ASSERT_NOT_NULL(TableAt(ptr, x, y));
      unsigned int width = CellWidth(TableAt(ptr, x, y));
      if (width > widths[x]) widths[x] = width;
    }
  }

  // print header

  // line 1 of header
  for (x = 0; x < TableNCols(ptr); ++x) {
    printSymbol(args, 1, (x == 0 ? "\u250C" : "\u252C"), '+');
    printSymbol(args, 2 + widths[x], "\u2500", '-');
  }
  printSymbol(args, 1, "\u2510", '+');
  fputc('\n', args->out);

  // line 2 of header
  for (int x = 0; x < TableNCols(ptr); ++x) {
    printSymbol(args, 1, "\u2502", '|');
    fputc(' ', args->out);
    CellPrint(RowAt(ptr->header, x), args);
    printSymbol(args, widths[x] - CellWidth(RowAt(ptr->header, x)), " ", ' ');
    fputc(' ', args->out);
  }
  printSymbol(args, 1, "\u2502", '|');
  fputc('\n', args->out);

  // line 3 of header
  for (int x = 0; x < TableNCols(ptr); ++x) {
    if (x == 0 && TableNRows(ptr) == 0) {
      printSymbol(args, 1, "\u2514", '+');
    } else if (x == 0) {
      printSymbol(args, 1, "\u251C", '+');
    } else {
      printSymbol(args, 1, (TableNRows(ptr) == 0 ? "\u2534" : "\u253C"), '+');
    }
    printSymbol(args, 2 + widths[x], "\u2500", '-');
  }

  printSymbol(args, 1, (TableNRows(ptr) == 0 ? "\u2518" : "\u2524"), '+');
  fputc('\n', args->out);

  // print body
  for (y = 0; y < TableNRows(ptr); ++y) {
    RowPtr row = TableRowAt(ptr, y);
    // line  of data
    for (x = 0; x < TableNCols(ptr); ++x) {
      CellPtr cell = RowAt(row, x);
      printSymbol(args, 1, "\u2502", '|');
      fputc(' ', args->out);
      CellPrint(cell, args);
      printSymbol(args, widths[x] - CellWidth(cell), " ", ' ');
      fputc(' ', args->out);
    }
    printSymbol(args, 1, "\u2502", '|');
    fputc('\n', args->out);
  }
  // last line
  if (TableNRows(ptr) > 0) {
    for (x = 0; x < TableNCols(ptr); ++x) {
      printSymbol(args, 1, (x == 0 ? "\u2514" : "\u2534"), '+');
      printSymbol(args, 2 + widths[x], "\u2500", '-');
    }
    printSymbol(args, 1, "\u2518", '+');
    fputc('\n', args->out);
  }
  fputc('\n', args->out);
  free(widths);
}

/**
 * This method is used to find if a dictionary contains two known contig
 * (name/length) in order to identify the build : hg19, hg38, etc...
 */
static int findContigs(bcf_hdr_t* hdr_in, const char* ctg1a, uint64_t len1,
                       const char* ctg2a, uint64_t len2) {
  char ctg1b[10];
  char ctg2b[10];
  // try to add a 'chr' prefix to the chromosome name
  sprintf(ctg1b, "chr%s", ctg1a);
  sprintf(ctg2b, "chr%s", ctg2a);
  int found = 0;
  int i, n_contigs = hdr_in->n[BCF_DT_CTG];
  for (i = 0; i < n_contigs && found < 2; i++) {
    uint64_t len;
    bcf_idpair_t c = hdr_in->id[BCF_DT_CTG][i];
    if (c.val == NULL) continue;
    //if (c.val->info == NULL) continue; no info is always an array
    len = c.val->info[0];
    const char* contig_name = c.key;
    if (len == len1 &&
        (strcmp(ctg1a, contig_name) == 0 || strcmp(ctg1b, contig_name) == 0)) {
      found++;
    } else if (len == len2 && (strcmp(ctg2a, contig_name) == 0 ||
                               strcmp(ctg2b, contig_name) == 0)) {
      found++;
    }
  }
  return found == 2;
}

/** return true if allele is of ATGC */
static int is_ATGC(const char* s) {
  char* p = (char*)s;
  if (*p == 0) return 0;
  while (*p != 0) {
    if (strchr("ATGCatgc", *p) == NULL) return 0;
    p++;
  }

  return 1;
}

/** insert a new hyperlink in the url table, ignore duplicates */
static void InsertHyperLink(const char* database, const char* label,
                            const char* url) {
  unsigned int i;
  for (i = 0; i < TableNRows(args.hyperlinksTable); ++i) {
    CellPtr cell = TableAt(args.hyperlinksTable, 2, i);
    if (strcmp(CellCStr(cell), url) == 0) return;
  }
  RowPtr row = TableNewRow(args.hyperlinksTable);
  RowSetText(row, 0, database);
  RowSetText(row, 1, label);
  RowSetText(row, 2, url);
}

const char* about(void) {
  return "Convert VCF to tables in the terminal.\n"
         "Author Pierre Lindenbaum PhD. Institut-du-Thorax. U1087. "
         "Nantes/France\n"
         "Options:\n"
         "   -h|--help (string)  help (this screen).\n"
         "   -x|--hide (string)  comma separated list of features to hide:\n"
         "                       HOM_REF or RR : genotypes with REF allele only\n"
         "                       HET or AR : heterozygous genotypes\n"
         "                       NO_CALL or MISSING : missing genotypes\n"
         "                       CSQ or VEP : VEP table\n"
         "                       SPLICEAI : SPLICEAI table\n"
         "                       ANN or SNPEFF : SNPEFF table\n"
         "                       LOF: SNPEFF LOF table\n"
         "                       VC: general table\n"
         "                       INFO: INFO table\n"
         "                       GT: Genotype table\n"
         "                       GTTYPES: Genotype count table\n"
         "                       URL: hyperlink table\n"
         "\nExample:\n"
         "$ wget -O - "
         "'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
         "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
         "ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased."
         "vcf.gz' |\\\n\tbcftools +vcf2table -i 'AC<10' -- --hide "
         "'HOM_REF,INFO,NO_CALL' \n"
         "(...)\n"
         "<<< 22:10714247:C (n. 446)\n"
         "\n"
         "# Variant\n"
         "+--------+----------+\n"
         "| KEY    | VALUE    |\n"
         "+--------+----------+\n"
         "| CHROM  | 22       |\n"
         "| POS    | 10714247 |\n"
         "| ID     | .        |\n"
         "| REF    | C        |\n"
         "| ALT    | G        |\n"
         "| QUAL   | .        |\n"
         "| FILTER | PASS     |\n"
         "+--------+----------+\n"
         "\n"
         "# GENOTYPE TYPES\n"
         "+-----------+-------+----------+\n"
         "| Type      | Count | %        |\n"
         "+-----------+-------+----------+\n"
         "| REF only  | 2545  | 99.8823  |\n"
         "| HET       | 3     | 0.117739 |\n"
         "+-----------+-------+----------+\n"
         "\n"
         "# GENOTYPES\n"
         "+---------+-------+-----+\n"
         "| SAMPLE  | GTYPE | GT  |\n"
         "+---------+-------+-----+\n"
         "| HG03136 | HET   | 1|0 |\n"
         "| HG03171 | HET   | 0|1 |\n"
         "| HG03270 | HET   | 0|1 |\n"
         "+---------+-------+-----+\n"
         "\n"
         ">>> 22:10714247:C (n. 446)\n"

      ;
}

/*
    Called once at startup, it initializes local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char** argv, bcf_hdr_t* hdr_in, bcf_hdr_t* out) {
  int c;
  memset((void*)&args, 0, sizeof(args_t));
  args.header = hdr_in;
  args.out = stdout;
  // initialize table that will not change
  args.annTable =
      TableNewStr("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                  "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
                  "Rank", "HGVS.c", "HGVS.p", "cDNA.pos/length",
                  "CDS.pos/length", "AA.pos/length", "Distance,Message", NULL);
  args.spliceaiTable =
      TableNewStr("ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL",
                  "DP_AG", "DP_AL", "DP_DG", "DP_DL", NULL);
  args.infoTable = TableNewStr("KEY", "IDX", "VALUE", NULL);
  args.vcTable = TableNewStr("KEY", "VALUE", NULL);
  args.lofTable =
      TableNewStr("Gene_Name", "Gene_ID", "Number_of_transcripts_in_gene",
                  "Percent_of_transcripts_affected", NULL);
  args.gtypeTable = TableNewStr("Type", "Count", "%", NULL);
  args.hyperlinksTable = TableNewStr("DB", "" /* empty/misc */, "URL", NULL);

  static struct option loptions[] =
        {
        {"help", required_argument, NULL, 'h'},
        {"hide", required_argument, NULL, 'x'},
        {0, 0, 0, 0}
        };

  while ((c = getopt_long(argc, argv, "hx:", loptions, NULL)) >= 0) {
    switch (c) {
      case 'x': {
        int i;
        StringListPtr hide = StringListMake(optarg, ',');
        for (i = 0; i < hide->size; ++i) {
          const char* hidden = StringListAt(hide, i);
          if (strcasecmp(hidden, "HOM_REF") == 0 ||
              strcasecmp(hidden, "RR") == 0) {
            args.hide_HOM_REF = 1;
          } else if (strcasecmp(hidden, "NO_CALL") == 0 ||
                     strcasecmp(hidden, "MISSING") == 0) {
            args.hide_NO_CALL = 1;
          } else if (strcasecmp(hidden, "HOM_VAR") == 0 ||
                     strcasecmp(hidden, "AA") == 0) {
            args.hide_HOM_VAR = 1;
          } else if (strcasecmp(hidden, "HET") == 0 ||
                     strcasecmp(hidden, "AR") == 0) {
            args.hide_HET = 1;
          } else if (strcasecmp(hidden, "OTHER") == 0) {
            args.hide_OTHER = 1;
          } else if (strcasecmp(hidden, "ANN") == 0 ||
                     strcasecmp(hidden, "SNPEFF") == 0) {
            args.hide_ANN_table = 1;
          } else if (strcasecmp(hidden, "CSQ") == 0 ||
                     strcasecmp(hidden, "VEP") == 0) {
            args.hide_VEP_table = 1;
          } else if (strcasecmp(hidden, "BCSQ") == 0 ||
                     strcasecmp(hidden, "BCFTOOLS") == 0) {
            args.hide_BCSQ_table = 1;
          } else if (strcasecmp(hidden, "SPLICEAI") == 0) {
            args.hide_SPLICEAI_table = 1;
          } else if (strcasecmp(hidden, "INFO") == 0) {
            args.hide_INFO_table = 1;
          } else if (strcasecmp(hidden, "VC") == 0) {
            args.hide_VC_table = 1;
          } else if (strcasecmp(hidden, "LOF") == 0) {
            args.hide_LOF_table = 1;
          } else if (strcasecmp(hidden, "GT") == 0 ||
                     strcasecmp(hidden, "GENOTYPES") == 0) {
            args.hide_GT_table = 1;
          } else if (strcasecmp(hidden, "GTTYPES") == 0) {
            args.hide_GTTYPE_table = 1;
          } else if (strcasecmp(hidden, "URL") == 0 ||
                     strcasecmp(hidden, "URLS") == 0) {
            args.hide_links = 1;
          }
        }
        StringListFree(hide);
        break;
      }
      case 'h':
        fputs(about(), stdout);
        exit(EXIT_SUCCESS);
        break;
      case '?':
      default:
        error("wrong arguments. Use option --help to get help.\n");
        exit(EXIT_FAILURE);
        break;
    }
  }

  if (!isatty(fileno((FILE*)stdout))) {
    args.ascii = 1;
  }

  if (args.ascii == 0) {
    if (setlocale(LC_CTYPE, "") == NULL) {
      fprintf(stderr, "setlocale failed. Switching to ascii\n");
      args.ascii = 1;
    }
  }
  /* guess the build by finding signature of chromosome/length ?*/
  if (findContigs(hdr_in, "1", 249250621, "2", 243199373)) {
    args.build = human_hg19;
  } else if (findContigs(hdr_in, "1", 248956422, "2", 242193529)) {
    args.build = human_hg38;
  } else if (findContigs(hdr_in, "RF01", 3302, "RF02", 2687)) {
    args.build = rotavirus_rf;
  } else {
    args.build = undefined;
  }

  /** find INFO/CSQ and decode it */
  bcf_hrec_t* hrec = bcf_hdr_get_hrec(hdr_in, BCF_HL_INFO, NULL, "CSQ", NULL);
  if (hrec != NULL) {
    int ret = bcf_hrec_find_key(hrec, "Description");
    char* format = ret < 0 ? NULL : strstr(hrec->vals[ret], "Format: ");
    if (format != NULL) {
      format += 8;
      char* vep_format = strdup(format);
      // remove trailing quote
      if (vep_format[strlen(vep_format) - 1] == '"') {
        vep_format[strlen(vep_format) - 1] = 0;
      }
      args.vepTokens = StringListMake(vep_format, '|');
      free(vep_format);
    }
  }
  /** find INFO/BCSQ */
  hrec = bcf_hdr_get_hrec(hdr_in, BCF_HL_INFO, NULL, "BCSQ", NULL);
  if (hrec != NULL) {
    int ret = bcf_hrec_find_key(hrec, "Description");
    char* format = ret < 0 ? NULL : strstr(hrec->vals[ret], "Format: ");
    if (format != NULL) {
      format += 8;
      args.bcsqTokens = StringListMake(format, '|');
    }
  }

  return 1;  // suppress VCF/BCF header
}

#define PRINT_HEADER                                            \
  switch (args.build) {                                         \
    case human_hg19:                                            \
      fputs(" GRCh37 : ", args.out);                            \
      break;                                                    \
    case human_hg38:                                            \
      fputs(" GRCh38 : ", args.out);                            \
      break;                                                    \
    case rotavirus_rf:                                          \
      fputs(" Rotavirus : ", args.out);                         \
      break;                                                    \
    default:                                                    \
      break;                                                    \
  }                                                             \
  fprintf(args.out, " %s:%s:%s (n. %ld)\n", tokens->strings[0], \
          tokens->strings[1], tokens->strings[3], args.n_variants)

/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t* process(bcf1_t* v) {
  TablePtr vepTable = NULL;
  TablePtr bcsqTable = NULL;
  hts_pos_t variant_end = v->pos + v->rlen;

  unsigned int i;
  args.n_variants++;
  /* instead of re-inventing the wheel: the conversion of v to text, let's use
   * htslib/vcf_format to convert the whole line to string */
  kstring_t vcf_line = KS_INITIALIZE;
  vcf_format(args.header, v, &vcf_line);
  // remove last CR/LF
  if (vcf_line.s[vcf_line.l - 1] == '\n') {
    vcf_line.s[vcf_line.l - 1] = 0;
    vcf_line.l--;
  }

  /* split the VCF line into a list of string */
  StringListPtr tokens = StringListMake(vcf_line.s, '\t');
  /* split the ALT alleles */
  StringListPtr alt_alleles = StringListMake(StringListAt(tokens, 4), ',');

  fputs("<<<", args.out);
  PRINT_HEADER;
  fputc('\n', args.out);

  RowPtr row = TableNewRow(args.vcTable);
  CellSetText(RowAt(row, 0), "CHROM");
  CellSetText(RowAt(row, 1), StringListAt(tokens, 0));

  row = TableNewRow(args.vcTable);
  CellSetText(RowAt(row, 0), "POS");
  CellSetText(RowAt(row, 1), StringListAt(tokens, 1));

  if (v->pos + 1 != variant_end) {
    row = TableNewRow(args.vcTable);
    CellSetText(RowAt(row, 0), "end");
    CellSetLL(RowAt(row, 1), variant_end);

    row = TableNewRow(args.vcTable);
    CellSetText(RowAt(row, 0), "length");
    CellSetLL(RowAt(row, 1), variant_end - v->pos);
  }

  row = TableNewRow(args.vcTable);
  CellSetText(RowAt(row, 0), "ID");
  CellSetText(RowAt(row, 1), StringListAt(tokens, 2));

  row = TableNewRow(args.vcTable);
  CellSetText(RowAt(row, 0), "REF");
  CellSetText(RowAt(row, 1), StringListAt(tokens, 3));

  row = TableNewRow(args.vcTable);
  CellSetText(RowAt(row, 0), "ALT");
  CellSetText(RowAt(row, 1), StringListAt(tokens, 4));

  row = TableNewRow(args.vcTable);
  CellSetText(RowAt(row, 0), "QUAL");
  CellSetText(RowAt(row, 1), tokens->strings[5]);

  row = TableNewRow(args.vcTable);
  CellSetText(RowAt(row, 0), "FILTER");
  CellSetText(RowAt(row, 1), StringListAt(tokens, 6));
  if (strcmp(StringListAt(tokens, 6), ".") == 0 ||
      strcmp(StringListAt(tokens, 6), "PASS") == 0) {
    RowAt(row, 1)->color = COLOR_GREEN;
  } else {
    RowAt(row, 1)->color = COLOR_RED;
  }

  if (!args.hide_VC_table) {
    fprintf(args.out, "# Variant\n");
    TablePrint(args.vcTable, &args);
  }

  /** fill URL */
  if (!args.hide_links) {
    hts_pos_t pos1 = v->pos + 1;
    kstring_t url = KS_INITIALIZE;

    if ((args.build == human_hg19 || args.build == human_hg38)) {
      ks_clear(&url);
      kputs("https://genome.ucsc.edu/cgi-bin/hgTracks?db=", &url);
      kputs((args.build == human_hg38 ? "hg38" : "hg19"), &url);
      kputs("&position=", &url);
      kputs(StringListAt(tokens, 0), &url);
      kputs("%3A", &url);
      kputll(pos1, &url);
      kputc('-', &url);
      kputll(variant_end, &url);
      InsertHyperLink("UCSC", "", ks_str(&url));
    }

    for (i = 0; i < alt_alleles->size; ++i) {
      if ((args.build == human_hg19 || args.build == human_hg38) &&
          is_ATGC(StringListAt(tokens, 3)) &&
          is_ATGC(StringListAt(alt_alleles, i))) {
        ks_clear(&url);
        kputs("https://gnomad.broadinstitute.org/variant/", &url);
        kputs(StringListAt(tokens, 0), &url);
        kputc('-', &url);
        kputll(pos1, &url);
        kputc('-', &url);
        kputs(StringListAt(tokens, 3), &url);
        kputc('-', &url);
        kputs(StringListAt(alt_alleles, i), &url);
        kputs((args.build == human_hg38 ? "?dataset=gnomad_r4"
                                        : "?dataset=gnomad_r2_1"),
              &url);
        InsertHyperLink("Gnomad", StringListAt(alt_alleles, i), ks_str(&url));

        ks_clear(&url);
        kputs("https://spliceailookup.broadinstitute.org/#variant=", &url);
        kputs(StringListAt(tokens, 0), &url);
        kputc('-', &url);
        kputll(pos1, &url);
        kputc('-', &url);
        kputs(StringListAt(tokens, 3), &url);
        kputc('-', &url);
        kputs(StringListAt(alt_alleles, i), &url);
        kputs((args.build == human_hg38 ? "?hg=38" : "?hg=19"), &url);
        InsertHyperLink("SpliceAI", StringListAt(alt_alleles, i), ks_str(&url));
      }
      if (args.build == human_hg38 && is_ATGC(StringListAt(tokens, 3)) &&
          is_ATGC(StringListAt(alt_alleles, i))) {
        ks_clear(&url);
        kputs("https://genetics.opentargets.org/variant/", &url);
        kputs(StringListAt(tokens, 0), &url);
        kputc('_', &url);
        kputll(pos1, &url);
        kputc('_', &url);
        kputs(StringListAt(tokens, 3), &url);
        kputc('_', &url);
        kputs(StringListAt(alt_alleles, i), &url);
        InsertHyperLink("OpenTargets", StringListAt(alt_alleles, i),
                        ks_str(&url));

        ks_clear(&url);
        kputs("https://afb.ukbiobank.ac.uk/variant/", &url);
        kputs(StringListAt(tokens, 0), &url);
        kputc('_', &url);
        kputll(pos1, &url);
        kputc('_', &url);
        kputs(StringListAt(tokens, 3), &url);
        kputc('_', &url);
        kputs(StringListAt(alt_alleles, i), &url);
        InsertHyperLink("AF.ukbiobank", StringListAt(alt_alleles, i),
                        ks_str(&url));

        ks_clear(&url);
        kputs("https://genebe.net/variant/hg38/", &url);
        kputs(StringListAt(tokens, 0), &url);
        kputc('-', &url);
        kputll(pos1, &url);
        kputc('-', &url);
        kputs(StringListAt(tokens, 3), &url);
        kputc('-', &url);
        kputs(StringListAt(alt_alleles, i), &url);
        InsertHyperLink("Genebe", StringListAt(alt_alleles, i), ks_str(&url));
      }
    }

    ks_free(&url);
  }

  /* parse values in the INFO column */
  if (tokens->size > 7 && strcmp(StringListAt(tokens, 7), ".") != 0) {
    /* split INFO by semicolon */
    StringListPtr infos = StringListMake(StringListAt(tokens, 7), ';');
    for (i = 0; i < infos->size; i++) {
      unsigned int j;
      const char* info = StringListAt(infos, i);
      char* eq = strchr(info, '=');
      if (eq == NULL || eq == info) continue;
      /* split multiple values for this info using commas */
      StringListPtr values = StringListMake(eq + 1, ',');
      for (j = 0; j < values->size; j++) {
        // skip CSQ
        if (args.vepTokens != NULL && strncmp(info, "CSQ=", 4) == 0) {
          if (args.hide_VEP_table) continue;
          unsigned int k;
          // build VEP table if needed
          if (vepTable == NULL) {
            vepTable = TableNew(0);
            for (k = 0; k < args.vepTokens->size; ++k) {
              TableAppendColumn(vepTable, StringListAt(args.vepTokens, k));
            }
          }
          // fill VEP table
          row = TableNewRow(vepTable);
          StringListPtr veps = StringListMake(StringListAt(values, j), '|');
          for (k = 0; k < args.vepTokens->size && k < veps->size; ++k) {
            CellSetText(RowAt(row, k), StringListAt(veps, k));
          }
          StringListFree(veps);
          continue;
        }

        // skip BCSQ
        if (args.bcsqTokens != NULL && strncmp(info, "BCSQ=", 5) == 0) {
          if (args.hide_BCSQ_table) continue;
          unsigned int k;
          // build BCSQ table if needed
          if (bcsqTable == NULL) {
            bcsqTable = TableNew(0);
            for (k = 0; k < args.bcsqTokens->size; ++k) {
              TableAppendColumn(bcsqTable, StringListAt(args.bcsqTokens, k));
            }
          }

          // fill BCSQ table
          row = TableNewRow(bcsqTable);
          StringListPtr bcsq = StringListMake(StringListAt(values, j), '|');
          for (k = 0; k < args.bcsqTokens->size && k < bcsq->size; ++k) {
            CellSetText(RowAt(row, k), StringListAt(bcsq, k));
          }
          StringListFree(bcsq);
          continue;
        }

        // skip SNPEFF/ANN
        if (strncmp(info, "ANN=", 4) == 0) {
          if (args.hide_ANN_table) continue;
          unsigned int k;
          // fill ANN table
          row = TableNewRow(args.annTable);
          StringListPtr ann = StringListMake(StringListAt(values, j), '|');
          for (k = 0; k < TableNCols(args.annTable) && k < ann->size; ++k) {
            RowSetText(row, k, StringListAt(ann, k));
          }
          StringListFree(ann);
          continue;
        }

        // skip SNPEFF/LOF
        if (strncmp(info, "LOF=", 4) == 0) {
          if (args.hide_LOF_table) continue;
          unsigned int k;
          char* copy = strdup(StringListAt(values, j));
          ASSERT_NOT_NULL(copy);
          // remove first & last char
          if (copy[0] == '(') memmove((void*)&copy[0], &copy[1], strlen(copy));
          if (copy[strlen(copy) - 1] == ')') copy[strlen(copy) - 1] = 0;

          // fill ANN table
          row = TableNewRow(args.lofTable);
          StringListPtr ann = StringListMake(copy, '|');
          for (k = 0; k < TableNCols(args.lofTable) && k < ann->size; ++k) {
            RowSetText(row, k, StringListAt(ann, k));
          }
          StringListFree(ann);
          free(copy);
          continue;
        }

        // skip SpliceAI
        if (strncmp(info, "SpliceAI=", 4) == 0) {
          if (args.hide_SPLICEAI_table) continue;
          unsigned int k;
          // fill ANN table
          row = TableNewRow(args.spliceaiTable);
          StringListPtr spliceai = StringListMake(StringListAt(values, j), '|');
          for (k = 0; k < TableNCols(args.spliceaiTable) && k < spliceai->size;
               ++k) {
            RowSetText(row, k, StringListAt(spliceai, k));
          }
          StringListFree(spliceai);
          continue;
        }

        row = TableNewRow(args.infoTable);
        CellAppendTextN(RowAt(row, 0), info, eq - info);
        if (values->size > 1) CellSetD(RowAt(row, 1), (int)(j + 1));
        RowSetText(row, 2, values->strings[j]);
      }
      StringListFree(values);
    }
    if (!args.hide_INFO_table && TableNRows(args.infoTable) > 0) {
      fprintf(args.out, "# INFO\n");
      TablePrint(args.infoTable, &args);
    }
    StringListFree(infos);
  }

  if (TableNRows(args.hyperlinksTable) > 0) {
    fprintf(args.out, "# HYPERLINKS\n");
    TablePrint(args.hyperlinksTable, &args);
  }

  if (vepTable != NULL && TableNRows(vepTable) > 0) {
    fprintf(args.out, "# VEP/CSQ\n");
    TableRemoveEmptyColumns(vepTable);
    TablePrint(vepTable, &args);
  }

  if (bcsqTable != NULL && TableNRows(bcsqTable) > 0) {
    fprintf(args.out, "# BCSQ\n");
    TableRemoveEmptyColumns(bcsqTable);
    TablePrint(bcsqTable, &args);
  }

  if (TableNRows(args.annTable) > 0) {
    fprintf(args.out, "# ANN/SNPEFF\n");
    // no keep it inummutable TableRemoveEmptyColumns(args.annTable);
    TablePrint(args.annTable, &args);
  }

  if (TableNRows(args.lofTable) > 0) {
    fprintf(args.out, "# LOF\n");
    TablePrint(args.lofTable, &args);
  }

  if (TableNRows(args.spliceaiTable) > 0) {
    fprintf(args.out, "# SpliceAI\n");
    // no keep it inummutable TableRemoveEmptyColumns(args.annTable);
    TablePrint(args.spliceaiTable, &args);
  }

  // is there any genotype here ?
  if (tokens->size > 9) {
    int count_hom_ref = 0;
    int count_het = 0;
    int count_hom_var = 0;
    int count_missing = 0;
    int count_other = 0;
    // column for genotype
    int gt_col = -1;
    // column for filter FT
    int ft_col = -1;
    StringListPtr formats = StringListMake(tokens->strings[8], ':');
    TablePtr genotypeTable = TableNewStr("SAMPLE", NULL);
    TableAppendColumn(genotypeTable, "GTYPE");

    for (i = 0; i < formats->size; i++) {
      TableAppendColumn(genotypeTable, formats->strings[i]);
      if (strcmp("GT", formats->strings[i]) == 0) gt_col = (int)i;
      if (strcmp("FT", formats->strings[i]) == 0) ft_col = (int)i;
    }

    for (i = 9; i < tokens->size; i++) {
      kstring_t gtype_name = KS_INITIALIZE;
      int count_allele_0 = 0;
      int count_allele_1 = 0;
      int count_allele_missing = 0;
      int count_allele_other = 0;
      int print_it = 1;
      // split Genotype components
      StringListPtr values = StringListMake(tokens->strings[i], ':');
      const char* color = COLOR_BLACK;
      unsigned int j;
      if (gt_col != -1 && gt_col < values->size) {
        // clone the GT value
        char* gt = strdup(StringListAt(values, gt_col));
        // remove phasing
        for (j = 0; gt[j] != 0; j++) {
          if (gt[j] == '|') gt[j] = '/';
        }
        // split the alleles in the GT
        StringListPtr alleles = StringListMake(gt, '/');
        for (j = 0; j < alleles->size; ++j) {
          char* allele = alleles->strings[j];
          if (strcmp(allele, "0") == 0)
            count_allele_0++;
          else if (strcmp(allele, "1") == 0)
            count_allele_1++;
          else if (strcmp(allele, ".") == 0)
            count_allele_missing++;
          else
            count_allele_other++;
        }

        if (alleles->size == 2) {
          if (count_allele_0 == 0 && count_allele_1 == 0 &&
              count_allele_other == 0) {
            kputs("NO_CALL", &gtype_name);
            if (args.hide_NO_CALL) print_it = 0;
            count_missing++;
          } else if (count_allele_0 == 2) {
            kputs("HOM_REF", &gtype_name);
            color = COLOR_GREEN;
            if (args.hide_HOM_REF) print_it = 0;
            count_hom_ref++;
          } else if (count_allele_missing == 0 &&
                     strcmp(StringListAt(alleles, 0),
                            StringListAt(alleles, 1)) == 0) {
            kputs("HOM_VAR", &gtype_name);
            color = COLOR_RED;
            if (args.hide_HOM_VAR) print_it = 0;
            count_hom_var++;
          } else if (count_allele_missing == 0 &&
                     strcmp(StringListAt(alleles, 0),
                            StringListAt(alleles, 1)) != 0) {
            kputs("HET", &gtype_name);
            color = COLOR_CYAN;
            count_het++;
            if (args.hide_HET) print_it = 0;
          } else {
            if (args.hide_OTHER) print_it = 0;
            count_other++;
          }
        } else if (alleles->size == 1) {
          if (count_allele_0 == 1) {
            kputs("REF", &gtype_name);
            color = COLOR_GREEN;
            if (args.hide_HOM_REF) print_it = 0;
            count_hom_ref++;
          } else if (count_allele_1 == 1) {
            kputs("ALT", &gtype_name);
            color = COLOR_RED;
            count_hom_var++;
          } else if (count_allele_missing == 1) {
            kputs("NO_CALL", &gtype_name);
            if (args.hide_NO_CALL) print_it = 0;
            count_missing++;
          } else {
            if (args.hide_OTHER) print_it = 0;
            count_other++;
          }
        } else {
          if (count_allele_0 == alleles->size) {
            kputs("HOM_REF", &gtype_name);
            color = COLOR_GREEN;
            if (args.hide_HOM_REF) print_it = 0;
            count_hom_ref++;
          } else if (count_allele_1 == alleles->size) {
            kputs("HOM_VAR", &gtype_name);
            color = COLOR_RED;
            if (args.hide_HOM_VAR) print_it = 0;
            count_hom_ref++;
          } else if (count_allele_missing == alleles->size) {
            kputs("NO_CALL", &gtype_name);
            if (args.hide_NO_CALL) print_it = 0;
            count_missing++;
          } else {
            if (args.hide_OTHER) print_it = 0;
            count_other++;
          }
        }
        StringListFree(alleles);
        free(gt);
      }

      if (print_it && !args.hide_GT_table) {
        row = TableNewRow(genotypeTable);
        CellSetText(RowAt(row, 0), args.header->samples[i - 9]);
        CellSetText(RowAt(row, 1), gtype_name.s);
        RowAt(row, 1)->color = color;

        for (j = 0; j < values->size; j++) {
          CellSetText(RowAt(row, j + 2), StringListAt(values, j));
          // color genotype if FORMAT/FT
          if (ft_col == j) {
            if (strcmp(StringListAt(values, j), "PASS") == 0 ||
                strcmp(StringListAt(values, j), ".") == 0) {
              RowAt(row, j + 2)->color = COLOR_GREEN;
            } else {
              RowAt(row, j + 2)->color = COLOR_RED;
            }
          }
        }
      }
      StringListFree(values);
      ks_free(&gtype_name);
    }
#define ADD_GT(LABEL, COUNT)                                   \
  if (COUNT > 0 && total > 0) {                                \
    row = TableNewRow(args.gtypeTable);                        \
    RowSetText(row, 0, LABEL);                                 \
    CellSetLL(RowAt(row, 1), COUNT);                           \
    CellSetD(RowAt(row, 2), 100.0 * (COUNT / ((float)total))); \
  }
    if (!args.hide_GTTYPE_table) {
      int total = count_hom_ref + count_het + count_hom_var + count_missing +
                  count_other;
      ADD_GT("REF only ", count_hom_ref)
      ADD_GT("HET", count_het)
      ADD_GT("ALT only", count_hom_var)
      ADD_GT("MISSING", count_missing)
      ADD_GT("OTHER", count_other)

      if (TableNRows(args.gtypeTable) > 0) {
        fprintf(args.out, "# GENOTYPE TYPES\n");
        TablePrint(args.gtypeTable, &args);
      }
    }
#undef ADD_GT

    if (!args.hide_GT_table && TableNRows(genotypeTable) > 0) {
      fprintf(args.out, "# GENOTYPES\n");
      TablePrint(genotypeTable, &args);
    }
    TableFree(genotypeTable);
    StringListFree(formats);
  }

  fputs(">>>", args.out);
  PRINT_HEADER;

  fputc('\n', args.out);

  /** final cleanup */
  ks_free(&vcf_line);
  StringListFree(tokens);
  StringListFree(alt_alleles);
  TableFree(bcsqTable);
  TableFree(vepTable);
  TableClear(args.annTable);
  TableClear(args.lofTable);
  TableClear(args.spliceaiTable);
  TableClear(args.infoTable);
  TableClear(args.vcTable);
  TableClear(args.gtypeTable);
  TableClear(args.hyperlinksTable);
  return NULL; /* suppress bcf output */
}

void destroy(void) {
  StringListFree(args.vepTokens);
  StringListFree(args.bcsqTokens);
  TableFree(args.spliceaiTable);
  TableFree(args.annTable);
  TableFree(args.lofTable);
  TableFree(args.infoTable);
  TableFree(args.vcTable);
  TableFree(args.gtypeTable);
  TableFree(args.hyperlinksTable);
}
