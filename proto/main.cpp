#include <iostream>
#include <string>

#define BLOSUM62_ROWS 33
#define BLOSUM62_COLS 80

// static const std::string blosum62 =
//   "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
//   "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n"
//   "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n"
//   "N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n"
//   "D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n"
//   "C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n"
//   "Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n"
//   "E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n"
//   "G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n"
//   "H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n"
//   "I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n"
//   "L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n"
//   "K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n"
//   "M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n"
//   "F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n"
//   "P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n"
//   "S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n"
//   "T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n"
//   "W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n"
//   "Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n"
//   "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n"
//   "B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n"
//   "Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n"
//   "X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n"
//   "* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

static const char blosum62[] = {
#if 0
  "#  Matrix made by matblas from blosum62.iij\n",
  "#  * column uses minimum score\n",
  "#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units\n",
  "#  Blocks Database = /data/blocks_5.0/blocks.dat\n",
  "#  Cluster Percentage: >= 62\n",
  "#  Entropy =   0.6979, Expected =  -0.5209\n",
#endif
  "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
  "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n"
  "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n"
  "N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n"
  "D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n"
  "C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n"
  "Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n"
  "E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n"
  "G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n"
  "H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n"
  "I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n"
  "L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n"
  "K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n"
  "M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n"
  "F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n"
  "P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n"
  "S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n"
  "T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n"
  "W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n"
  "Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n"
  "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n"
  "B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n"
  "Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n"
  "X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n"
  "* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n"
};

#ifndef int2
typedef int int2[2];
#endif

int getResName(char *resName) {
    int name = 0;
    for (int i = 0; i < 3; i++){
        name = (name << 8) | resName[i];
    }
    return name;
}

const char *ParseNextLine(const char *p)
{
  char ch;
  const char mask = -16;        /* 0xF0 */
  while((mask & p[0]) && (mask & p[1]) && (mask & p[2]) && (mask & p[3]))       /* trusting short-circuit to avoid overrun */
    p += 4;
  while((ch = *p)) {
    p++;
    if(ch == 0xD) {             /* Mac or PC */
      if((*p) == 0xA)           /* PC */
        return p + 1;
      return p;
    } else if(ch == 0xA) {      /* Unix */
      return p;
    }
  }
  return p;
}

const char *ParseWordCopy(char *q, const char *p, int n)
{                               /* word copy */
  while(*p) {
    if((*p == 0xD) || (*p == 0xA))      /* don't skip end of lines */
      break;
    if(*p <= 32)
      p++;
    else
      break;
  }
  while(*p) {
    if(*p <= 32)
      break;
    if(!n) {
      while(*p > 32)            /* finish scanning word, but don't copy into field */
        p++;
      break;
    }
    if((*p == 0xD) || (*p == 0xA))      /* don't copy end of lines */
      break;
    *(q++) = *(p++);
    n--;
  }
  *q = 0;
  return p;
}

class Protein {
public:
    int *rnames;
    int *rcodes;
    int seq_size;
    Protein(int size);
    ~Protein();
    void addResidues(char residues[][4]);
    void printValues();
};

Protein::Protein(int size) {
    rnames = (int *) calloc(size, sizeof(int));
    rcodes = (int *) calloc(size, sizeof(int));
    seq_size = size;
}

Protein::~Protein() {
    free(rnames);
    free(rcodes);
    rnames = NULL;
    rcodes = NULL;
}
void Protein::addResidues(char residues[][4]) {
    int a, b, c;
    // get integral values for the residue names
    for(a = 0; a < seq_size; a++) {
        b = 0;
        for(c = 0; c < 3; c++) {
            b = (b << 8) | residues[a][c];
        }
        rnames[a] = b; // store
    }
}

void Protein::printValues() {
    printf("\nrname: ");
    for(int i= 0; i < seq_size; i++) {
        printf("%i, ", rnames[i]);
    }
    printf("\nrcode: ");
    for(int i= 0; i < seq_size; i++) {
        printf("%i, ", rcodes[i]);
    }
}

void MatchResidueToCode(Protein &p) {
    int ok = true;
    int a, b, c;
    int found;
    int *trg;
    char res[][4] = {

        "HOH", "O",                 /* water */

        /* using numbers to prevent nucleic acids from getting confounded with
        protein scores */

        "A", "1", 
        "DA", "1",
        "ADE", "1",

        "C", "2",
        "DC", "2",
        "CYT", "2",

        "G", "3",
        "DG", "3",
        "GUA", "3",

        "T", "4",
        "DT", "4",
        "THY", "4",

        "U", "4",
        "URA", "4",

        /* these should correspond to the matrix being read */

        "ALA", "A",
        "CYS", "C",
        "ASP", "D",
        "GLU", "E",
        "PHE", "F",

        "GLY", "G",
        "HIS", "H",
        "ILE", "I",
        "LYS", "K",
        "LEU", "L",

        "MET", "M",
        "MSE", "M",                 /* selenomet */

        "ASN", "N",
        "PRO", "P",
        "GLN", "Q",
        "ARG", "R",

        "SER", "S",
        "SEP", "S",                 /* phosphoserine */
        "THR", "T",
        "TPO", "T",                 /* phosphothreonine */
        "VAL", "V",

        "TRP", "W",
        "TYR", "Y",
        "PTR", "Y",                 /* phosphotyrosine */
        "TYS", "Y"                  /* sulfotyrosine */
    };

    const int cNRES = (sizeof res) / 8;
    int rcode[cNRES], rname[cNRES];

    /* get integral values for the residue names */

    for(a = 0; a < cNRES; a++) {
        b = 0;
        for(c = 0; c < 3; c++) {
        b = (b << 8) | res[a * 2][c];
        }
        rname[a] = b;
        rcode[a] = res[a * 2 + 1][0];
    }

    int n = p.seq_size;

    // get one letter code for rcodes
    for(b = 0; b < n; b++) {
        found = 0;
        for(a = 0; a < cNRES; a++) {
            if(rname[a] == p.rnames[b]) {
                found = true;
                p.rcodes[b] = rcode[a];
                break;
            }
        }
        if(!found) {
            // codes for unknown residues are three byte (mask 0xFFFFFF80)
            p.rcodes[b] = p.rcodes[b] << 8;
        }
    }
}



int main(int, char**) {
    char res_seq_1[][4] =  {"GLY", "PRO", "SER", "GLN", "ALA", "ILE", "LYS", "CYS", "VAL", "VAL", "VAL", "GLY", "ASP"};
    char res_seq_2[][4] =  {"MET", "ARG", "GLU", "TYR", "LYS", "VAL", "VAL", "VAL", "LEU", "GLY", "SER", "GLY", "GLY",};

    Protein p1(sizeof(res_seq_1)/sizeof(res_seq_1[0]));
    Protein p2(sizeof(res_seq_2)/sizeof(res_seq_2[0]));

    int na = p1.seq_size;
    int nb = p2.seq_size;

    for(int i = 0; i < na; i++) {
        p1.addResidues(res_seq_1);
    }
    for(int i = 0; i < nb; i++) {
        p2.addResidues(res_seq_2);
    }

    // MatchNew ======================================================================================
    // int seqSize1 = 167;
    // int seqSize2 = 183;
    float mat[na][nb];
    float smat[128][128];

    /* lay down an 10/-1 identity matrix to cover matches for known
       residues other than amino acids (dna, rna, as 1,2,3,4 etc.)
       these values will be overwritten by the matrix */
    for(int i = 0; i < 128; i++) {
      for(int j = 0; j < 128; j++) {
        smat[i][j] = -1.0F; 
      }
    }
    for(int i = 0; i < 128; i++) {
      smat[i][i] = 10.0F; /* these values will be overwritten by BLOSUM, etc. */
    }
    smat['O']['O'] = -1.0F;
    // printf("%f %f %f %f", smat[0][0], smat[1][1], smat[2][2], smat[3][3]);

    p1.printValues();
    p2.printValues();

    MatchResidueToCode(p1);
    MatchResidueToCode(p2);

    p1.printValues();
    p2.printValues();

    // MatchFile ================================================================================
    int ok = 1;
    std::string buffer;
    const char *p;
    char cc[255];
    char *code = NULL;
    unsigned int x, y;
    int a;
    int n_entry;

    buffer = blosum62;

    // if(ok && !buffer.empty()) {

    /* count codes */

    p = buffer.c_str();
    n_entry = 0;
    while(*p && ok) {
        switch (*p) {
            case '#':
                break;
            default:
                if((*p) > 32)
                n_entry++;
                break;
        }
        p = ParseNextLine(p);
    }

    // if(!n_entry)
    //     ok = false;
    // else {
    code = (char *) calloc(n_entry, sizeof(int));
    // code = (char*) pymol::calloc<int>(n_entry);

    /* read codes */

    p = buffer.c_str();
    n_entry = 0;
    while(*p && ok) {
        switch (*p) {
        case '#':
        break;
        default:
        if((*p) > 32) {
            code[n_entry] = *p;
            n_entry++;
        }
        break;
        }
        p = ParseNextLine(p);
    }

    /* read values */

    p = buffer.c_str();
    while((*p) && ok) {
        switch (*p) {
        case '#':
        break;
        default:
        if((*p) > 32) {
            x = *(p++);
            for(a = 0; a < n_entry; a++) {
            p = ParseWordCopy(cc, p, 255);
            y = (unsigned int) code[a];
            ok = sscanf(cc, "%f", &smat[x][y]);
            }
        }
        break;
        }
        if(!ok)
        break;
        p = ParseNextLine(p);
    }
    // }
    // }
    free(code);
    code = NULL;

    // MatchPreScore =================================================================================
    int b;

    for(a = 0; a < na; a++) {
        std::cout << std::endl;
        for(b = 0; b < nb; b++) {
        // codes for known   residues are one   byte (mask 0x0000007F)
        // codes for unknown residues are three byte (mask 0xFFFFFF80)
        // This allows for exact match of unknown three-letter codes.
        // Fallback for unknown residues with no exact match is 'X'
            int code1 = p1.rcodes[a];
            int code2 = p1.rcodes[b];
            if (code1 & 0xFFFFFF80) {
                if (code1 == code2) {
                mat[a][b] = 5.F;
                continue;
                }
                code1 = 'X';
            }
            if (code2 & 0xFFFFFF80) {
                code2 = 'X';
            }
            mat[a][b] = smat[code1][code2];
            std::cout << mat[a][b] << " ";
        }
    }

    // MatchFile =================================================================================
    PyMOLGlobals *G = I->G;

    int ok = 1;
    std::string buffer;
    const char *p;
    char cc[255];
    char *code = NULL;
    unsigned int x, y;
    int a;
    int n_entry;

    if(fname && fname[0]
  #ifdef _PYMOL_NOPY
      /* if Python is absent, then use the hardcoded copy of BLOSUM62 */
      && !(strcmp(fname, "BLOSUM62") == 0)
  #endif
      ) {
      try {
        buffer = pymol::file_get_contents(fname);
      } catch (...) {
        PRINTFB(G, FB_Match, FB_Errors)
          " Match-Error: unable to open matrix file '%s'.\n", fname ENDFB(G);
        ok = false;
      }
    } else {
      buffer = blosum62;
    }

    if(ok && !buffer.empty()) {

      /* count codes */

      p = buffer.c_str();
      n_entry = 0;
      while(*p && ok) {
        switch (*p) {
        case '#':
          break;
        default:
          if((*p) > 32)
            n_entry++;
          break;
        }
        p = ParseNextLine(p);
      }

      if(!n_entry)
        ok = false;
      else {
        code = (char*) pymol::calloc<int>(n_entry);

        /* read codes */

        p = buffer.c_str();
        n_entry = 0;
        while(*p && ok) {
          switch (*p) {
          case '#':
            break;
          default:
            if((*p) > 32) {
              code[n_entry] = *p;
              n_entry++;
            }
            break;
          }
          p = ParseNextLine(p);
        }

        /* read values */

        p = buffer.c_str();
        while((*p) && ok) {
          switch (*p) {
          case '#':
            break;
          default:
            if((*p) > 32) {
              x = *(p++);
              for(a = 0; a < n_entry; a++) {
                p = ParseWordCopy(cc, p, 255);
                y = (unsigned int) code[a];
                ok = sscanf(cc, "%f", &I->smat[x][y]);
              }
            }
            break;
          }
          if(!ok)
            break;
          p = ParseNextLine(p);
        }
      }
    }
    if(ok) {
      if(!quiet) {
        PRINTFB(G, FB_Match, FB_Details)
          " Match: read scoring matrix.\n" ENDFB(G);
      }
    }
    FreeP(code);
    return (ok);
    return 1;
}
