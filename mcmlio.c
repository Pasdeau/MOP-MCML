/***********************************************************
 *  Input/output of data.
 *  Refinements:
 *   - ShowVersion prototype matches header (const char*)
 *   - CheckChar uses unsigned char to avoid UB and removes
 *     tautological comparisons; keeps "Non-ASCII file" check.
 *   - Ensure trailing newline at EOF.
 ****/

#include "mcml.h"

 /***********************************************************
  * Name list for duplicate file checks.
  ****/
struct NameList {
    char name[STRLEN];
    struct NameList* next;
};
typedef struct NameList NameNode;
typedef NameNode* NameLink;

/***********************************************************
 * Center a string according to the column width.
 ****/
static char* CenterStr(short Wid, char* InStr, char* OutStr)
{
    size_t nspaces = (Wid - (short)strlen(InStr)) / 2;
    if ((short)nspaces < 0) nspaces = 0;

    OutStr[0] = '\0';
    while (nspaces--) strcat(OutStr, " ");
    strcat(OutStr, InStr);
    return OutStr;
}

/***********************************************************
 * Print banner.
 ****/
#define COLWIDTH 80
void ShowVersion(const char* version)
{
    char str[STRLEN];

    CenterStr(COLWIDTH, "MCML - Monte Carlo Simulation of Multi-layered Turbid Media", str);
    puts(str); puts("");

    CenterStr(COLWIDTH, "Lihong Wang, Ph. D.", str); puts(str);
    CenterStr(COLWIDTH, "Steven L. Jacques, Ph. D.", str); puts(str);

    CenterStr(COLWIDTH, "Laser Biology Research Laboratory - Box 17", str); puts(str);
    CenterStr(COLWIDTH, "University of Texas / M.D. Anderson Cancer Center", str); puts(str);
    CenterStr(COLWIDTH, "Houston, TX 77030", str); puts(str);
    CenterStr(COLWIDTH, "Phone: (713) 792-3664", str); puts(str);
    CenterStr(COLWIDTH, "Fax:   (713) 792-3995", str); puts(str);

    CenterStr(COLWIDTH, "e-mails: lihong@laser.mda.uth.tmc.edu", str); puts(str);
    CenterStr(COLWIDTH, "         or slj@laser.mda.uth.tmc.edu", str); puts(str);
    puts("");

    CenterStr(COLWIDTH, (char*)version, str); puts(str);
    puts("\n\n\n\n");
}
#undef COLWIDTH

/***********************************************************
 * Get a filename and open it for reading.
 ****/
FILE* GetFile(char* Fname)
{
    FILE* file = NULL;
    Boolean firsttime = true;

    do {
        if (firsttime && Fname[0] != '\0') {
            firsttime = false;
        }
        else {
            printf("Input filename(or . to exit):");
            if (scanf("%255s", Fname) != 1) exit(1);
            firsttime = false;
        }

        if (strlen(Fname) == 1 && Fname[0] == '.')
            exit(1);

        file = fopen(Fname, "r");
    } while (file == NULL);

    return file;
}

/***********************************************************
 * Kill the ith char and shift left.
 ****/
static void KillChar(size_t i, char* Str)
{
    size_t sl = strlen(Str);
    for (; i < sl; i++) Str[i] = Str[i + 1];
}

/***********************************************************
 * Eliminate non-printing chars (or non-ASCII).
 * Return 1 if any bad char found (as original semantic).
 ****/
Boolean CheckChar(char* Str)
{
    Boolean found = 0;
    size_t sl = strlen(Str);
    size_t i = 0;

    while (i < sl) {
        unsigned char c = (unsigned char)Str[i];
        /* Keep original policy: reject non-ASCII input lines */
        if (c > 127) {
            nrerror("Non-ASCII file\n");
        }
        else if (isprint(c) || isspace(c)) {
            i++;
        }
        else {
            found = 1;
            KillChar(i, Str);
            sl--;
        }
    }
    return found;
}

/***********************************************************
 * Return 1 if comment line or blank line.
 ****/
static Boolean CommentLine(char* Buf)
{
    size_t spn = strspn(Buf, " \t");
    size_t cspn = strcspn(Buf, "#\n");
    return (spn == cspn) ? 1 : 0;
}

/***********************************************************
 * Skip space or comment lines and return a data line only.
 ****/
char* FindDataLine(FILE* File_Ptr)
{
    static char buf[STRLEN];

    buf[0] = '\0';
    do {
        if (fgets(buf, 255, File_Ptr) == NULL) {
            printf("Incomplete data\n");
            buf[0] = '\0';
            break;
        }
        else {
            CheckChar(buf);
        }
    } while (CommentLine(buf));

    return buf;
}

/***********************************************************
 * Read number of runs (skip file version line).
 ****/
short ReadNumRuns(FILE* File_Ptr)
{
    char buf[STRLEN];
    short n = 0;

    FindDataLine(File_Ptr); /* skip version */
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading number of runs\n");
    sscanf(buf, "%hd", &n);
    return n;
}

/***********************************************************
 * Read output file name and format.
 ****/
void ReadFnameFormat(FILE* File_Ptr, InputStruct* In_Ptr)
{
    char buf[STRLEN];
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading file name and format.\n");
    sscanf(buf, "%255s %c", In_Ptr->out_fname, &(In_Ptr->out_fformat));
    if (toupper(In_Ptr->out_fformat) != 'B') In_Ptr->out_fformat = 'A';
}

/***********************************************************
 * Read number of photons.
 ****/
void ReadNumPhotons(FILE* File_Ptr, InputStruct* In_Ptr)
{
    char buf[STRLEN];
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading number of photons.\n");
    sscanf(buf, "%ld", &In_Ptr->num_photons);
    if (In_Ptr->num_photons <= 0) nrerror("Nonpositive number of photons.\n");
}

/***********************************************************
 * Read dz, dr.
 ****/
void ReadDzDr(FILE* File_Ptr, InputStruct* In_Ptr)
{
    char buf[STRLEN];
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading dz, dr.\n");
    sscanf(buf, "%lf%lf", &In_Ptr->dz, &In_Ptr->dr);
    if (In_Ptr->dz <= 0) nrerror("Nonpositive dz.\n");
    if (In_Ptr->dr <= 0) nrerror("Nonpositive dr.\n");
}

/***********************************************************
 * Read nz, nr, na.
 ****/
void ReadNzNrNa(FILE* File_Ptr, InputStruct* In_Ptr)
{
    char buf[STRLEN];
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading number of dz, dr, da's.\n");
    sscanf(buf, "%hd%hd%hd", &In_Ptr->nz, &In_Ptr->nr, &In_Ptr->na);
    if (In_Ptr->nz <= 0) nrerror("Nonpositive number of dz's.\n");
    if (In_Ptr->nr <= 0) nrerror("Nonpositive number of dr's.\n");
    if (In_Ptr->na <= 0) nrerror("Nonpositive number of da's.\n");
    In_Ptr->da = 0.5 * PI / In_Ptr->na;
}

/***********************************************************
 * Read number of layers.
 ****/
void ReadNumLayers(FILE* File_Ptr, InputStruct* In_Ptr)
{
    char buf[STRLEN];
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading number of layers.\n");
    sscanf(buf, "%hd", &In_Ptr->num_layers);
    if (In_Ptr->num_layers <= 0) nrerror("Nonpositive number of layers.\n");
}

/***********************************************************
 * Read PD parameters.
 ****/
void ReadPD(FILE* File_Ptr, InputStruct* In_Ptr) {
    char buf[STRLEN];
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading PD.\n");
    sscanf(buf, "%lf%lf%lf%lf%lf%lf",
        &In_Ptr->PD_Rx, &In_Ptr->PD_Ry, &In_Ptr->PD_Rl,
        &In_Ptr->PD_Tx, &In_Ptr->PD_Ty, &In_Ptr->PD_Tl);
    if (In_Ptr->PD_Rl <= 0) nrerror("Nonpositive PD_Rl.\n");
    if (In_Ptr->PD_Tl <= 0) nrerror("Nonpositive PD_Tl.\n");
}

/***********************************************************
 * Read light source parameters.
 ****/
void ReadLight(FILE* File_Ptr, InputStruct* In_Ptr) {
    char buf[STRLEN];
    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading light.\n");
    sscanf(buf, "%hd%lf%lf%lf",
        &In_Ptr->lightType, &In_Ptr->light_x, &In_Ptr->light_y, &In_Ptr->light_l);
    if (In_Ptr->light_l <= 0) nrerror("Nonpositive light_l.\n");
}

/***********************************************************
 * Read ambient refractive index n.
 ****/
void ReadAmbient(FILE* File_Ptr, LayerStruct* Layer_Ptr, char* side)
{
    (void)side;
    char buf[STRLEN];
    double n;

    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') nrerror("Reading ambient n.\n");

    sscanf(buf, "%lf", &n);
    if (n <= 0) nrerror("Wrong n.\n");
    Layer_Ptr->n = n;
}

/***********************************************************
 * Read one layer.
 ****/
Boolean ReadOneLayer(FILE* File_Ptr, LayerStruct* Layer_Ptr, double* Z_Ptr)
{
    char buf[STRLEN];
    double d, n, mua, mus, g;

    strcpy(buf, FindDataLine(File_Ptr));
    if (buf[0] == '\0') return 1;

    sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &d);
    if (d < 0 || n <= 0 || mua < 0 || mus < 0 || g < 0 || g > 1) return 1;

    Layer_Ptr->n = n;
    Layer_Ptr->mua = mua;
    Layer_Ptr->mus = mus;
    Layer_Ptr->g = g;
    Layer_Ptr->z0 = *Z_Ptr;
    *Z_Ptr += d;
    Layer_Ptr->z1 = *Z_Ptr;

    return 0;
}

/***********************************************************
 * Read layer specs.
 ****/
void ReadLayerSpecs(FILE* File_Ptr, short Num_Layers, LayerStruct** Layerspecs_PP)
{
    char msg[STRLEN];
    short i = 0;
    double z = 0.0;

    *Layerspecs_PP = (LayerStruct*)malloc((unsigned)(Num_Layers + 2) * sizeof(LayerStruct));
    if (!(*Layerspecs_PP)) nrerror("allocation failure in ReadLayerSpecs()");

    ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "top");
    for (i = 1; i <= Num_Layers; i++)
        if (ReadOneLayer(File_Ptr, &((*Layerspecs_PP)[i]), &z)) {
            sprintf(msg, "Error reading %hd of %hd layers\n", i, Num_Layers);
            nrerror(msg);
        }
    ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "bottom");
}

/***********************************************************
 * Compute critical cosines for total internal reflection.
 ****/
void CriticalAngle(short Num_Layers, LayerStruct** Layerspecs_PP)
{
    short i = 0;
    double n1, n2;

    for (i = 1; i <= Num_Layers; i++) {
        n1 = (*Layerspecs_PP)[i].n;
        n2 = (*Layerspecs_PP)[i - 1].n;
        (*Layerspecs_PP)[i].cos_crit0 = n1 > n2 ? sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;

        n2 = (*Layerspecs_PP)[i + 1].n;
        (*Layerspecs_PP)[i].cos_crit1 = n1 > n2 ? sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;
    }
}

/***********************************************************
 * Read all input parameters for one run.
 ****/
void ReadParm(FILE* File_Ptr, InputStruct* In_Ptr)
{
    In_Ptr->Wth = WEIGHT;

    ReadFnameFormat(File_Ptr, In_Ptr);
    ReadNumPhotons(File_Ptr, In_Ptr);
    ReadDzDr(File_Ptr, In_Ptr);
    ReadNzNrNa(File_Ptr, In_Ptr);
    ReadNumLayers(File_Ptr, In_Ptr);

    ReadLayerSpecs(File_Ptr, In_Ptr->num_layers, &In_Ptr->layerspecs);
    CriticalAngle(In_Ptr->num_layers, &In_Ptr->layerspecs);
    ReadPD(File_Ptr, In_Ptr);
    ReadLight(File_Ptr, In_Ptr);
}

/***********************************************************
 * Duplicate output name checks.
 ****/
static Boolean NameInList(char* Name, NameLink List)
{
    while (List != NULL) {
        if (strcmp(Name, List->name) == 0) return 1;
        List = List->next;
    }
    return 0;
}

static void AddNameToList(char* Name, NameLink* List_Ptr)
{
    NameLink list = *List_Ptr;

    if (list == NULL) {
        *List_Ptr = list = (NameLink)malloc(sizeof(NameNode));
        strcpy(list->name, Name);
        list->next = NULL;
    }
    else {
        while (list->next != NULL) list = list->next;
        list->next = (NameLink)malloc(sizeof(NameNode));
        list = list->next;
        strcpy(list->name, Name);
        list->next = NULL;
    }
}

static Boolean FnameTaken(char* fname, NameLink* List_Ptr)
{
    if (NameInList(fname, *List_Ptr)) return 1;
    AddNameToList(fname, List_Ptr);
    return 0;
}

static void FreeFnameList(NameLink List)
{
    NameLink next;
    while (List != NULL) {
        next = List->next;
        free(List);
        List = next;
    }
}

/***********************************************************
 * Validate parameters for each run (filenames duplicate).
 ****/
void CheckParm(FILE* File_Ptr, InputStruct* In_Ptr)
{
    short i_run;
    short num_runs;
    NameLink head = NULL;
    Boolean name_taken;
    char msg[STRLEN];

    num_runs = ReadNumRuns(File_Ptr);
    for (i_run = 1; i_run <= num_runs; i_run++) {
        printf("Checking input data for run %hd\n", i_run);
        ReadParm(File_Ptr, In_Ptr);

        name_taken = FnameTaken(In_Ptr->out_fname, &head);
        if (name_taken)
            sprintf(msg, "file name %s duplicated.\n", In_Ptr->out_fname);

        free(In_Ptr->layerspecs);
        if (name_taken) nrerror(msg);
    }
    FreeFnameList(head);
    rewind(File_Ptr);
}

/***********************************************************
 * Init output arrays for a run.
 ****/
void InitOutputData(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    short nz = In_Parm.nz, nr = In_Parm.nr, na = In_Parm.na, nl = In_Parm.num_layers;

    if (nz <= 0 || nr <= 0 || na <= 0 || nl <= 0)
        nrerror("Wrong grid parameters.\n");

    Out_Ptr->Rsp = 0.0;
    Out_Ptr->Rd = 0.0;
    Out_Ptr->A = 0.0;
    Out_Ptr->Tt = 0.0;

    Out_Ptr->Rd_ra = AllocMatrix(0, nr - 1, 0, na - 1);
    Out_Ptr->Rd_r = AllocVector(0, nr - 1);
    Out_Ptr->Rd_a = AllocVector(0, na - 1);

    Out_Ptr->A_rz = AllocMatrix(0, nr - 1, 0, nz - 1);
    Out_Ptr->A_z = AllocVector(0, nz - 1);
    Out_Ptr->A_l = AllocVector(0, nl + 1);

    Out_Ptr->Tt_ra = AllocMatrix(0, nr - 1, 0, na - 1);
    Out_Ptr->Tt_r = AllocVector(0, nr - 1);
    Out_Ptr->Tt_a = AllocVector(0, na - 1);

    Out_Ptr->OP = AllocMatrix(0, nr - 1, 0, nz - 1);

    Out_Ptr->photonsnbrR = 0;
    Out_Ptr->photonsnbrT = 0;
    Out_Ptr->nbrphotons = 0;
}

/***********************************************************
 * Free output arrays.
 ****/
void FreeData(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    short nz = In_Parm.nz, nr = In_Parm.nr, na = In_Parm.na, nl = In_Parm.num_layers;

    free(In_Parm.layerspecs);

    FreeMatrix(Out_Ptr->Rd_ra, 0, nr - 1, 0, na - 1);
    FreeVector(Out_Ptr->Rd_r, 0, nr - 1);
    FreeVector(Out_Ptr->Rd_a, 0, na - 1);

    FreeMatrix(Out_Ptr->A_rz, 0, nr - 1, 0, nz - 1);
    FreeVector(Out_Ptr->A_z, 0, nz - 1);
    FreeVector(Out_Ptr->A_l, 0, nl + 1);

    FreeMatrix(Out_Ptr->Tt_ra, 0, nr - 1, 0, na - 1);
    FreeVector(Out_Ptr->Tt_r, 0, nr - 1);
    FreeVector(Out_Ptr->Tt_a, 0, na - 1);

    FreeMatrix(Out_Ptr->OP, 0, nr - 1, 0, nz - 1);
}

/***********************************************************
 * Summations.
 ****/
void Sum2DRd(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    short nr = In_Parm.nr, na = In_Parm.na;
    short ir, ia;
    double sum;

    for (ir = 0; ir < nr; ir++) {
        sum = 0.0;
        for (ia = 0; ia < na; ia++) sum += Out_Ptr->Rd_ra[ir][ia];
        Out_Ptr->Rd_r[ir] = sum;
    }

    for (ia = 0; ia < na; ia++) {
        sum = 0.0;
        for (ir = 0; ir < nr; ir++) sum += Out_Ptr->Rd_ra[ir][ia];
        Out_Ptr->Rd_a[ia] = sum;
    }

    sum = 0.0;
    for (ir = 0; ir < nr; ir++) sum += Out_Ptr->Rd_r[ir];
    Out_Ptr->Rd = sum;
}

short IzToLayer(short Iz, InputStruct In_Parm)
{
    short i = 1;
    short num_layers = In_Parm.num_layers;
    double dz = In_Parm.dz;

    while ((Iz + 0.5) * dz >= In_Parm.layerspecs[i].z1 && i < num_layers) i++;
    return i;
}

void Sum2DA(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    short nz = In_Parm.nz, nr = In_Parm.nr;
    short iz, ir;
    double sum;

    for (iz = 0; iz < nz; iz++) {
        sum = 0.0;
        for (ir = 0; ir < nr; ir++) sum += Out_Ptr->A_rz[ir][iz];
        Out_Ptr->A_z[iz] = sum;
    }

    sum = 0.0;
    for (iz = 0; iz < nz; iz++) {
        sum += Out_Ptr->A_z[iz];
        Out_Ptr->A_l[IzToLayer(iz, In_Parm)] += Out_Ptr->A_z[iz];
    }
    Out_Ptr->A = sum;
}

void Sum2DTt(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    short nr = In_Parm.nr, na = In_Parm.na;
    short ir, ia;
    double sum;

    for (ir = 0; ir < nr; ir++) {
        sum = 0.0;
        for (ia = 0; ia < na; ia++) sum += Out_Ptr->Tt_ra[ir][ia];
        Out_Ptr->Tt_r[ir] = sum;
    }

    for (ia = 0; ia < na; ia++) {
        sum = 0.0;
        for (ir = 0; ir < nr; ir++) sum += Out_Ptr->Tt_ra[ir][ia];
        Out_Ptr->Tt_a[ia] = sum;
    }

    sum = 0.0;
    for (ir = 0; ir < nr; ir++) sum += Out_Ptr->Tt_r[ir];
    Out_Ptr->Tt = sum;
}

/***********************************************************
 * Scaling to physical units.
 ****/
void ScaleRdTt(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    short nr = In_Parm.nr, na = In_Parm.na;
    double dr = In_Parm.dr, da = In_Parm.da;
    short ir, ia;
    double scale1, scale2;

    scale1 = 4.0 * PI * PI * dr * sin(da / 2) * dr * In_Parm.num_photons;

    for (ir = 0; ir < nr; ir++)
        for (ia = 0; ia < na; ia++) {
            scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
            Out_Ptr->Rd_ra[ir][ia] *= scale2;
            Out_Ptr->Tt_ra[ir][ia] *= scale2;
        }

    scale1 = 2.0 * PI * dr * dr * In_Parm.num_photons;
    for (ir = 0; ir < nr; ir++) {
        scale2 = 1.0 / ((ir + 0.5) * scale1);
        Out_Ptr->Rd_r[ir] *= scale2;
        Out_Ptr->Tt_r[ir] *= scale2;
    }

    scale1 = 2.0 * PI * da * In_Parm.num_photons;
    for (ia = 0; ia < na; ia++) {
        scale2 = 1.0 / (sin((ia + 0.5) * da) * scale1);
        Out_Ptr->Rd_a[ia] *= scale2;
        Out_Ptr->Tt_a[ia] *= scale2;
    }

    scale2 = 1.0 / (double)In_Parm.num_photons;
    Out_Ptr->Rd *= scale2;
    Out_Ptr->Tt *= scale2;
}

void ScaleA(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    short nz = In_Parm.nz, nr = In_Parm.nr, nl = In_Parm.num_layers;
    double dz = In_Parm.dz, dr = In_Parm.dr;
    short iz, ir, il;
    double scale1;

    scale1 = 2.0 * PI * dr * dr * dz * In_Parm.num_photons;
    for (iz = 0; iz < nz; iz++)
        for (ir = 0; ir < nr; ir++)
            Out_Ptr->A_rz[ir][iz] /= (ir + 0.5) * scale1;

    scale1 = 1.0 / (dz * In_Parm.num_photons);
    for (iz = 0; iz < nz; iz++)
        Out_Ptr->A_z[iz] *= scale1;

    scale1 = 1.0 / (double)In_Parm.num_photons;
    for (il = 0; il <= nl + 1; il++)
        Out_Ptr->A_l[il] *= scale1;

    Out_Ptr->A *= scale1;
}

/***********************************************************
 * Sum & scale results of current run and write to disk.
 ****/
static void WriteVersion(FILE* file, const char* Version)
{
    fprintf(file, "%s \t# Version number of the file format.\n\n", Version);
    fprintf(file, "####\n# Data categories include: \n");
    fprintf(file, "# InParm, RAT, \n");
    fprintf(file, "# A_l, A_z, Rd_r, Rd_a, Tt_r, Tt_a, \n");
    fprintf(file, "# A_rz, Rd_ra, Tt_ra \n####\n\n");
}

static void WriteInParm(FILE* file, InputStruct In_Parm)
{
    short i;

    fprintf(file, "InParm \t\t\t# Input parameters. cm is used.\n");
    fprintf(file, "%s \tA\t\t# output file name, ASCII.\n", In_Parm.out_fname);
    fprintf(file, "%ld \t\t\t# No. of photons\n", In_Parm.num_photons);
    fprintf(file, "%G\t%G\t\t# dz, dr [cm]\n", In_Parm.dz, In_Parm.dr);
    fprintf(file, "%hd\t%hd\t%hd\t# No. of dz, dr, da.\n\n", In_Parm.nz, In_Parm.nr, In_Parm.na);

    fprintf(file, "%hd\t\t\t\t\t# Number of layers\n", In_Parm.num_layers);
    fprintf(file, "#n\tmua\tmus\tg\td\t# One line for each layer\n");
    fprintf(file, "%G\t\t\t\t\t# n for medium above\n", In_Parm.layerspecs[0].n);
    for (i = 1; i <= In_Parm.num_layers; i++) {
        LayerStruct s = In_Parm.layerspecs[i];
        fprintf(file, "%G\t%G\t%G\t%G\t%G\t# layer %hd\n",
            s.n, s.mua, s.mus, s.g, s.z1 - s.z0, i);
    }
    fprintf(file, "%G\t\t\t\t\t# n for medium below\n\n", In_Parm.layerspecs[i].n);

    fprintf(file, "%G\t%G\t\t%G\t\t%G\t%G\t\t%G\t\t# PD : (Rx, Ry), Rl, (Tx, Ty), Tl (in cm)\n",
        In_Parm.PD_Rx, In_Parm.PD_Ry, In_Parm.PD_Rl,
        In_Parm.PD_Tx, In_Parm.PD_Ty, In_Parm.PD_Tl);

    fprintf(file, "%hd\t\t%G\t%G\t\t%G\t\t# light source : 1.point 2.Gaussian 3.Flat\n",
        In_Parm.lightType, In_Parm.light_x, In_Parm.light_y, In_Parm.light_l);
}

static void WriteRAT(FILE* file, OutStruct Out_Parm)
{
    fprintf(file, "RAT #Reflectance, absorption, transmission. \n");
    fprintf(file, "%-14.6G \t#Specular reflectance [-]\n", Out_Parm.Rsp);
    fprintf(file, "%-14.6G \t#Diffuse reflectance [-]\n", Out_Parm.Rd);
    fprintf(file, "%-14.6G \t#Absorbed fraction [-]\n", Out_Parm.A);
    fprintf(file, "%-14.6G \t#Transmittance [-]\n", Out_Parm.Tt);
    fprintf(file, "\n");
}

static void WriteNbrPhotons(FILE* file, OutStruct Out_Parm)
{
    fprintf(file, "nbr photons photodetecteur R, T et OP \n");
    fprintf(file, "%d \t#NphR\n", Out_Parm.photonsnbrR);
    fprintf(file, "%d \t#NphT\n", Out_Parm.photonsnbrT);
    fprintf(file, "%d \t#NphOP\n", Out_Parm.nbrphotons);
    fprintf(file, "\n");
}

static void WriteA_layer(FILE* file, short Num_Layers, OutStruct Out_Parm)
{
    short i;
    fprintf(file, "A_l #Absorption as a function of layer. [-]\n");
    for (i = 1; i <= Num_Layers; i++)
        fprintf(file, "%12.4G\n", Out_Parm.A_l[i]);
    fprintf(file, "\n");
}

static void WriteRd_ra(FILE* file, short Nr, short Na, OutStruct Out_Parm)
{
    short ir, ia;
    fprintf(file,
        "%s\n%s\n%s\n%s\n%s\n%s\n",
        "# Rd[r][angle]. [1/(cm2sr)].",
        "# Rd[0][0], [0][1],..[0][na-1]",
        "# Rd[1][0], [1][1],..[1][na-1]",
        "# ...",
        "# Rd[nr-1][0], [nr-1][1],..[nr-1][na-1]",
        "Rd_ra");

    for (ir = 0; ir < Nr; ir++)
        for (ia = 0; ia < Na; ia++) {
            fprintf(file, "%12.4E ", Out_Parm.Rd_ra[ir][ia]);
            if ((ir * Na + ia + 1) % 5 == 0) fprintf(file, "\n");
        }
    fprintf(file, "\n");
}

static void WriteRd_r(FILE* file, short Nr, OutStruct Out_Parm)
{
    short ir;
    fprintf(file, "Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm2]\n");
    for (ir = 0; ir < Nr; ir++)
        fprintf(file, "%12.4E\n", Out_Parm.Rd_r[ir]);
    fprintf(file, "\n");
}

static void WriteRd_a(FILE* file, short Na, OutStruct Out_Parm)
{
    short ia;
    fprintf(file, "Rd_a #Rd[0], [1],..Rd[na-1]. [sr-1]\n");
    for (ia = 0; ia < Na; ia++)
        fprintf(file, "%12.4E\n", Out_Parm.Rd_a[ia]);
    fprintf(file, "\n");
}

static void WriteTt_ra(FILE* file, short Nr, short Na, OutStruct Out_Parm)
{
    short ir, ia;
    fprintf(file,
        "%s\n%s\n%s\n%s\n%s\n%s\n",
        "# Tt[r][angle]. [1/(cm2sr)].",
        "# Tt[0][0], [0][1],..[0][na-1]",
        "# Tt[1][0], [1][1],..[1][na-1]",
        "# ...",
        "# Tt[nr-1][0], [nr-1][1],..[nr-1][na-1]",
        "Tt_ra");

    for (ir = 0; ir < Nr; ir++)
        for (ia = 0; ia < Na; ia++) {
            fprintf(file, "%12.4E ", Out_Parm.Tt_ra[ir][ia]);
            if ((ir * Na + ia + 1) % 5 == 0) fprintf(file, "\n");
        }
    fprintf(file, "\n");
}

static void WriteA_rz(FILE* file, short Nr, short Nz, OutStruct Out_Parm)
{
    short iz, ir;
    fprintf(file,
        "%s\n%s\n%s\n%s\n%s\n%s\n",
        "# A[r][z]. [1/cm3]",
        "# A[0][0], [0][1],..[0][nz-1]",
        "# A[1][0], [1][1],..[1][nz-1]",
        "# ...",
        "# A[nr-1][0], [nr-1][1],..[nr-1][nz-1]",
        "A_rz");

    for (ir = 0; ir < Nr; ir++)
        for (iz = 0; iz < Nz; iz++) {
            fprintf(file, "%12.4E ", Out_Parm.A_rz[ir][iz]);
            if ((ir * Nz + iz + 1) % 5 == 0) fprintf(file, "\n");
        }
    fprintf(file, "\n");
}

static void WriteOP(FILE* file, short Nr, short Nz, OutStruct Out_Parm)
{
    short iz, ir;
    fprintf(file,
        "%s\n%s\n%s\n%s\n%s\n%s\n",
        "# OP[r][z]. [1/cm3]",
        "# OP[0][0], [0][1],..[0][nz-1]",
        "# OP[1][0], [1][1],..[1][nz-1]",
        "# ...",
        "# OP[nr-1][0], [nr-1][1],..[nr-1][nz-1]",
        "OP");

    for (ir = 0; ir < Nr; ir++)
        for (iz = 0; iz < Nz; iz++) {
            fprintf(file, "%12.4E ", Out_Parm.OP[ir][iz]);
            if ((ir * Nz + iz + 1) % 5 == 0) fprintf(file, "\n");
        }
    fprintf(file, "\n");
}

static void WriteA_z(FILE* file, short Nz, OutStruct Out_Parm)
{
    short iz;
    fprintf(file, "A_z #A[0], [1],..A[nz-1]. [1/cm]\n");
    for (iz = 0; iz < Nz; iz++)
        fprintf(file, "%12.4E\n", Out_Parm.A_z[iz]);
    fprintf(file, "\n");
}

static void WriteTt_r(FILE* file, short Nr, OutStruct Out_Parm)
{
    short ir;
    fprintf(file, "Tt_r #Tt[0], [1],..Tt[nr-1]. [1/cm2]\n");
    for (ir = 0; ir < Nr; ir++)
        fprintf(file, "%12.4E\n", Out_Parm.Tt_r[ir]);
    fprintf(file, "\n");
}

static void WriteTt_a(FILE* file, short Na, OutStruct Out_Parm)
{
    short ia;
    fprintf(file, "Tt_a #Tt[0], [1],..Tt[na-1]. [sr-1]\n");
    for (ia = 0; ia < Na; ia++)
        fprintf(file, "%12.4E\n", Out_Parm.Tt_a[ia]);
    fprintf(file, "\n");
}

void SumScaleResult(InputStruct In_Parm, OutStruct* Out_Ptr)
{
    Sum2DRd(In_Parm, Out_Ptr);
    Sum2DA(In_Parm, Out_Ptr);
    Sum2DTt(In_Parm, Out_Ptr);

    ScaleRdTt(In_Parm, Out_Ptr);
    ScaleA(In_Parm, Out_Ptr);
}

void WriteResult(InputStruct In_Parm, OutStruct Out_Parm, char* TimeReport)
{
    FILE* file = fopen(In_Parm.out_fname, "w");
    if (file == NULL) nrerror("Cannot open file to write.\n");

    if (toupper(In_Parm.out_fformat) == 'A') WriteVersion(file, "A1");
    else                                      WriteVersion(file, "B1");

    fprintf(file, "# %s", TimeReport);
    fprintf(file, "\n");

    WriteInParm(file, In_Parm);
    WriteRAT(file, Out_Parm);
    WriteNbrPhotons(file, Out_Parm);

    WriteA_layer(file, In_Parm.num_layers, Out_Parm);
    WriteA_z(file, In_Parm.nz, Out_Parm);
    WriteRd_r(file, In_Parm.nr, Out_Parm);
    WriteRd_a(file, In_Parm.na, Out_Parm);
    WriteTt_r(file, In_Parm.nr, Out_Parm);
    WriteTt_a(file, In_Parm.na, Out_Parm);

    WriteA_rz(file, In_Parm.nr, In_Parm.nz, Out_Parm);
    WriteRd_ra(file, In_Parm.nr, In_Parm.na, Out_Parm);
    WriteTt_ra(file, In_Parm.nr, In_Parm.na, Out_Parm);
    WriteOP(file, In_Parm.nr, In_Parm.nz, Out_Parm);

    fclose(file);
}
