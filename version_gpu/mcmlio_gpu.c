/***********************************************************
 *  Input/output of data for GPU version (Flattened arrays).
 ***********************************************************/

#include "mcml_gpu.h"

/***********************************************************
 *  Report error and exit.
 ***********************************************************/
void nrerror(char *error_text) {
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

/***********************************************************
 * Name list for duplicate file checks.
 ***********************************************************/
struct NameList {
  char name[STRLEN];
  struct NameList *next;
};
typedef struct NameList NameNode;
typedef NameNode *NameLink;

static char *CenterStr(short Wid, char *InStr, char *OutStr) {
  size_t nspaces = (Wid - (short)strlen(InStr)) / 2;
  if ((short)nspaces < 0)
    nspaces = 0;
  OutStr[0] = '\0';
  while (nspaces--)
    strcat(OutStr, " ");
  strcat(OutStr, InStr);
  return OutStr;
}

#define COLWIDTH 80
void ShowVersion(const char *version) {
  char str[STRLEN];
  CenterStr(COLWIDTH,
            "MCML - Monte Carlo Simulation of Multi-layered Turbid Media (GPU "
            "Version)",
            str);
  puts(str);
  puts("");
  CenterStr(COLWIDTH, "GPU Porting", str);
  puts(str);
  puts("");
  CenterStr(COLWIDTH, (char *)version, str);
  puts(str);
  puts("\n\n\n\n");
}
#undef COLWIDTH

FILE *GetFile(char *Fname) {
  FILE *file = NULL;
  static bool firsttime = true;

  do {
    if (firsttime && Fname[0] != '\0') {
      firsttime = false;
    } else {
      printf("Input filename(or . to exit):");
      if (scanf("%255s", Fname) != 1)
        exit(1);
      firsttime = false;
    }
    if (strlen(Fname) == 1 && Fname[0] == '.')
      exit(1);
    file = fopen(Fname, "r");
  } while (file == NULL);
  return file;
}

static void KillChar(size_t i, char *Str) {
  size_t sl = strlen(Str);
  for (; i < sl; i++)
    Str[i] = Str[i + 1];
}

bool CheckChar(char *Str) {
  bool found = false;
  size_t sl = strlen(Str);
  size_t i = 0;
  while (i < sl) {
    unsigned char c = (unsigned char)Str[i];
    if (c > 127)
      nrerror("Non-ASCII file\n");
    else if (isprint(c) || isspace(c))
      i++;
    else {
      found = true;
      KillChar(i, Str);
      sl--;
    }
  }
  return found;
}

static bool CommentLine(char *Buf) {
  size_t spn = strspn(Buf, " \t");
  size_t cspn = strcspn(Buf, "#\n");
  return (spn == cspn) ? 1 : 0;
}

char *FindDataLine(FILE *File_Ptr) {
  static char buf[STRLEN];
  buf[0] = '\0';
  do {
    if (fgets(buf, 255, File_Ptr) == NULL) {
      printf("Incomplete data\n");
      buf[0] = '\0';
      break;
    } else {
      CheckChar(buf);
    }
  } while (CommentLine(buf));
  return buf;
}

short ReadNumRuns(FILE *File_Ptr) {
  char buf[STRLEN];
  short n = 0;
  FindDataLine(File_Ptr);
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading number of runs\n");
  sscanf(buf, "%hd", &n);
  return n;
}

void ReadFnameFormat(FILE *File_Ptr, InputStruct *In_Ptr) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading file name and format.\n");
  sscanf(buf, "%255s %c", In_Ptr->out_fname, &(In_Ptr->out_fformat));
  if (toupper(In_Ptr->out_fformat) != 'B')
    In_Ptr->out_fformat = 'A';
}

void ReadNumPhotons(FILE *File_Ptr, InputStruct *In_Ptr) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading number of photons.\n");
  sscanf(buf, "%ld", &In_Ptr->num_photons);
}

void ReadDzDr(FILE *File_Ptr, InputStruct *In_Ptr) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading dz, dr.\n");
  sscanf(buf, "%lf%lf", &In_Ptr->dz, &In_Ptr->dr);
}

void ReadNzNrNa(FILE *File_Ptr, InputStruct *In_Ptr) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading info.\n");
  sscanf(buf, "%hd%hd%hd", &In_Ptr->nz, &In_Ptr->nr, &In_Ptr->na);
  In_Ptr->da = 0.5 * PI / In_Ptr->na;
}

void ReadNumLayers(FILE *File_Ptr, InputStruct *In_Ptr) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading number of layers.\n");
  sscanf(buf, "%hd", &In_Ptr->num_layers);
}

void ReadPD(FILE *File_Ptr, InputStruct *In_Ptr) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading PD.\n");
  sscanf(buf, "%lf%lf%lf%lf%lf%lf", &In_Ptr->PD_Rx, &In_Ptr->PD_Ry,
         &In_Ptr->PD_Rl, &In_Ptr->PD_Tx, &In_Ptr->PD_Ty, &In_Ptr->PD_Tl);
}

void ReadLight(FILE *File_Ptr, InputStruct *In_Ptr) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading light.\n");
  sscanf(buf, "%hd%lf%lf%lf", &In_Ptr->lightType, &In_Ptr->light_x,
         &In_Ptr->light_y, &In_Ptr->light_l);
}

void ReadAmbient(FILE *File_Ptr, LayerStruct *Layer_Ptr, char *side) {
  char buf[STRLEN];
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    nrerror("Reading ambient n.\n");
  sscanf(buf, "%lf", &Layer_Ptr->n);
}

bool ReadOneLayer(FILE *File_Ptr, LayerStruct *Layer_Ptr, double *Z_Ptr) {
  char buf[STRLEN];
  double d, n, mua, mus, g;
  strcpy(buf, FindDataLine(File_Ptr));
  if (buf[0] == '\0')
    return 1;
  sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &d);
  Layer_Ptr->n = n;
  Layer_Ptr->mua = mua;
  Layer_Ptr->mus = mus;
  Layer_Ptr->g = g;
  Layer_Ptr->z0 = *Z_Ptr;
  *Z_Ptr += d;
  Layer_Ptr->z1 = *Z_Ptr;
  return 0;
}

void ReadLayerSpecs(FILE *File_Ptr, short Num_Layers,
                    LayerStruct **Layerspecs_PP) {
  char msg[STRLEN];
  short i = 0;
  double z = 0.0;
  *Layerspecs_PP =
      (LayerStruct *)malloc((unsigned)(Num_Layers + 2) * sizeof(LayerStruct));
  if (!(*Layerspecs_PP))
    nrerror("allocation failure");
  ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "top");
  for (i = 1; i <= Num_Layers; i++)
    if (ReadOneLayer(File_Ptr, &((*Layerspecs_PP)[i]), &z)) {
      sprintf(msg, "Error reading %hd of %hd layers\n", i, Num_Layers);
      nrerror(msg);
    }
  ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "bottom");
}

void CriticalAngle(short Num_Layers, LayerStruct **Layerspecs_PP) {
  short i = 0;
  double n1, n2;
  for (i = 1; i <= Num_Layers; i++) {
    n1 = (*Layerspecs_PP)[i].n;
    n2 = (*Layerspecs_PP)[i - 1].n;
    (*Layerspecs_PP)[i].cos_crit0 =
        n1 > n2 ? sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;
    n2 = (*Layerspecs_PP)[i + 1].n;
    (*Layerspecs_PP)[i].cos_crit1 =
        n1 > n2 ? sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;
  }
}

void ReadParm(FILE *File_Ptr, InputStruct *In_Ptr) {
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

static bool NameInList(char *Name, NameLink List) {
  while (List != NULL) {
    if (strcmp(Name, List->name) == 0)
      return true;
    List = List->next;
  }
  return false;
}

static void AddNameToList(char *Name, NameLink *List_Ptr) {
  NameLink list = *List_Ptr;
  if (list == NULL) {
    *List_Ptr = list = (NameLink)malloc(sizeof(NameNode));
    strcpy(list->name, Name);
    list->next = NULL;
  } else {
    while (list->next != NULL)
      list = list->next;
    list->next = (NameLink)malloc(sizeof(NameNode));
    list = list->next;
    strcpy(list->name, Name);
    list->next = NULL;
  }
}

static bool FnameTaken(char *fname, NameLink *List_Ptr) {
  if (NameInList(fname, *List_Ptr))
    return true;
  AddNameToList(fname, List_Ptr);
  return false;
}

static void FreeFnameList(NameLink List) {
  NameLink next;
  while (List != NULL) {
    next = List->next;
    free(List);
    List = next;
  }
}

void CheckParm(FILE *File_Ptr, InputStruct *In_Ptr) {
  short i_run, num_runs;
  NameLink head = NULL;
  bool name_taken;
  char msg[STRLEN];
  num_runs = ReadNumRuns(File_Ptr);
  for (i_run = 1; i_run <= num_runs; i_run++) {
    printf("Checking input data for run %hd\n", i_run);
    ReadParm(File_Ptr, In_Ptr);
    name_taken = FnameTaken(In_Ptr->out_fname, &head);
    if (name_taken)
      sprintf(msg, "file name %s duplicated.\n", In_Ptr->out_fname);
    free(In_Ptr->layerspecs);
    if (name_taken)
      nrerror(msg);
  }
  FreeFnameList(head);
  rewind(File_Ptr);
}

/***********************************************************
 * Init output arrays (Flattened).
 ***********************************************************/
void InitOutputData(InputStruct In_Parm, OutStruct *Out_Ptr) {
  short nz = In_Parm.nz, nr = In_Parm.nr, na = In_Parm.na,
        nl = In_Parm.num_layers;

  if (nz <= 0 || nr <= 0 || na <= 0 || nl <= 0)
    nrerror("Wrong grid parameters.\n");

  Out_Ptr->Rsp = 0.0;
  Out_Ptr->Rd = 0.0;
  Out_Ptr->A = 0.0;
  Out_Ptr->Tt = 0.0;

  Out_Ptr->Rd_ra = (double *)calloc(nr * na, sizeof(double));
  Out_Ptr->Rd_r = (double *)calloc(nr, sizeof(double));
  Out_Ptr->Rd_a = (double *)calloc(na, sizeof(double));

  Out_Ptr->A_rz = (double *)calloc(nr * nz, sizeof(double));
  Out_Ptr->A_z = (double *)calloc(nz, sizeof(double));
  Out_Ptr->A_l = (double *)calloc(nl + 2, sizeof(double));

  Out_Ptr->Tt_ra = (double *)calloc(nr * na, sizeof(double));
  Out_Ptr->Tt_r = (double *)calloc(nr, sizeof(double));
  Out_Ptr->Tt_a = (double *)calloc(na, sizeof(double));

  Out_Ptr->OP = (double *)calloc(nr * nz, sizeof(double));

  Out_Ptr->photonsnbrR = 0;
  Out_Ptr->photonsnbrT = 0;
  Out_Ptr->nbrphotons = 0;
}

void FreeData(InputStruct In_Parm, OutStruct *Out_Ptr) {
  free(In_Parm.layerspecs);
  free(Out_Ptr->Rd_ra);
  free(Out_Ptr->Rd_r);
  free(Out_Ptr->Rd_a);
  free(Out_Ptr->A_rz);
  free(Out_Ptr->A_z);
  free(Out_Ptr->A_l);
  free(Out_Ptr->Tt_ra);
  free(Out_Ptr->Tt_r);
  free(Out_Ptr->Tt_a);
  free(Out_Ptr->OP);
}

/***********************************************************
 * Summations (Looping over flat arrays).
 ***********************************************************/
void Sum2DRd(InputStruct In_Parm, OutStruct *Out_Ptr) {
  short nr = In_Parm.nr, na = In_Parm.na;
  short ir, ia;
  double sum;

  for (ir = 0; ir < nr; ir++) {
    sum = 0.0;
    for (ia = 0; ia < na; ia++)
      sum += Out_Ptr->Rd_ra[ir * na + ia];
    Out_Ptr->Rd_r[ir] = sum;
  }

  for (ia = 0; ia < na; ia++) {
    sum = 0.0;
    for (ir = 0; ir < nr; ir++)
      sum += Out_Ptr->Rd_ra[ir * na + ia];
    Out_Ptr->Rd_a[ia] = sum;
  }

  sum = 0.0;
  for (ir = 0; ir < nr; ir++)
    sum += Out_Ptr->Rd_r[ir];
  Out_Ptr->Rd = sum;
}

static short IzToLayer(short Iz, InputStruct In_Parm) {
  short i = 1;
  short num_layers = In_Parm.num_layers;
  double dz = In_Parm.dz;
  while ((Iz + 0.5) * dz >= In_Parm.layerspecs[i].z1 && i < num_layers)
    i++;
  return i;
}

void Sum2DA(InputStruct In_Parm, OutStruct *Out_Ptr) {
  short nz = In_Parm.nz, nr = In_Parm.nr;
  short iz, ir;
  double sum;

  for (iz = 0; iz < nz; iz++) {
    sum = 0.0;
    for (ir = 0; ir < nr; ir++)
      sum += Out_Ptr->A_rz[ir * nz + iz];
    Out_Ptr->A_z[iz] = sum;
  }

  sum = 0.0;
  for (iz = 0; iz < nz; iz++) {
    sum += Out_Ptr->A_z[iz];
    Out_Ptr->A_l[IzToLayer(iz, In_Parm)] += Out_Ptr->A_z[iz];
  }
  Out_Ptr->A = sum;
}

void Sum2DTt(InputStruct In_Parm, OutStruct *Out_Ptr) {
  short nr = In_Parm.nr, na = In_Parm.na;
  short ir, ia;
  double sum;

  for (ir = 0; ir < nr; ir++) {
    sum = 0.0;
    for (ia = 0; ia < na; ia++)
      sum += Out_Ptr->Tt_ra[ir * na + ia];
    Out_Ptr->Tt_r[ir] = sum;
  }

  for (ia = 0; ia < na; ia++) {
    sum = 0.0;
    for (ir = 0; ir < nr; ir++)
      sum += Out_Ptr->Tt_ra[ir * na + ia];
    Out_Ptr->Tt_a[ia] = sum;
  }

  sum = 0.0;
  for (ir = 0; ir < nr; ir++)
    sum += Out_Ptr->Tt_r[ir];
  Out_Ptr->Tt = sum;
}

/***********************************************************
 * Scaling.
 ***********************************************************/
void ScaleRdTt(InputStruct In_Parm, OutStruct *Out_Ptr) {
  short nr = In_Parm.nr, na = In_Parm.na;
  double dr = In_Parm.dr, da = In_Parm.da;
  short ir, ia;
  double scale1, scale2;

  scale1 = 4.0 * PI * PI * dr * sin(da / 2) * dr * In_Parm.num_photons;

  for (ir = 0; ir < nr; ir++)
    for (ia = 0; ia < na; ia++) {
      scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
      Out_Ptr->Rd_ra[ir * na + ia] *= scale2;
      Out_Ptr->Tt_ra[ir * na + ia] *= scale2;
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

void ScaleA(InputStruct In_Parm, OutStruct *Out_Ptr) {
  short nz = In_Parm.nz, nr = In_Parm.nr, nl = In_Parm.num_layers;
  double dz = In_Parm.dz, dr = In_Parm.dr;
  short iz, ir, il;
  double scale1;

  scale1 = 2.0 * PI * dr * dr * dz * In_Parm.num_photons;
  for (iz = 0; iz < nz; iz++)
    for (ir = 0; ir < nr; ir++)
      Out_Ptr->A_rz[ir * nz + iz] /= (ir + 0.5) * scale1;

  scale1 = 1.0 / (dz * In_Parm.num_photons);
  for (iz = 0; iz < nz; iz++)
    Out_Ptr->A_z[iz] *= scale1;

  scale1 = 1.0 / (double)In_Parm.num_photons;
  for (il = 0; il <= nl + 1; il++)
    Out_Ptr->A_l[il] *= scale1;

  Out_Ptr->A *= scale1;
}

/***********************************************************
 * Write to disk.
 ***********************************************************/
static void WriteVersion(FILE *file, const char *Version) {
  fprintf(file, "%s \t# Version number of the file format.\n\n", Version);
  fprintf(file, "####\n# Data categories include: \n");
  fprintf(file, "# InParm, RAT, \n");
  fprintf(file, "# A_l, A_z, Rd_r, Rd_a, Tt_r, Tt_a, \n");
  fprintf(file, "# A_rz, Rd_ra, Tt_ra \n####\n\n");
}

static void WriteInParm(FILE *file, InputStruct In_Parm) {
  short i;
  fprintf(file, "InParm \t\t\t# Input parameters. cm is used.\n");
  fprintf(file, "%s \tA\t\t# output file name, ASCII.\n", In_Parm.out_fname);
  fprintf(file, "%ld \t\t\t# No. of photons\n", In_Parm.num_photons);
  fprintf(file, "%G\t%G\t\t# dz, dr [cm]\n", In_Parm.dz, In_Parm.dr);
  fprintf(file, "%hd\t%hd\t%hd\t# No. of dz, dr, da.\n\n", In_Parm.nz,
          In_Parm.nr, In_Parm.na);
  fprintf(file, "%hd\t\t\t\t\t# Number of layers\n", In_Parm.num_layers);
  fprintf(file, "#n\tmua\tmus\tg\td\t# One line for each layer\n");
  fprintf(file, "%G\t\t\t\t\t# n for medium above\n", In_Parm.layerspecs[0].n);
  for (i = 1; i <= In_Parm.num_layers; i++) {
    LayerStruct s = In_Parm.layerspecs[i];
    fprintf(file, "%G\t%G\t%G\t%G\t%G\t# layer %hd\n", s.n, s.mua, s.mus, s.g,
            s.z1 - s.z0, i);
  }
  fprintf(file, "%G\t\t\t\t\t# n for medium below\n\n",
          In_Parm.layerspecs[i].n);
  fprintf(file,
          "%G\t%G\t\t%G\t\t%G\t%G\t\t%G\t\t# PD : (Rx, Ry), Rl, (Tx, Ty), Tl "
          "(in cm)\n",
          In_Parm.PD_Rx, In_Parm.PD_Ry, In_Parm.PD_Rl, In_Parm.PD_Tx,
          In_Parm.PD_Ty, In_Parm.PD_Tl);
  fprintf(file,
          "%hd\t\t%G\t%G\t\t%G\t\t# light source : 1.point 2.Gaussian 3.Flat\n",
          In_Parm.lightType, In_Parm.light_x, In_Parm.light_y, In_Parm.light_l);
}

static void WriteRAT(FILE *file, OutStruct Out_Parm) {
  fprintf(file, "RAT #Reflectance, absorption, transmission. \n");
  fprintf(file, "%-14.6G \t#Specular reflectance [-]\n", Out_Parm.Rsp);
  fprintf(file, "%-14.6G \t#Diffuse reflectance [-]\n", Out_Parm.Rd);
  fprintf(file, "%-14.6G \t#Absorbed fraction [-]\n", Out_Parm.A);
  fprintf(file, "%-14.6G \t#Transmittance [-]\n", Out_Parm.Tt);
  fprintf(file, "\n");
}

static void WriteNbrPhotons(FILE *file, OutStruct Out_Parm) {
  fprintf(file, "nbr photons photodetecteur R, T et OP \n");
  fprintf(file, "%d \t#NphR\n", Out_Parm.photonsnbrR);
  fprintf(file, "%d \t#NphT\n", Out_Parm.photonsnbrT);
  fprintf(file, "%d \t#NphOP\n", Out_Parm.nbrphotons);
  fprintf(file, "\n");
}

static void WriteA_layer(FILE *file, short Num_Layers, OutStruct Out_Parm) {
  short i;
  fprintf(file, "A_l #Absorption as a function of layer. [-]\n");
  for (i = 1; i <= Num_Layers; i++)
    fprintf(file, "%12.4G\n", Out_Parm.A_l[i]);
  fprintf(file, "\n");
}

static void WriteRd_ra(FILE *file, short Nr, short Na, OutStruct Out_Parm) {
  short ir, ia;
  fprintf(file, "Rd_ra\n");
  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++) {
      fprintf(file, "%12.4E ", Out_Parm.Rd_ra[ir * Na + ia]);
      if ((ir * Na + ia + 1) % 5 == 0)
        fprintf(file, "\n");
    }
  fprintf(file, "\n");
}

static void WriteRd_r(FILE *file, short Nr, OutStruct Out_Parm) {
  short ir;
  fprintf(file, "Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm2]\n");
  for (ir = 0; ir < Nr; ir++)
    fprintf(file, "%12.4E\n", Out_Parm.Rd_r[ir]);
  fprintf(file, "\n");
}

static void WriteRd_a(FILE *file, short Na, OutStruct Out_Parm) {
  short ia;
  fprintf(file, "Rd_a #Rd[0], [1],..Rd[na-1]. [sr-1]\n");
  for (ia = 0; ia < Na; ia++)
    fprintf(file, "%12.4E\n", Out_Parm.Rd_a[ia]);
  fprintf(file, "\n");
}

static void WriteTt_ra(FILE *file, short Nr, short Na, OutStruct Out_Parm) {
  short ir, ia;
  fprintf(file, "Tt_ra\n");
  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++) {
      fprintf(file, "%12.4E ", Out_Parm.Tt_ra[ir * Na + ia]);
      if ((ir * Na + ia + 1) % 5 == 0)
        fprintf(file, "\n");
    }
  fprintf(file, "\n");
}

static void WriteTt_r(FILE *file, short Nr, OutStruct Out_Parm) {
  short ir;
  fprintf(file, "Tt_r #Tt[0], [1],..Tt[nr-1]. [1/cm2]\n");
  for (ir = 0; ir < Nr; ir++)
    fprintf(file, "%12.4E\n", Out_Parm.Tt_r[ir]);
  fprintf(file, "\n");
}

static void WriteTt_a(FILE *file, short Na, OutStruct Out_Parm) {
  short ia;
  fprintf(file, "Tt_a #Tt[0], [1],..Tt[na-1]. [sr-1]\n");
  for (ia = 0; ia < Na; ia++)
    fprintf(file, "%12.4E\n", Out_Parm.Tt_a[ia]);
  fprintf(file, "\n");
}

static void WriteA_rz(FILE *file, short Nr, short Nz, OutStruct Out_Parm) {
  short iz, ir;
  fprintf(file, "A_rz\n");
  for (ir = 0; ir < Nr; ir++)
    for (iz = 0; iz < Nz; iz++) {
      fprintf(file, "%12.4E ", Out_Parm.A_rz[ir * Nz + iz]);
      if ((ir * Nz + iz + 1) % 5 == 0)
        fprintf(file, "\n");
    }
  fprintf(file, "\n");
}

static void WriteA_z(FILE *file, short Nz, OutStruct Out_Parm) {
  short iz;
  fprintf(file, "A_z #A[0], [1],..A[nz-1]. [1/cm]\n");
  for (iz = 0; iz < Nz; iz++)
    fprintf(file, "%12.4E\n", Out_Parm.A_z[iz]);
  fprintf(file, "\n");
}

static void WriteOP(FILE *file, short Nr, short Nz, OutStruct Out_Parm) {
  short iz, ir;
  fprintf(file, "OP\n");
  for (ir = 0; ir < Nr; ir++)
    for (iz = 0; iz < Nz; iz++) {
      fprintf(file, "%12.4E ", Out_Parm.OP[ir * Nz + iz]);
      if ((ir * Nz + iz + 1) % 5 == 0)
        fprintf(file, "\n");
    }
  fprintf(file, "\n");
}

void SumScaleResult(InputStruct In_Parm, OutStruct *Out_Ptr) {
  Sum2DRd(In_Parm, Out_Ptr);
  Sum2DA(In_Parm, Out_Ptr);
  Sum2DTt(In_Parm, Out_Ptr);
  ScaleRdTt(In_Parm, Out_Ptr);
  ScaleA(In_Parm, Out_Ptr);
}

void WriteResult(InputStruct In_Parm, OutStruct Out_Parm, char *TimeReport) {
  FILE *file = fopen(In_Parm.out_fname, "w");
  if (file == NULL)
    nrerror("Cannot open file to write.\n");
  if (toupper(In_Parm.out_fformat) == 'A')
    WriteVersion(file, "A1");
  else
    WriteVersion(file, "B1");

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
