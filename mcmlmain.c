/***********************************************************
 *  Main program for MCML (refined for macOS/Clang).
 ***********************************************************/

#define THINKCPROFILER 0
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif

#include "mcml.h"

/* Forward declarations used in main */
FILE* GetFile(char*);
short ReadNumRuns(FILE*);
void ReadParm(FILE*, InputStruct*);
void CheckParm(FILE*, InputStruct*);
void InitOutputData(InputStruct, OutStruct*);
void FreeData(InputStruct, OutStruct*);
double Rspecular(LayerStruct*);
void LaunchPhoton(double Rspecular, LayerStruct* Layerspecs_Ptr, PhotonStruct* Photon_Ptr,
                  short lightType, double light_x, double light_y, double light_l);
void HopDropSpin(InputStruct*, PhotonStruct*, OutStruct*, int*, double*, short*, short*);
void SumScaleResult(InputStruct, OutStruct*);
void WriteResult(InputStruct, OutStruct, char*);

/* 改动：DoOneRun 现在接收一个 summary CSV 的 FILE* 指针 */
static double sOut[60000] = { 0 };
static short  irOut[60000] = { 0 };
static short  izOut[60000] = { 0 };

/***********************************************************
 * Timing helpers.
 ***********************************************************/
time_t PunchTime(char F, char* Msg)
{
#if GNUCC
    (void)F; (void)Msg;
    return 0;
#else
    static clock_t ut0;
    static time_t  rt0;
    double secs;
    char s[STRLEN];

    if (F == 0) {
        ut0 = clock();
        rt0 = time(NULL);
        return 0;
    } else if (F == 1) {
        secs = (clock() - (double)ut0) / (double)CLOCKS_PER_SEC;
        if (secs < 0) secs = 0;
        sprintf(s, "User time: %8.0lf sec = %8.2lf hr. %s\n",
                secs, secs / 3600.0, Msg);
        puts(s);
        strcpy(Msg, s);
        return (difftime(time(NULL), rt0));
    } else if (F == 2) {
        return (difftime(time(NULL), rt0));
    } else {
        return 0;
    }
#endif
}

void PredictDoneTime(long P1, long Pt)
{
    time_t now, done_time;
    struct tm* date;
    char s[80];

    now = time(NULL);
    date = localtime(&now);
    strftime(s, 80, "%H:%M %x", date);
    printf("Start at %s, ", s);

    done_time = now + (time_t)(PunchTime(2, "") / (double)P1 * (Pt - (double)P1));
    date = localtime(&done_time);
    strftime(s, 80, "%H:%M %x", date);
    printf("End at %s\n", s);
}

void ReportResult(InputStruct In_Parm, OutStruct Out_Parm)
{
    char time_report[STRLEN];
    strcpy(time_report, "Simulation time of this run.");
    PunchTime(1, time_report);

    /* 注意：这里对的是 Out_Parm 的副本 */
    SumScaleResult(In_Parm, &Out_Parm);
    WriteResult(In_Parm, Out_Parm, time_report);
}

/***********************************************************
 * Get input file name from argv (if any).
 ***********************************************************/
void GetFnameFromArgv(int argc, char* argv[], char* input_filename)
{
    if (argc >= 2) strcpy(input_filename, argv[1]);
    else input_filename[0] = '\0';
}

/***********************************************************
 * Execute one independent run.
 * 改动：新增参数 summary_fp 用于写入 CSV 汇总。
 ***********************************************************/
static void DoOneRun(short NumRuns, InputStruct* In_Ptr, FILE* summary_fp)
{
    register long i_photon;
    OutStruct out_parm;
    PhotonStruct photon;
    long num_photons = In_Ptr->num_photons, photon_rep = 10;
    int cpt = 0;

    out_parm.nbrphotons  = 0;
    out_parm.photonsnbrR = 0;
    out_parm.photonsnbrT = 0;

    long double OP_totR = 0.0L, MOP_R = 0.0L, OP_totT = 0.0L, MOP_T = 0.0L;
    long double w_totR = 0.0L, w_totT = 0.0L;
    long double Reflectance, Transmittance, Reflectance_SL, Transmittance_SL;

#if THINKCPROFILER
    InitProfile(200, 200);
    cecho2file("prof.rpt", 0, stdout);
#endif

    InitOutputData(*In_Ptr, &out_parm);
    out_parm.Rsp = Rspecular(In_Ptr->layerspecs);
    i_photon = num_photons;
    PunchTime(0, "");

    double epsilon = 0.0001;

    do {
        photon.OP = 0.0;

        if (num_photons - i_photon == photon_rep) {
            printf("%ld photons & %hd runs left, ", i_photon, NumRuns);
            PredictDoneTime(num_photons - i_photon, num_photons);
            photon_rep *= 10;
        }

        LaunchPhoton(out_parm.Rsp, In_Ptr->layerspecs, &photon,
                     In_Ptr->lightType, In_Ptr->light_x, In_Ptr->light_y, In_Ptr->light_l);
        do {
            HopDropSpin(In_Ptr, &photon, &out_parm, &cpt, sOut, irOut, izOut);
        } while (!photon.dead);

        /* 反射 PD（上表面窗口） */
        if (photon.z <= epsilon &&
            photon.x <= (In_Ptr->PD_Rx + In_Ptr->PD_Rl / 2) &&
            photon.x >= (In_Ptr->PD_Rx - In_Ptr->PD_Rl / 2) &&
            photon.y <= (In_Ptr->PD_Ry + In_Ptr->PD_Rl / 2) &&
            photon.y >= (In_Ptr->PD_Ry - In_Ptr->PD_Rl / 2))
        {
            out_parm.photonsnbrR += 1;
            out_parm.nbrphotons  += 1;
            if (cpt > 0) w_totR += sOut[cpt - 1];
            else         w_totR += 1 - out_parm.Rsp;
            OP_totR += photon.OP;

            for (int i = 0; i < cpt; i++) {
                short ir_temp = irOut[i];
                short iz_temp = izOut[i];
                double s_temp = sOut[i];
                out_parm.OP[ir_temp][iz_temp] += s_temp;
                irOut[i] = 0; izOut[i] = 0; sOut[i] = 0.0;
            }
        }

        /* 透射 PD（下表面窗口） */
        if (photon.z >= (In_Ptr->dz * In_Ptr->nz - epsilon) &&
            photon.x <= (In_Ptr->PD_Tx + In_Ptr->PD_Tl / 2) &&
            photon.x >= (In_Ptr->PD_Tx - In_Ptr->PD_Tl / 2) &&
            photon.y <= (In_Ptr->PD_Ty + In_Ptr->PD_Tl / 2) &&
            photon.y >= (In_Ptr->PD_Ty - In_Ptr->PD_Tl / 2))
        {
            out_parm.photonsnbrT += 1;
            out_parm.nbrphotons  += 1;

            if (cpt > 0) w_totT += sOut[cpt - 1];
            else         w_totT += 1 - out_parm.Rsp;

            OP_totT += photon.OP;

            for (int i = 0; i < cpt; i++) {
                short ir_temp = irOut[i];
                short iz_temp = izOut[i];
                double s_temp = sOut[i];
                out_parm.OP[ir_temp][iz_temp] += s_temp;
                irOut[i] = 0; izOut[i] = 0; sOut[i] = 0.0;
            }
        }

        cpt = 0;

    } while (--i_photon);

    /* 这些是 PD 口径的度量，仅打印参考 */
    Reflectance      = (out_parm.photonsnbrR > 0) ? (w_totR / out_parm.photonsnbrR) : 0.0L;
    Transmittance    = (out_parm.photonsnbrT > 0) ? (w_totT / out_parm.photonsnbrT) : 0.0L;
    Reflectance_SL   = (num_photons > 0) ? (w_totR / num_photons) : 0.0L;
    Transmittance_SL = (num_photons > 0) ? (w_totT / num_photons) : 0.0L;
    MOP_R            = (out_parm.photonsnbrR > 0) ? (OP_totR / out_parm.photonsnbrR) : 0.0L;
    MOP_T            = (out_parm.photonsnbrT > 0) ? (OP_totT / out_parm.photonsnbrT) : 0.0L;

    printf("Reflectance_CHatterjee = %1.16Lf;\nTransmittance_CHatterjee = %1.16Lf;\n",
           Reflectance, Transmittance);
    printf("Reflectance_SL = %1.16Lf;\nTransmittance_SL = %1.16Lf;\n",
           Reflectance_SL, Transmittance_SL);
    printf("Reflectance : MOP = %f cm;\n", (double)MOP_R);
    printf("Transmittance : MOP = %f cm.\n", (double)MOP_T);

    /* 写单个 run 的结果文件（包含 RAT） */
    ReportResult(*In_Ptr, out_parm);

    /* 关键：在当前 out_parm（本体）上再做一次 scale，
       这样我们拿到的 out_parm.Rd / out_parm.Tt 与文件里 RAT 对齐 */
    SumScaleResult(*In_Ptr, &out_parm);

    /* 追加到 summary.csv，并在终端打印一行 */
    if (summary_fp) {
        /* CSV: output,Rd,Tt */
        fprintf(summary_fp, "%s,%.10g,%.10g\n", In_Ptr->out_fname, out_parm.Rd, out_parm.Tt);
        fflush(summary_fp);
    }
    printf("[SUMMARY] %s -> Rd=%.6g, Tt=%.6g\n", In_Ptr->out_fname, out_parm.Rd, out_parm.Tt);

    FreeData(*In_Ptr, &out_parm);
}

/***********************************************************
 * main
 * 改动：打开 summary.csv，循环中把文件指针传给 DoOneRun。
 ***********************************************************/
int main(int argc, char* argv[])
{
    char input_filename[STRLEN];
    FILE* input_file_ptr;
    short num_runs;
    InputStruct in_parm;

    ShowVersion("Version LIP6, 2025");
    GetFnameFromArgv(argc, argv, input_filename);
    input_file_ptr = GetFile(input_filename);
    CheckParm(input_file_ptr, &in_parm);
    num_runs = ReadNumRuns(input_file_ptr);

    /* 打开汇总 CSV（覆盖写） */
    const char* summary_name = "summary.csv";
    FILE* summary_fp = fopen(summary_name, "w");
    if (!summary_fp) {
        fprintf(stderr, "Warning: cannot open %s for writing. Summary CSV will be skipped.\n", summary_name);
    } else {
        /* 表头 */
        fprintf(summary_fp, "output,Rd,Tt\n");
        fflush(summary_fp);
    }

    while (num_runs--) {
        ReadParm(input_file_ptr, &in_parm);
        DoOneRun(num_runs, &in_parm, summary_fp);
    }

    if (summary_fp) {
        fclose(summary_fp);
        printf("All runs done. Summary written to %s\n", summary_name);
    }

    fclose(input_file_ptr);
    return 0;
}
