#define ZCOM_PICK
#define ZCOM_CAGO
#define ZCOM_AV
#define ZCOM_HIST
#define ZCOM_ARGOPT
#include "zcom.h"
#include "avb.h"

/* Ca-Go model */
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = 0.5f;
real nbe = 1.f;
real nbc = 4.f;
real rcc = 7.f;

int nsteps = 1000000000;
real mddt = 0.002f;
real tp = 1.5f;
real thermdt = 0.05f;

char *fnpdb = "pdb/1VII.pdb";
real emin = 0.f;
real emax = 140.f;
real edel = 2.0f;
/* parameters for the potential energy histogram */
real vmin = -300.f;
real vmax = 300.f;
real vdel = 0.1f;
int nevery = 100000;
int nreport = 2000000;


/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-f", NULL, &fnpdb,    "PDB file");
  argopt_add(ao, "-T", "%r", &tp,       "default temperature");
  argopt_add(ao, "-1", "%d", &nsteps,   "number of simulation steps");
  argopt_add(ao, "-d", "%r", &mddt,     "time step for molecular dynamics");
  argopt_add(ao, "-q", "%r", &thermdt,  "time step for mc-vrescaling thermostat");
  argopt_add(ao, "-c", "%r", &rcc,      "cutoff distance for selecting contacts");
  argopt_add(ao, "--emin", "%r", &emin,     "minimal total energy");
  argopt_add(ao, "--emax", "%r", &emax,     "maximal total energy");
  argopt_add(ao, "--edel", "%r", &edel,     "total energy interval");
  argopt_add(ao, "-:",     "%r", &edel,     "total energy interval");
  argopt_add(ao, "--vmin", "%r", &vmin,     "minimal potential energy");
  argopt_add(ao, "--vmax", "%r", &vmax,     "maximal potential energy");
  argopt_add(ao, "--vdel", "%r", &vdel,     "potential energy interval");
  argopt_add(ao, "--every", "%d", &nevery,   "print message every this # of steps");
  argopt_add(ao, "--report", "%d", &nreport, "save data every this # of steps");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  argopt_dump(ao);
  argopt_close(ao);
}


static void domd(void)
{
  cago_t *w;
  int it, tacc = 0;
  static av_t avep[1], avet[1];
  real tpe;
  avb_t *avb;
  hist_t *hs;

  w = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc);
  cago_initmd(w, 0.1, 0.0);
  printf("epot %g, epotref %g\n", w->epot, w->epotref);

  avb = avb_open(emin, emax, edel, 1.0/tp, w->dof);
  hs = hs_open(1, vmin, vmax, vdel);

/* // 2CI2
  avb = avb_open(-50, 100, 0.5, w->dof);
  hs = hs_open(1, -160.0, 320.0, 0.01);
*/

  for (it = 1; it <= nsteps; it++) {
    cago_vv(w, 1.0f, mddt);
    cago_rmcom(w, w->x, w->v);
    w->etot = w->epot + w->ekin;
    //avb_add(avb, w->etot, w->ekin, 0.0);
    cago_rotfit(w, w->x, NULL); /* compute rmsd */
    avb_addbet(avb, w->etot, (.5*avb->dof - 1)/w->ekin, w->rmsd, 0.0);
    hs_add1(hs, 0, w->epot, 1.0, HIST_VERBOSE); /* add to histogram */

    tpe = 1.0/avb_getbet(avb, w->etot);
    
    //tacc += avb_mcvrescale(avb, (real *) w->v, w->n * 3, w->dof, 
    //    thermdt, w->epot, &w->ekin, &w->tkin);    
    md_vrescale3d(w->v, w->n, w->dof, tpe, thermdt, &w->ekin, &w->tkin);
    av_add(avep, w->epot);
    av_add(avet, w->etot);

    if (it % nevery == 0) {
      printf("t %d, ep %g/%g, etot %g/%g, T %g/%g, rmsd %g, tacc %g\n", it, w->epot,
        av_getave(avep), w->etot, av_getave(avet), w->tkin, tpe, w->rmsd, 1.0*tacc/it);
    }
    if (it % nreport == 0) {
      avb_write(avb, "goavb.dat");
      cago_writepos(w, w->x, w->v, "go.pos");
      hs_save(hs, "epot.his", HIST_ADDAHALF);
    }
  }
  cago_close(w);
}

int main(int argc, char **argv)
{
  doargs(argc, argv);
  domd();
  return 0;
}
