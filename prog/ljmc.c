/* Lennard-Jones monte carlo for entropic sampling */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_AV
#define ZCOM_HIST
#define ZCOM_ARGOPT
#include "zcom.h"
#include "avb.h"

#define D 3

int N = 108;
real rho = 0.7f;
real rcdef = 2.5f; /* cut-off distance */
real tp = 1.0f;
real thermdt = 0.3f; /* v-rescaling thermostat step size */
int nsteps = 100000000;
int method = 10;  /* thermostat method */
real hooverQ = 100.f;  /* mass of the Hoover thermostat */
int usesw = 0;
real rshift = 2.0f;
real mcamp = 0.1f;
int betmeth = 0;
real emin = -4.f;
real emax =  2.f;
real edel = 0.1f;
real vmin = -7.f;
real vmax =  0.f;
real vdel = 0.01f;
int nevery = 1000000;
int nreport = 10000000;

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &N,        "number of particles");
  argopt_add(ao, "-T", "%r", &tp,       "temperature");
  argopt_add(ao, "-r", "%r", &rho,      "density");
  argopt_add(ao, "-c", "%r", &rcdef,    "cutoff distance");
  argopt_add(ao, "-1", "%d", &nsteps,   "number of simulation steps");
  argopt_add(ao, "-a", "%r", &mcamp,    "step size for Metropolics algorithm");
  argopt_add(ao, "-q", "%r", &thermdt,  "time step for the mc thermostat");
  argopt_add(ao, "-m", "%d", &method,   "0: mc samp");
  argopt_add(ao, "-M", "%d", &betmeth,  "0: plain average; 1: ratio average");  
  argopt_add(ao, "-Q", "%r", &hooverQ,  "mass of the thermostat");
  argopt_add(ao, "-w", "%b", &usesw,    "use switched potential");
  argopt_add(ao, "-s", "%r", &rshift,   "potential shift distance");
  argopt_add(ao, "--emin", "%r", &emin,     "minimal total energy");
  argopt_add(ao, "--emax", "%r", &emax,     "maximal total energy");
  argopt_add(ao, "--edel", "%r", &edel,     "total energy interval");
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

INLINE int lj_microcan(lj_t *lj, real amp)
{
  int i;
  real xi[3], du = 0.f, dvir = 0.f, epo, epn;

  i = lj_randmv3d(lj, xi, amp);
  epo = lj->epot;
  du = lj_depotlj3d(lj, i, xi, &dvir);
  epn = epo + du;
  if (epn > lj->etot) return 0;
  /* accept by metropolis algorithm  */
  if (du < 0. || metroacc1((lj->dof*.5f - 1) * log((lj->etot - epo)/(lj->etot - epn)), 1.f)) {
  //if (metroacc1(du, (lj->dof*.5 -1)/(lj->etot - epo))) {
    lj_commitlj3d(lj, i, xi, du, dvir);
    return 1;
  }
  return 0;
}

/* entropic monte carlo */
static void domc(lj_t *lj)
{
  int t, tacc = 0, acc = 0;
  double Emin, Emax;
  static av_t ave[1];
  avb_t *avb;
  hist_t *hs;

  lj->dof = 3*lj->n;

  Emin = emin*N;
  Emax = emax*N;
  avb = avb_open(Emin, Emax, edel*N, 1.0/tp, lj->dof);
  printf("%d bins\n", avb->n);
  /* histogram for the potential energy */
  hs = hs_open(1, vmin*N, vmax*N, vdel*N);
  
  lj->etot = .5f*(avb->emax + avb->emin); /* put it in the middle */

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    acc += lj_microcan(lj, mcamp);
    
    lj->ekin = lj->etot - lj->epot;
    if (betmeth == 0) {
      avb_add(avb, lj->etot, lj->ekin, 0.f);
    } else { /* use the biased ratio method */
      avb_addbetrat(avb, lj->etot, .5*lj->dof/lj->ekin, lj->ekin, 0.f);
    }
    
    /* exact MC sampling */
    tacc += avb_mcvrescale(avb, lj->v, lj->n*lj->d,
        thermdt, lj->epot, &lj->ekin, &lj->tkin);
    lj->etot = lj->ekin + lj->epot;

    av_add(ave, lj->etot);
    hs_add1(hs, 0, lj->epot, 1.0, HIST_VERBOSE); /* add to histogram */
    if (t % nevery == 0) printf("t %d, u %g, ek %g, etot %g, bet %g, tacc %g, acc %g\n",
      t, lj->epot/lj->n, lj->ekin/lj->n, lj->etot/lj->n, 1.f/lj->tkin, 1.0*tacc/t, 1.0*acc/t);

    if (t % nreport == 0) {
      avb_write(avb, "avb.dat");
      hs_save(hs, "epot.his", HIST_ADDAHALF);
    }
  }
  avb_write(avb, "avb.dat");
  avb_close(avb);
  hs_close(hs);
  printf("etot %g\n", av_getave(ave)/lj->n);
}

int main(int argc, char **argv)
{
  lj_t *lj;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rcdef);
  if (usesw) lj_initsw(lj, rshift);

  domc(lj);

  lj_close(lj);
  return 0;
}
