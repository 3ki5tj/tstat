/* Lennard-Jones molecular dynamics volume */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_AV
#define ZCOM_HIST
#define ZCOM_ARGOPT
#include "zcom.h"
#include "avp.h"

#define D 3

int N = 108;
real rho = (real) 0.5;
real rhomin = (real) 0.05;
real rhomax = (real) 0.75;
real rhodel = (real) 0.01;
real rcdef = (real) 1000.0; /* cut-off distance, should be half box */
real tp = (real) 1.2;
real pressure = (real) 0.1; /* default pressure */
real mddt = (real) 0.001; /* time step for molecular dynamics */
real thermdt = (real) 0.01; /* v-rescaling thermostat step size */
real barodt = (real) 0.0001; /* step size for the barostat */
real baroamp = (real) 0.05; /* step size for changing the volume in MC-stat */
int barofreq = 1; /* frequency for Monte-Carlo barostat */
int ensexp = 2; /* ensemble is defined as dV / V^ensexp */
real dlnvmax = (real) 0.1; /* maximal percentage in a Langevin volumn move */
real nhQ = 300; /* mass for the Nose-Hoover thermostat */
real nhW = 1000; /* mass for the Nose-Hoover barostat */
int nsteps = 10000000;
int method = 10;  /* thermostat method */
int nevery = 100000;
int nreport = 2000000;

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &N,        "number of particles");
  argopt_add(ao, "-T", "%r", &tp,       "temperature");
  argopt_add(ao, "-r", "%r", &rho,      "density");
  argopt_add(ao, "--rmin", "%r", &rhomin,   "minimal density");
  argopt_add(ao, "--rmax", "%r", &rhomax,   "maximal density");
  argopt_add(ao, "--rdel", "%r", &rhodel,   "delta density");
  argopt_add(ao, "-:", "%r", &rhodel,   "delta density");
  argopt_add(ao, "-1", "%d", &nsteps,   "number of simulation steps");
  argopt_add(ao, "-d", "%r", &mddt,     "time step for molecular dynamics");
  argopt_add(ao, "-q", "%r", &thermdt,  "time step for vrescaling thermostat");
  argopt_add(ao, "-m", "%d", &method,   "0: mc samp; 1: p-langevin; 2: Nose-Hoover");
  argopt_add(ao, "-a", "%r", &baroamp,  "amplitude for change volume by MC");
  argopt_add(ao, "-f", "%d", &barofreq, "mc change vol. every this # of steps");
  argopt_add(ao, "-b", "%r", &barodt,   "step size for changing volume by Langevin");
  argopt_add(ao, "-x", "%d", &ensexp,   "ensemble exponent");
  argopt_add(ao, "-Q", "%r", &nhQ,      "mass for the Nose-Hoover thermostat");
  argopt_add(ao, "-W", "%r", &nhW,      "mass for the Nose-Hoover barostat");
  argopt_add(ao, "--every", "%d", &nevery,   "print message every this # of steps");
  argopt_add(ao, "--report", "%d", &nreport, "save data every this # of steps");   
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  argopt_dump(ao);
  argopt_close(ao);
}

/* exact entropic Monte Carlo like move
 * r = r*s, p = p/s; */
INLINE int mcprescaleX(avp_t *avp, lj_t *lj, real baroamp, real tp)
{
  int i, acc = 0;
  real vo, vn, lnvo, lnvn, lo = lj->l, s, epo, bet = 1.f/tp;
  double dW, dex;
  lj_t *lj1;

  vo = lj->vol;
  epo = lj->epot;
  lj1 = lj_clone(lj, LJ_CPF); /* make a copy */
  lnvo = (real) log(vo);
  /* compute the internal pressure */
  lnvn = lnvo + baroamp * (2.f * rnd0() - 1.f);
  vn = (real) exp(lnvn);
  if (vn < avp->vmin || vn >= avp->vmax) return 0;
  lj_setrho(lj, lj->n/vn);
  lj_force(lj); /* we change force here */
  dW = avp_getdW(avp, vo, vn);
  dex = bet*(lj->epot - epo + dW)
      + bet*(pow(vo/vn, 2.0/lj->d) - 1)*lj->ekin
      + (lnvo - lnvn) * (1 - ensexp);
  if (metroacc1(dex, 1.f)) { /* scale the velocities */
    s = lo/lj->l;
    for (i = 0; i < lj->d * lj->n; i++) lj->v[i] *= s;
    lj->ekin *= s*s;
    lj->tkin *= s*s;
    acc = 1;
  } else {
    lj_copy(lj, lj1, LJ_CPF); /* restore force etc. */
  }
  lj_close(lj1);
  return acc;
}


/* molecular dynamics simulation */
static void domd(lj_t *lj)
{
  int t, pacc = 0, ptot = 0, tacc = 0;
  static av_t avpr[1], avv[1];
  avp_t *avp;
  real pr, prv, eta = 0, zeta = 0, Vmin, Vmax;
  hist_t *hs;

  avp = avp_open(rhomin, rhomax, rhodel, N, lj->dof, ensexp, pressure, 100.0);
  Vmin = N / avp->rhomax;
  Vmax = N / avp->rhomin;
  printf("%d bins\n", avp->n);

  /* histogram for the potential energy */
  hs = hs_open(1, 0.0, 2.0*N/rhomin, 1.0);

/* add the instaneneous pressure (pr) and get the thermostatic average (prv) */
#define GETPRV() \
  pr = lj_calcpk(lj); \
  avp_add(avp, lj->vol, pr, 0.f); \
  prv = avp_getprv(avp, lj->vol)

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    if (method % 10 == 2) { /* Nose-Hoover */
      lj->ekin = md_ekin(lj->v, lj->n*lj->d, lj->dof, &lj->tkin);
      GETPRV();
      lj_hoovertp(lj, .5f*mddt, tp, prv, &zeta, &eta, nhQ, nhW, ensexp);
    }

    if (method % 10 == 2) { /* Nose-Hoover */
      lj_vv_hoovertp(lj, mddt, eta);
    } else {
      lj_vv(lj, mddt);
    }
    lj_shiftcom(lj, lj->v);
    /* angular momenta are not conserved because of PBC */
    /* lj_shiftang(lj, lj->x, lj->v);  */
    lj->ekin = md_ekin(lj->v, lj->n*lj->d, lj->dof, &lj->tkin);
    
    GETPRV();

    if (method == 0 || method == 10) {
      if (t % barofreq == 0)
      {
        ptot += 1;
        if (method == 0) {
          pacc += lj_mcprescale(lj, baroamp, tp, prv, Vmin, Vmax, ensexp);
        } else {
          pacc += mcprescaleX(avp, lj, baroamp, tp);
        }
      }
    } else if (method == 1) {
      lj_prescale(lj, barodt, tp, prv, Vmin, Vmax, dlnvmax, ensexp);
    //} else if (method == 11) { /* totally wrong */
    //  lj_pberendsen(lj, barodt, tp, prv);
    } else if (method % 10 == 2) { /* Nose-Hoover */
      lj_hoovertp(lj, .5f*mddt, tp, prv, &zeta, &eta, nhQ, nhW, ensexp);
    }

    if (method % 10 != 2) /* unless for Nose-Hoover, T/P-stats are indep. */
      lj_vrescalex(lj, tp, thermdt);

    /* basic statistics, only for monitoring the progress */
    av_add(avpr, pr);
    av_add(avv, lj->vol);
    hs_add1(hs, 0, lj->vol, 1.0, HIST_VERBOSE); /* add to histogram */
    if (t % nevery == 0)
        printf("t %d, u %g, ek %g, bet %g, rho %g(%g), pr %g, prv %g, pacc %g, tacc %g\n",
          t, lj->epot/lj->n, lj->ekin/lj->n, 1.f/lj->tkin, lj->n/lj->vol, 
          lj->n/av_getave(avv), pr, prv, 1.0*pacc/(1e-6 + ptot), 1.0*tacc/t);
    if (t % nreport == 0) {
      avp_write(avp, "avp.dat");
      hs_save(hs, "vol.his", HIST_ADDAHALF);
    }
  }
  avp_close(avp);
  hs_close(hs);
}

int main(int argc, char **argv)
{
  lj_t *lj;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rcdef);

  domd(lj);

  lj_close(lj);
  return 0;
}
