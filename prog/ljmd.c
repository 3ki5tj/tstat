/* Lennard-Jones molecular dynamics */
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
real tpmin = 0.4f;
real mddt = 0.001f; /* time step for molecular dynamics */
real thermdt = 0.01f; /* v-rescaling thermostat step size */
int nsteps = 10000000;
int method = 10;  /* thermostat method */
int betmeth = 0; /* method for computing the temperature beta */
real hooverQ = 100.f;  /* mass of the Hoover thermostat */
int nhM = 5; /* number Nose-Hoover chain variables */
int usesw = 0;
real rshift = 2.0f;
real emin = -4.f;
real emax =  2.f;
real edel = 0.02f;
real vmin = -7.f;
real vmax =  0.f;
real vdel = 0.01f;
int nevery = 100000;
int nreport = 2000000;

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &N,        "number of particles");
  argopt_add(ao, "-T", "%r", &tp,       "temperature");
  argopt_add(ao, "-r", "%r", &rho,      "density");
  argopt_add(ao, "-c", "%r", &rcdef,    "cutoff distance");
  argopt_add(ao, "-1", "%d", &nsteps,   "number of simulation steps");
  argopt_add(ao, "-d", "%r", &mddt,     "time step for molecular dynamics");
  argopt_add(ao, "-q", "%r", &thermdt,  "time step for vrescaling/mc thermostat");
  argopt_add(ao, "-m", "%d", &method,   "0: mc samp; 1: vrescale; 2: hoover; 3: andersen; 4: langevin; see code for more");
  argopt_add(ao, "-M", "%d", &betmeth,  "0: plain average (bet = <(N/2-1)/K>); 1: ratio average (bet = N/2/<K>)");
  argopt_add(ao, "-Y", "%r", &tpmin,    "minimal temperature to control the freq. for andersen T-stat");
  argopt_add(ao, "-Q", "%r", &hooverQ,  "mass of the thermostat");
  argopt_add(ao, "-R", "%d", &nhM,      "number of Nose-Hoover chain variables");
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

/* multicanonical molecular dynamics simulation */
static void domd(lj_t *lj)
{
  int i, t, tacc = 0;
  static av_t ave[1];
  avb_t *avb;
  double tpe = tp, zeta = 0, Emin, Emax;
  real *nhQ, *nhv;
  hist_t *hs;

  xnew(nhQ, nhM);
  xnew(nhv, nhM);
  for (i = 0; i < nhM; i++) {
    nhQ[i] = (i == 0 ? hooverQ : 1.f);
    nhv[i] = (real) ( sqrt(tp/nhQ[i]) * grand0() );
  }

  if (method % 10 == 3 || method % 10 == 4) { /* allow all degrees of freedom */
    lj->dof = 3*lj->n;
  }

  Emin = emin*N;
  Emax = emax*N;
  avb = avb_open(Emin, Emax, edel*N, 1.0/tp, lj->dof);
  printf("%d bins\n", avb->n);
  /* histogram for the potential energy */
  hs = hs_open(1, vmin*N, vmax*N, vdel*N);

#define GETTPE() \
    lj->etot = lj->epots + lj->ekin; \
    if (betmeth == 0) avb_add(avb, lj->etot, lj->ekin, 0.f); \
    else avb_addbetrat(avb, lj->etot, .5*lj->dof/lj->ekin, lj->ekin, 0.f); \
    tpe = 1.0/avb_getbet(avb, lj->etot);

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    if (method % 2 == 0) {
      lj->ekin = md_ekin(lj->v, lj->n*lj->d, lj->dof, &lj->tkin);
      GETTPE();

      if (method == 102) { /* bad Nose-Hoover */
        md_hoover(lj->v, lj->n*lj->d, lj->dof, tpe, .5f*mddt, &zeta, hooverQ, &lj->ekin, &lj->tkin);
      } else if (method <= 22) { /* good (scaled) Nose-Hoover */
        md_hoover(lj->v, lj->n*lj->d, lj->dof, tpe, .5f*mddt, &zeta, hooverQ*tpe/tp, &lj->ekin, &lj->tkin);
      } else if (method == 202) { /* bad Nose-Hoover chain */
        md_nhchain(lj->v, lj->n*lj->d, lj->dof, tpe, 1.f, .5f*mddt, nhv, nhQ, nhM, &lj->ekin, &lj->tkin);
      } else if (method == 32) { /* Nose-Hoover chain */
        md_nhchain(lj->v, lj->n*lj->d, lj->dof, tp, tp/tpe, .5f*mddt, nhv, nhQ, nhM, &lj->ekin, &lj->tkin);
      }

    }
    
    if (method == 22) {
      lj_vv(lj, mddt/tpe);
    } else {
      lj_vv(lj, mddt);
    }
    
    if (method % 10 != 3 && method % 10 != 4) {
      lj_shiftcom(lj, lj->v);
      /* angular momentum can be broken by PBC */
      /* lj_shiftang(lj, lj->x, lj->v); */
      lj->ekin = md_ekin(lj->v, lj->n*lj->d, lj->dof, &lj->tkin);
    }

    GETTPE();
    
    if (method == 0) { /* sampling the kinetic energy */
      tacc += lj_mcvrescale(lj, tpe, thermdt);
    } else if (method == 10){ /* sampling the kinetic energy, exact */
      tacc += avb_mcvrescale(avb, lj->v, lj->n*lj->d,
        thermdt, lj->epots, &lj->ekin, &lj->tkin, 0, 0, 0);

    } else if (method == 1) { /* velocity rescaling (simple) */
      md_vrescale(lj->v, lj->n * lj->d, lj->dof, tpe, thermdt/tpe, &lj->ekin, &lj->tkin);
    } else if (method == 11) { /* velocity rescaling (exact) */
      md_vrescalex(lj->v, lj->n * lj->d, lj->dof, tpe, thermdt/tpe, &lj->ekin, &lj->tkin);
    } else if (method == 101) { /* wrong velocity rescaling */
      md_vrescalex(lj->v, lj->n * lj->d, lj->dof, tpe, thermdt, &lj->ekin, &lj->tkin);


    } else if (method == 2 || method == 22) { /* Nose-Hoover */
      md_hoover(lj->v, lj->n*lj->d, lj->dof, tpe, .5f*mddt, &zeta, hooverQ*tpe/tp, &lj->ekin, &lj->tkin);
    } else if (method == 102) { /* wrong Nose-Hoover */
      md_hoover(lj->v, lj->n*lj->d, lj->dof, tpe, .5f*mddt, &zeta, hooverQ, &lj->ekin, &lj->tkin);
    } else if (method == 32) { /* Nose-Hoover chain */
      md_nhchain(lj->v, lj->n*lj->d, lj->dof, tp, tp/tpe, .5f*mddt, nhv, nhQ, nhM, &lj->ekin, &lj->tkin);
    } else if (method == 202) { /* wrong Nose-Hoover chain */
      md_nhchain(lj->v, lj->n*lj->d, lj->dof, tpe, 1.f, .5f*mddt, nhv, nhQ, nhM, &lj->ekin, &lj->tkin);


    } else if (method % 10 == 3) { /* Andersen */
      if (method == 23) { /* exact Andersen */
        tacc += avb_mcandersen(avb, lj->v, lj->n, lj->d, 
          lj->epots, &lj->ekin, &lj->tkin, 0);
      } else { /* approximate Andersen */ 
        if (rnd0() < tpmin/tpe || method == 13) /* 13: wrong Andersen thermostat */
          md_andersen(lj->v, lj->n, lj->d, tpe);
      }


    } else if (method == 4) { /* Langevin equation */
      md_langevin(lj->v, lj->n, lj->d, tpe, thermdt/tpe);
    } else if (method == 14) { /* wrong langevin dynamics */
      md_langevin(lj->v, lj->n, lj->d, tpe, thermdt);
    }

    av_add(ave, lj->etot);
    hs_add1(hs, 0, lj->epots, 1.0, HIST_VERBOSE); /* add to histogram */
    if (t % nevery == 0) printf("t %d, u %g, ek %g, etot %g, bet %g, tacc %g\n",
      t, lj->epots/lj->n, lj->ekin/lj->n, lj->etot/lj->n, 1.f/lj->tkin, 1.0*tacc/t);

    if (t % nreport == 0) {
      avb_write(avb, "avb.dat");
      hs_save(hs, "epot.his", HIST_ADDAHALF);
    }
  }
  printf("etot %g\n", av_getave(ave)/lj->n);
  avb_write(avb, "avb.dat");
  avb_close(avb);
  hs_close(hs);
  free(nhQ); free(nhv);
}

int main(int argc, char **argv)
{
  lj_t *lj;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rcdef);
  if (usesw) lj_initsw(lj, rshift);

  domd(lj);

  lj_close(lj);
  return 0;
}
