#ifndef AVB_H__
#define AVB_H__

/* routines for sampling a flat energy distribution
  Copyright (c) 2012 Cheng Zhang

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/

typedef struct {
  double emin, emax, edel;
  double bet0; /* default temperature */
  int n, dof;
  av_t *av, *av0, *avu, *avdb, *avdb0;
} avb_t;

INLINE avb_t *avb_open(double emin, double emax, double edel,
    double bet0, int dof)
{
  avb_t *avb;
  int i;

  xnew(avb, 1);
  avb->emin = emin;
  avb->emax = emax;
  avb->n = (int)((emax - emin)/edel + .5);
  avb->bet0 = bet0;
  avb->dof = dof;
  avb->edel = (emax - emin)/avb->n;
  xnew(avb->av,  avb->n + 1);
  xnew(avb->av0, avb->n + 1);
  xnew(avb->avu, avb->n + 1);
  xnew(avb->avdb, avb->n + 1);
  xnew(avb->avdb0, avb->n + 1);
  for (i = 0; i <= avb->n; i++) {
    av_clear( &(avb->av[i]) );
    av_clear( &(avb->av0[i]) );
    av_clear( &(avb->avu[i]) );
    av_clear( &(avb->avdb[i]) );
    av_clear( &(avb->avdb0[i]) );
  }
  return avb;
}

INLINE void avb_close(avb_t *avb)
{
  free(avb->av);
  free(avb->av0);
  free(avb->avu);
  free(avb->avdb);
  free(avb->avdb0);
  free(avb);
}

/* return 0 if successful
   for total-energy entropic sampling */
INLINE int avb_add(avb_t *avb, double etot, double ekin, double ngam)
{
  int i;
  real fhm1;

  if (etot < avb->emin || etot >= avb->emax) return 1;
  i = (int)((etot - avb->emin) / avb->edel);
  die_if(i < 0 || i >= avb->n, "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  fhm1 = .5f*avb->dof - 1;
  av_gadd( &(avb->av[i]), fhm1/ekin, ngam);
  av_add( &(avb->av0[i]), fhm1/ekin);
  av_add( &(avb->avu[i]), etot - ekin);
  av_gadd( &(avb->avdb[i]), fhm1/(ekin*ekin), ngam); 
  av_add( &(avb->avdb0[i]), fhm1/(ekin*ekin)); 
  return 0;
}

/* return 0 if successful
   for potential-energy entropic sampling */
INLINE int avb_addbet(avb_t *avb, double etot, double bet, double f2, double ngam)
{
  int i;
  real ekin;

  if (etot < avb->emin || etot >= avb->emax) return 1;
  i = (int)((etot - avb->emin) / avb->edel);
  die_if(i < 0 || i >= avb->n,
      "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  av_gadd( &(avb->av[i]), bet, ngam);
  av_add( &(avb->av0[i]), bet);
  av_add( &(avb->avu[i]), f2);
  ekin = (.5f*avb->dof - 1)/bet;
  av_gadd( &(avb->avdb[i]), bet/ekin, ngam);
  av_add( &(avb->avdb0[i]), bet/ekin);
  return 0;
}

/* return 0 if successful (ratio method)
   for the formula bet = < Lap U > / < F.F >
   can also be < lap K > / < p.p > */
INLINE int avb_addbetrat(avb_t *avb, double etot, double bet, double f2, double ngam)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax) return 1;
  i = (int)((etot - avb->emin) / avb->edel);
  die_if(i < 0 || i >= avb->n, "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  av_gaddw( &(avb->av[i]), bet, f2, ngam);
  av_addw( &(avb->av0[i]), bet, f2);
  av_add( &(avb->avu[i]), f2);
  /* we need derivative of Omega_N' here, which is done as an approximation */
  /* TODO */
  return 0;
}

#define AVB_MINDATA 100.0 /* for beta */
#define AVB_MINDATA2 1000.0 /* for variance etc */

/* we check avu[i].s, because av0[i].s may be weighted
 * in the ratio estiamtor, which does not give the # of visits
 * avdb0[i].s may not be set properly */
#define AVB_GETBET(avb, i) ( avb->avu[i].s < AVB_MINDATA ? avb->bet0 : av_getave( &(avb->av[i])) )

/* get the local temperature */
INLINE double avb_getbet(avb_t *avb, double etot)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax) return avb->bet0;
  i = (int)( (etot - avb->emin) / avb->edel );
  die_if(i < 0 || i >= avb->n,
      "getbet: index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  return AVB_GETBET(avb, i);
}

/* get the local dU/dE
 * return 0 if successful */
INLINE int avb_getdUdE(avb_t *avb, double etot, double *dUdE)
{
  int i;
  double bet, dbdE;
  
  if (etot < avb->emin || etot >= avb->emax)
    return 1;
  i = (int)( (etot - avb->emin) / avb->edel );
  die_if(i < 0 || i >= avb->n,
      "getdUdE: index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  if (avb->avu[i].s < AVB_MINDATA2) return 2;
  bet = AVB_GETBET(avb, i);
  dbdE = av_getvar( &(avb->av[i]) ) - av_getave( &(avb->avdb[i]) );
  /* the correction must be positive */
  *dUdE = dblmax( 1 + .5*avb->dof*dbdE/(bet*bet), 1e-10 );
  return 0;
}

/* return the correction of bet for a flat potential-energy histogram */
INLINE double avb_getbetcor(avb_t *avb, double etot)
{
  int i1, i2;
  double bet1, bet2, dbde1, dbde2, dude1, dude2, cor;
  
  if (etot < avb->emin || etot >= avb->emax)
    return 0;
  i1 = (int)( (etot - avb->emin) / avb->edel );
  /* NOTE: we differentiate i1 and i1+1, so i1 == n - 1 is unacceptable */
  if (i1 < 0 || i1 >= avb->n - 1) return 0;

  if (avb->avu[i1].s < AVB_MINDATA2) return 0;
  bet1 = AVB_GETBET(avb, i1);
  dbde1 = av_getvar( &(avb->av[i1]) ) - av_getave( &(avb->avdb[i1]) );
  /* the correction must be positive */
  dude1 = dblmax( 1 + .5*avb->dof*dbde1/(bet1*bet1), 1e-10 );

  /* although we can do the thermal averaging
   * we use discrete differentiation, so we need the second bin */
  i2 = i1+1;
  if (avb->avu[i2].s < AVB_MINDATA2) return 0;
  bet2 = AVB_GETBET(avb, i2);
  dbde2 = av_getvar( &(avb->av[i2]) ) - av_getave( &(avb->avdb[i2]) );
  dude2 = dblmax( 1 + .5*avb->dof*dbde2/(bet2*bet2), 1e-10 );

  /* correction */
  if (dude1 > 0 && dude2 > 0)
    cor = log(dude1/dude2);
  else
    cor = 2.*(dude1 - dude2)/(dude1 + dude2);
  cor /= avb->edel;
  /* we confine the maximal correction
   * strictly, shouldn't do this, but more stable this way */
  cor = dblconfine(cor, -bet1*.1, bet1*.1);
  return cor;
}

/* get the entropy difference: dS = S(e2) - S(e1)
 * *dScor is a correction intended to get the flat potential-energy histogram
 * but it makes unweighted averaging approximate */ 
INLINE double avb_getdS(avb_t *avb, double e1, double e2, double *dScor)
{
  int i1, i2, i, sgn = 1, a1, a2;
  double dS = 0, dude1 = 0, dude2 = 0;
 
  if (dScor) *dScor = 0; /* zero the correction by default */

  if (e1 > e2) { sgn = -1; dblswap(e1, e2); }
  if (e2 < avb->emin || e1 >= avb->emax)
    return (e2 - e1) * sgn * avb->bet0;

  /* make no correction for E outside of the range */
  if (e1 < avb->emin) { dS += (avb->emin - e1) * avb->bet0; e1 = avb->emin; }
  if (e2 >= avb->emax) { dS += (e2 - avb->emax) * avb->bet0; e2 = avb->emax - 1e-8; }
  
  i1 = (int)( (e1 - avb->emin) / avb->edel );
  i2 = (int)( (e2 - avb->emin) / avb->edel );
  die_if(i1 < 0 || i1 >= avb->n || i2 < 0 || i2 >= avb->n,
    "i1 %d, i2 %d out of range, n %d, e1 %g, e2 %g\n", i1, i2, avb->n, e1, e2);

  if (i1 == i2) {
    dS += (e2 - e1) * AVB_GETBET(avb, i1);
    /* no correction here */
  } else {
    dS += ((avb->emin + (i1+1)*avb->edel) - e1) * AVB_GETBET(avb, i1);
    dS += (e2 - (avb->emin + i2*avb->edel)) * AVB_GETBET(avb, i2);
    for (i = i1 + 1; i < i2; i++)
      dS += AVB_GETBET(avb, i) * avb->edel;

    /* correction */
    if (dScor != NULL
     && avb_getdUdE(avb, e1, &dude1) == 0
     && avb_getdUdE(avb, e2, &dude2) == 0) {
      if (dude1 > 0 && dude2 > 0)
        *dScor = log(dude1/dude2);
      else
        *dScor = 2.*(dude1 - dude2)/(dude1 + dude2);
      dS += *dScor;
    }
  }
  return sgn * dS;
}



/* write the output */
INLINE int avb_write(avb_t *avb, const char *fn)
{
  int i;
  FILE *fp;
  double bet, cnt, bet0, cnt0, uav, tot, scal, db, db0, devb, devb0;

  xfopen(fp, fn, "w", return -1);
  /* count the total */
  for (tot = 0, i = 0; i < avb->n; i++) {
    tot += avb->avu[i].s;
  }
  scal = 1.0/(tot*avb->edel);
  fprintf(fp, "# %d %g %g %g %d %g\n", avb->n, avb->emin, avb->edel,
      tot, avb->dof, avb->bet0);
  
  for (i = 0; i < avb->n; i++) {
    bet = av_getave( &(avb->av[i]) );
    devb = av_getdev( &(avb->av[i]) );
    db = devb*devb - av_getave( &(avb->avdb[i]) );
    cnt = avb->av[i].s;
    bet0 = av_getave( &(avb->av0[i]) );
    devb0 = av_getdev( &(avb->av0[i]) );
    db0 = devb0*devb0 - av_getave( &(avb->avdb0[i]) );
    cnt0 = avb->avu[i].s;
    uav = av_getave( &(avb->avu[i]) );
    fprintf(fp, "%g %g %g %g %g %g %g %g %g %g %g\n",
      avb->emin + i * avb->edel,
      bet, cnt, bet0, cnt0*scal, uav, cnt0, devb, db, devb0, db0);
  }
  fclose(fp);
  return 0;
}

/* change the kinetic energy by a Monte-Carlo move (exact)
 * to avoid the hard total energy limits set etotmin = etotmax = 0
 * epcor: aim at producing a flat potential-energy histogram */
INLINE int avb_mcvrescale(avb_t *avb, real *v, int nd,
  real amp, real ep, real *ekin, real *tkin,
  real etotmin, real etotmax, int epcor)
{
  int i;
  real ek1 = *ekin, s;
  double logek1, logek2, ek2, r, dS, dSc = 0, etot1, etot2;

  logek1 = log(ek1);
  logek2 = logek1 + amp*(2.f*rnd0() - 1);
  ek2 = exp(logek2);
  etot1 = ek1 + ep;
  etot2 = ek2 + ep;

  /* causes troubles in LJ with the following limitation
   * but might be okay in other cases */
  if (etotmin < etotmax
   && (etot1 > etotmin && etot1 < etotmax)
   && (etot2 > etotmax || etot2 < etotmin))
    return 0;

  dS = avb_getdS(avb, etot1, etot2, epcor ? &dSc : 0);
  r = dS - .5*avb->dof*(logek2 - logek1);
  
  if (r <= 0 || rnd0() < exp(-r)) {
    s = (real) sqrt(ek2/ek1);
    for (i = 0; i < nd; i++)
      v[i] *= s;
    *ekin = ek2;
    if (tkin) *tkin *= s*s;
    return 1;
  } else { /* do nothing otherwise */
    return 0;
  }
}

/* Exact Andersen theromstat
 * epcor: aim at producing a flat potential-energy histogram */
INLINE int avb_mcandersen(avb_t *avb, real *v, int n, int d,
  real ep, real *ekin, real *tkin, int epcor)
{
  int i, j;
  real ek1 = *ekin, ek2, eki1, eki2, vi[3];
  double r, sqtp, dS, dSc = 0, bet1, bet2, etot1, etot2, lnrose;

  i = (int)(rnd0() * n);
  etot1 = ek1 + ep;
  bet1 = avb_getbet(avb, etot1);
  sqtp = 1.0/sqrt(bet1);
  for (eki1 = eki2 = 0, j = 0; j < d; j++) {
    eki1 += .5 * v[i*d + j] * v[i*d + j];
    vi[j] = (real) (sqtp * grand0());
    eki2 += .5 * vi[j] * vi[j];
  }
  ek2 = ek1 + eki2 - eki1;
  etot2 = ek2 + ep;
  
  bet2 = avb_getbet(avb, etot2);
  lnrose = .5*d*log(bet1/bet2);
  dS = avb_getdS(avb, etot1, etot2, epcor ? &dSc : 0);
  r = bet2 * eki1 - bet1 * eki2 + dS + lnrose;

  if (r <= 0 || rnd0() < exp(-r)) {
    for (j = 0; j < d; j++)
      v[i*d + j] = vi[j];
    *ekin = ek2;
    *tkin = 2.*ek2/avb->dof; 
    return 1;
  } else {
    return 0;
  }
}

#endif

