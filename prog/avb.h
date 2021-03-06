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
  av_t *av, *av0, *avu, *avek;
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
  xnew(avb->avek, avb->n + 1);
  for (i = 0; i <= avb->n; i++) {
    av_clear( &(avb->av[i]) );
    av_clear( &(avb->av0[i]) );
    av_clear( &(avb->avu[i]) );
    av_clear( &(avb->avek[i]) );
  }
  return avb;
}

INLINE void avb_close(avb_t *avb)
{
  free(avb->av);
  free(avb->av0);
  free(avb->avu);
  free(avb->avek);
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
  av_gadd( &(avb->avek[i]), ekin, ngam);
  av_add( &(avb->av0[i]), fhm1/ekin);
  av_add( &(avb->avu[i]), etot - ekin);
  return 0;
}

/* return 0 if successful
   for potential-energy entropic sampling */
INLINE int avb_addbet(avb_t *avb, double etot, double bet, double f2, double ngam)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax) return 1;
  i = (int)((etot - avb->emin) / avb->edel);
  die_if(i < 0 || i >= avb->n,
      "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  av_gadd( &(avb->av[i]), bet, ngam);
  av_gadd( &(avb->avek[i]), (0.5*avb->dof - 1)/bet, ngam);
  av_add( &(avb->av0[i]), bet);
  av_add( &(avb->avu[i]), f2);
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
 * in the ratio estiamtor, which does not give the # of visits */
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

INLINE double avb_getdude_low(const av_t *avbet, const av_t *avek, real dof, double *pbet)
{
  double bet, k;

  bet = av_getave( avbet );
  if (pbet) *pbet = bet;
  k = av_getave( avek );
  return dblmax(k * bet - (.5*dof - 1), 1e-10);
}

/* get the local dU/dE
 * return 0 if successful */
INLINE int avb_getdude(avb_t *avb, double etot, double *dude)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax)
    return 1;
  i = (int)( (etot - avb->emin) / avb->edel );
  die_if(i < 0 || i >= avb->n,
      "getdude: index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  if (avb->avu[i].s < AVB_MINDATA2) return 2;
  *dude = avb_getdude_low( &(avb->av[i]), &(avb->avek[i]), avb->dof, NULL );
  return 0;
}

/* return the correction of bet for a flat potential-energy histogram */
INLINE double avb_getbetcor(avb_t *avb, double etot)
{
  int i1, i2;
  double bet1, dude1, dude2, cor;

  if (etot < avb->emin || etot >= avb->emax)
    return 0;
  i1 = (int)( (etot - avb->emin) / avb->edel );
  /* NOTE: we differentiate i1 and i1+1, so i1 == n - 1 is unacceptable */
  if (i1 < 0 || i1 >= avb->n - 1) return 0;

  if (avb->avu[i1].s < AVB_MINDATA2) return 0;
  dude1 = avb_getdude_low( &(avb->av[i1]), &(avb->avek[i1]), avb->dof, &bet1);

  /* although we can do the thermal averaging
   * we use discrete differentiation, so we need the second bin */
  i2 = i1+1;
  if (avb->avu[i2].s < AVB_MINDATA2) return 0;
  dude2 = avb_getdude_low( &(avb->av[i2]), &(avb->avek[i2]), avb->dof, NULL );

  /* correction */
  if (dude1 > 0 && dude2 > 0)
    cor = log(dude1/dude2);
  else
    cor = 2.*(dude1 - dude2)/(dude1 + dude2);
  cor /= avb->edel;
  /* we confine the maximal correction
   * strictly, shouldn't do this, but more stable this way */
  cor = dblconfine(cor, -bet1*.3, bet1*.3);
  return cor;
}

/* get the entropy difference: dS = S(e2) - S(e1)
 * *dscor is a correction intended to get the flat potential-energy histogram
 * but it makes unweighted averaging approximate */
INLINE double avb_getds(avb_t *avb, double e1, double e2, double *dscor)
{
  int i1, i2, i, sgn = 1;
  double ds = 0, dude1 = 0, dude2 = 0;

  if (dscor) *dscor = 0; /* zero the correction by default */

  if (e1 > e2) { sgn = -1; dblswap(e1, e2); }
  if (e2 < avb->emin || e1 >= avb->emax)
    return (e2 - e1) * sgn * avb->bet0;

  /* make no correction for E outside of the range */
  if (e1 < avb->emin) { ds += (avb->emin - e1) * avb->bet0; e1 = avb->emin; }
  if (e2 >= avb->emax) { ds += (e2 - avb->emax) * avb->bet0; e2 = avb->emax - 1e-8; }

  i1 = (int)( (e1 - avb->emin) / avb->edel );
  i2 = (int)( (e2 - avb->emin) / avb->edel );
  die_if(i1 < 0 || i1 >= avb->n || i2 < 0 || i2 >= avb->n,
    "i1 %d, i2 %d out of range, n %d, e1 %g, e2 %g\n", i1, i2, avb->n, e1, e2);

  if (i1 == i2) {
    ds += (e2 - e1) * AVB_GETBET(avb, i1);
    /* no correction here */
  } else {
    ds += ((avb->emin + (i1+1)*avb->edel) - e1) * AVB_GETBET(avb, i1);
    ds += (e2 - (avb->emin + i2*avb->edel)) * AVB_GETBET(avb, i2);
    for (i = i1 + 1; i < i2; i++)
      ds += AVB_GETBET(avb, i) * avb->edel;

    /* correction */
    if (dscor != NULL
     && avb_getdude(avb, e1, &dude1) == 0
     && avb_getdude(avb, e2, &dude2) == 0) {
      if (dude1 > 0 && dude2 > 0)
        *dscor = log(dude1/dude2);
      else
        *dscor = 2.*(dude1 - dude2)/(dude1 + dude2);
      ds += *dscor;
    }
  }
  return sgn * ds;
}



/* write the output */
INLINE int avb_write(avb_t *avb, const char *fn)
{
  int i;
  FILE *fp;
  double bet, cnt, bet0, cnt0, uav, tot, ek, dude, scal;

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
    cnt = avb->av[i].s;
    bet0 = av_getave( &(avb->av0[i]) );
    cnt0 = avb->avu[i].s;
    uav = av_getave( &(avb->avu[i]) );
    ek = av_getave( &(avb->avek[i]) );
    dude = avb_getdude_low( &(avb->av[i]), &(avb->avek[i]), avb->dof, NULL);
    fprintf(fp, "%g %g %g %g %g %g %g %g %g\n",
      avb->emin + i * avb->edel,
      bet, cnt, bet0, cnt0*scal, uav, cnt0, ek, dude);
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
  int i, acc = 0;
  real ek1 = *ekin, s;
  double logek1, logek2, ek2, r, ds, dsc = 0, etot1, etot2;

  logek1 = log(ek1);
  logek2 = logek1 + amp*(2.f*rnd0() - 1);
  ek2 = exp(logek2);
  etot1 = ek1 + ep;
  etot2 = ek2 + ep;

  /* causes troubles in LJ with the following limitation
   * but might be okay in other cases */
  if (etotmin < etotmax) {
    if ( (etot1 > etotmin && etot1 < etotmax) /* previously in range */
      && (etot2 > etotmax || etot2 < etotmin)) { /* now not */
      return 0;
    } else if (etot1 < etotmin) {
      acc = (etot2 > etot1);
    } else if (etot1 > etotmin) {
      acc = (etot2 < etot1);
    }
    if (acc) goto ACC;
  }

  ds = avb_getds(avb, etot1, etot2, epcor ? &dsc : 0);
  r = ds - .5*avb->dof*(logek2 - logek1);

  if (r <= 0 || rnd0() < exp(-r)) {
ACC:
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
  double r, sqtp, ds, dsc = 0, bet1, bet2, etot1, etot2, lnrose;

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
  ds = avb_getds(avb, etot1, etot2, epcor ? &dsc : 0);
  r = bet2 * eki1 - bet1 * eki2 + ds + lnrose;

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

