void intdeoini(int lenaw, double tiny, double eps, double *aw)
{
    /* ---- adjustable parameter ---- */
    int lmax = 5;
    double efs = 0.1, enoff = 0.40, pqoff = 2.9, ppoff = -0.72;
    /* ------------------------------ */
    int noff0, nk0, noff, k, nk, j;
    double pi4, tinyln, epsln, frq4, per2, pp, pq, ehp, ehm, h, t, 
        ep, em, tk, xw, wg, xa;
    
    pi4 = atan(1.0);
    tinyln = -log(tiny);
    epsln = 1 - log(efs * eps);
    frq4 = 1 / (2 * pi4);
    per2 = 4 * pi4;
    pq = pqoff / epsln;
    pp = ppoff - log(pq * pq * frq4);
    ehp = exp(2 * pq);
    ehm = 1 / ehp;
    aw[3] = lmax;
    aw[4] = eps;
    aw[5] = sqrt(efs * eps);
    noff0 = 6;
    nk0 = 1 + (int) (enoff * epsln);
    aw[1] = nk0;
    noff = 2 * nk0 + noff0;
    wg = 0;
    xw = 1;
    for (k = 1; k <= nk0; k++) {
        wg += xw;
        aw[noff - 2 * k] = wg;
        aw[noff - 2 * k + 1] = xw;
        xw = xw * (nk0 - k) / k;
    }
    wg = per2 / wg;
    for (k = noff0; k <= noff - 2; k += 2) {
        aw[k] *= wg;
        aw[k + 1] *= wg;
    }
    xw = exp(pp - 2 * pi4);
    aw[noff] = sqrt(xw * (per2 * 0.5));
    aw[noff + 1] = xw * pq;
    aw[noff + 2] = per2 * 0.5;
    h = 2;
    nk = 0;
    k = noff + 3;
    do {
        t = h * 0.5;
        do {
            em = exp(2 * pq * t);
            ep = pi4 * em;
            em = pi4 / em;
            tk = t;
            j = k;
            do {
                xw = exp(pp - ep - em);
                wg = sqrt(frq4 * xw + tk * tk);
                xa = xw / (tk + wg);
                wg = (pq * xw * (ep - em) + xa) / wg;
                aw[j] = xa;
                aw[j + 1] = xw * pq;
                aw[j + 2] = wg;
                ep *= ehp;
                em *= ehm;
                tk += 1;
                j += 3;
            } while (ep < tinyln && j <= lenaw - 3);
            t += h;
            k += nk;
        } while (t < 1);
        h *= 0.5;
        if (nk == 0) {
            if (j > lenaw - 6) j -= 3;
            nk = j - noff;
            k += nk;
            aw[2] = nk;
        }
    } while (2 * k - noff - 3 <= lenaw);
    aw[0] = k - 3;
}

void intdeo(double (*f)(double, double), double a, double omega, double *aw, 
    double *i, double *err, double x1)
{
    int lenawm, nk0, noff0, nk, noff, lmax, m, k, j, jm, l;
    double eps, per, perw, w02, ir, h, iback, irback, t, tk, 
        xa, fm, fp, errh, s0, s1, s2, errd;
 
	jm = 0; fm = 0.0; fp = 0.0; errh = 0.0; 
	 
    lenawm = (int) (aw[0] + 0.5);
    nk0 = (int) (aw[1] + 0.5);
    noff0 = 6;
    nk = (int) (aw[2] + 0.5);
    noff = 2 * nk0 + noff0;
    lmax = (int) (aw[3] + 0.5);
    eps = aw[4];
    per = 1 / fabs(omega);
    w02 = 2 * aw[noff + 2];
    perw = per * w02;
    *i = (*f)(a + aw[noff] * per, x1);
    ir = *i * aw[noff + 1];
    *i *= aw[noff + 2];
    *err = fabs(*i);
    h = 2;
    m = 1;
    k = noff;
    do {
        iback = *i;
        irback = ir;
        t = h * 0.5;
        do {
            if (k == noff) {
                tk = 1;
                k += nk;
                j = noff;
                do {
                    j += 3;
                    xa = per * aw[j];
                    fm = (*f)(a + xa, x1);
                    fp = (*f)(a + xa + perw * tk, x1);
                    ir += (fm + fp) * aw[j + 1];
                    fm *= aw[j + 2];
                    fp *= w02 - aw[j + 2];
                    *i += fm + fp;
                    *err += fabs(fm) + fabs(fp);
                    tk += 1;
                } while (aw[j] > eps && j < k);
                errh = *err * aw[5];
                *err *= eps;
                jm = j - noff;
            } else {
                tk = t;
                for (j = k + 3; j <= k + jm; j += 3) {
                    xa = per * aw[j];
                    fm = (*f)(a + xa, x1);
                    fp = (*f)(a + xa + perw * tk, x1);
                    ir += (fm + fp) * aw[j + 1];
                    fm *= aw[j + 2];
                    fp *= w02 - aw[j + 2];
                    *i += fm + fp;
                    tk += 1;
                }
                j = k + jm;
                k += nk;
            }
            while (fabs(fm) > *err && j < k) {
                j += 3;
                fm = (*f)(a + per * aw[j], x1);
                ir += fm * aw[j + 1];
                fm *= aw[j + 2];
                *i += fm;
            }
            fm = (*f)(a + perw * tk, x1);
            s2 = w02 * fm;
            *i += s2;
            if (fabs(fp) > *err || fabs(s2) > *err) {
                l = 0;
                for (;;) {
                    l++;
                    s0 = 0;
                    s1 = 0;
                    s2 = fm * aw[noff0 + 1];
                    for (j = noff0 + 2; j <= noff - 2; j += 2) {
                        tk += 1;
                        fm = (*f)(a + perw * tk, x1);
                        s0 += fm;
                        s1 += fm * aw[j];
                        s2 += fm * aw[j + 1];
                    }
                    if (s2 <= *err || l >= lmax) break;
                    *i += w02 * s0;
                }
                *i += s1;
                if (s2 > *err) *err = s2;
            }
            t += h;
        } while (t < 1);
        if (m == 1) {
            errd = 1 + 2 * errh;
        } else {
            errd = h * (fabs(*i - 2 * iback) + fabs(ir - 2 * irback));
        }
        h *= 0.5;
        m *= 2;
    } while (errd > errh && 2 * k - noff <= lenawm);
    *i *= h * per;
    if (errd > errh) {
        *err = -errd * per;
    } else {
        *err *= per * m * 0.5;
    }
}
