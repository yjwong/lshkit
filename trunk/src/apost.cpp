/* 
    Copyright (C) 2009 Wei Dong <wdong@princeton.edu>. All Rights Reserved.
  
    This file is part of LSHKIT.
  
    LSHKIT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LSHKIT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LSHKIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <queue>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <boost/progress.hpp>
#include <lshkit/apost.h>

namespace lshkit
{

// Gaussian model of nearest neighbor hash components
// this is the gaussian model for one training query point
struct ExampleModel {
    std::vector<float> H;
    std::vector<float> mean;
    std::vector<float> var;

    void estimate (const APostLsh &lsh, const APostExample &example)
    {
        // calculate H values
        const float *query = example.query;
        H.resize(lsh.M);
        for (unsigned i = 0; i < lsh.M; ++i) {
            H[i] = lsh.apply1(query, i);
        }

        // mean & cov matrix of examples
        const std::vector<const float *> &results = example.results;
        std::vector<float> M;
        std::vector<std::vector<float> > cov; // lower matrix
        unsigned dim = lsh.dim;

        BOOST_VERIFY(results.size() > 1);
        // calculate mean of NNs
        M.resize(dim);
        std::fill(M.begin(), M.end(), 0.0);
        BOOST_FOREACH(const float *v, results) {
            for (unsigned i = 0; i < dim; ++i) {
                M[i] += v[i];
            }
        }
        for (unsigned i = 0; i < dim; ++i) {
            M[i] /= results.size();
        }

        // calculate cov of NNs
        cov.resize(dim);
        for (unsigned i = 0; i < dim; i++) {
            cov[i].resize(i+1);
            std::fill(cov[i].begin(), cov[i].end(), 0.0);
        }

        BOOST_FOREACH(const float *v, results) {
            for (unsigned i = 0; i < dim; i++) {
                for (unsigned j = 0; j <= i; j++) {
                    cov[i][j] += (v[i] - M[i]) * (v[j] - M[j]);
                }
            }
        }
        for (unsigned i = 0; i < dim; i++) {
            for (unsigned j = 0; j <= i; j++) {
                cov[i][j] /= (results.size() - 1);
            }
        }

        // mean & std of hash values
        mean.resize(lsh.M);
        var.resize(lsh.M);
        for (unsigned i = 0; i < lsh.M; i++) { // for each hash component
            mean[i] = lsh.apply1(&M[0], i);
            var[i] = 0;
            for (unsigned ii = 0; ii < lsh.M; ii++) {
                unsigned jj;
                for (jj = 0; jj <= ii; jj++) {
                    var[i] += lsh.a[i][ii] * cov[ii][jj] * lsh.a[i][jj];
                }
                for (; jj < lsh.M; jj++) {
                    var[i] += lsh.a[i][ii] * cov[jj][ii] * lsh.a[i][jj];
                }
            }
            var[i] /= sqr(lsh.W);
        }
    }
};

class GaussianHashModel
{
    std::vector<ExampleModel> M;
    float sigma;
public:
    GaussianHashModel (const APostLsh &lsh, const std::vector<APostExample> &examples, float s = 0.0) {
        if (s == 0.0) {
            s = 1.0 / 5;
            // I believe there's a bug in paper.
            // Right below equation 16, the paper says the default
            // sigma is W/5.  However,
            // the h values are already divided by W (h = (a v + b) / W)
            // and it makes no sense to divide it by W again.
        }
        sigma = s;

        M.resize(examples.size());
        for (unsigned i = 0; i < examples.size(); ++i) {
            M[i].estimate(lsh, examples[i]);
        }
    }

    void estimate (unsigned m, float h, float *mean, float *std) { // mth hash component
        float mm = 0, vv = 0, ss = 0, k;
        for (unsigned i = 0; i < M.size(); i++) { // for each example
            k = gsl_ran_gaussian_pdf(M[i].H[m] - h, sigma);
            mm += k * M[i].mean[m];
            vv += k * M[i].var[m];
            ss += k;
        }
        mm /= ss;
        vv /= ss;
        *mean = mm;
        *std = sqrt(vv);
    }
};

static inline float GaussianInterval (float mean, float std, float l, float u) {
    return gsl_cdf_gaussian_P(u - mean, std)
        - gsl_cdf_gaussian_P(l - mean, std);
}

void APostModel::train (const APostLsh &lsh,
                        const std::vector<APostExample> &examples,
                        unsigned N, float k_sigma, float expand) {
    Nz = N;
    ex = expand;
    GaussianHashModel parzen(lsh, examples, k_sigma);
    umin = lsh.umin;
    umax = lsh.umax;

    lookup.resize(lsh.M);
    means.resize(lsh.M);
    stds.resize(lsh.M);
    // for each component h
    boost::progress_display progress(lsh.M * Nz);
    for (unsigned m = 0; m < lsh.M; ++m) {
        {
            float ex = expand * (umax[m] - umin[m]);
            umin[m] -= ex;
            umax[m] += ex;
        }

        float delta = (umax[m] - umin[m]) / Nz;
        int minh(std::floor(umin[m]));
        int maxh(std::floor(umax[m])); // min & max inclusive
        unsigned size = maxh - minh + 1;

        lookup[m].resize(Nz);
        means[m].resize(Nz);
        stds[m].resize(Nz);
        // for each quantum of h(q)
        for (unsigned n = 0; n < Nz; ++n) {
            float mean, std;
            parzen.estimate(m, umin[m] + (n + 0.5) * delta, &mean, &std);
            means[m][n] = mean;
            stds[m][n] = std;
            std::vector<PrH> &lmn = lookup[m][n];
            lmn.resize(size);
            for (int h = 0; h < size; ++h) {
                lmn[h].h = (unsigned)(minh + h);
                lmn[h].pr = GaussianInterval(mean, std, minh + h, minh + h + 1);
//            std::cout << "<" << lmn[h].h << ", " << lmn[h].pr << ">" << std::endl;
            }
            std::sort(lmn.begin(), lmn.end());
            ++progress;
        }
    }
}

struct PrC {
    unsigned m;
    const std::vector<PrH> *prh;
    friend bool operator < (const PrC &c1, const PrC &c2) {
        return c1.prh->at(0).pr > c2.prh->at(0).pr;
    }
};

struct Probe {
    std::vector<unsigned> off;
    unsigned last;
    float pr;

    Probe () {
        BOOST_VERIFY(0);
    }

    Probe (unsigned m) {
        off.resize(m);
        std::fill(off.begin(), off.end(), 0);
        last = 0;
    }

    bool canShift () const {
        if (off[last] != 1) return false;
        if (last + 1 >= off.size()) return false;
        return true;
    }

    void shift () {
        off[last] = 0;
        last++;
        off[last] = 1;
    }

    bool canExpand () const {
        if (last + 1 >= off.size()) return false;
        return true;
    }

    void expand () {
        last++;
        off[last] = 1;
    }

    void extend () {
        off[last]++;
    }

    friend bool operator < (const Probe &p1, const Probe &p2) {
        return p1.pr < p2.pr;
    }

    void setPr (const std::vector<PrC> &pl) {
        pr = 1.0;
        for (unsigned i = 0; i < off.size(); ++i) {
            pr *= pl[i].prh->at(off[i]).pr;
        }
//        std::cerr << "PR = " << pr << std::endl;
    }

    unsigned hash (const APostLsh &lsh,
                    const std::vector<PrC> &pl) {
        unsigned r = 0;
        for (unsigned i = 0; i < lsh.M; ++i) {
            r += lsh.c[pl[i].m] * pl[i].prh->at(off[i]).h;
        }
        return r % lsh.H;
    }

    unsigned print (const std::vector<PrC> &pl) {
        BOOST_FOREACH(unsigned v, off) {
            std::cout << ' ' << v;
        }
        std::cout << std::endl;
    }
};

void APostModel::genProbeSequence (const APostLsh &lsh,
                                APostLsh::Domain query,
                                float recall, unsigned T,
                                std::vector<unsigned> *seq) const
{
    std::cout << "Gaussian:" << std::endl;
    std::vector<PrC> pl(lsh.M);
    for (unsigned i = 0; i < lsh.M; ++i) {
        pl[i].m = i;

        float h = lsh.apply1(query, i);

        if (h < umin[i]) {
            std::cerr << "hash[" << i << "] out of range " << h << " < umin = " << umin[i] << std::endl;
            h = umin[i];
        }
        else if (h > umax[i]) {
            std::cerr << "hash[" << i << "] out of range " << h << " > umax = " << umax[i] << std::endl;
            h = umax[i];
        }

        unsigned qh = (h - umin[i]) * Nz / ( umax[i] - umin[i]);
        BOOST_VERIFY(qh < lookup[i].size());

        pl[i].prh = &lookup[i][qh];

        std::cout << i << ':' << h << ':' << qh << '\t' << umin[i] << ':' << umax[i] << '\t' << means[i][qh]
            << ':' << stds[i][qh] << std::endl;
    }

    std::cout << "PROB:" << std::endl;
    std::sort(pl.begin(), pl.end());

    BOOST_FOREACH(const PrC &prc, pl) {
        std::cout << prc.m << ":";
        BOOST_FOREACH(PrH prh, *prc.prh) {
            if (prh.pr == 0) break;
            std::cout << ' ' << prh.pr << '@' << (int)prh.h;
        }
        std::cout << std::endl;
    }

    std::cout << "SEQ" << std::endl;
    // generate probe sequence
    seq->clear();

    std::vector<Probe> heap;
    heap.push_back(Probe(lsh.M));

    heap.back().setPr(pl);
    float pr = heap.back().pr;
    seq->push_back(heap.back().hash(lsh, pl));

    if (pr >= recall) return;
    if (seq->size() >= T) return;

    heap.back().off[0] = 1;
    heap.back().setPr(pl);

    for (;;) {
        if (pr >= recall) break;
        if (seq->size() >= T) break;
        if (heap.empty()) break;

        pop_heap(heap.begin(), heap.end());

        seq->push_back(heap.back().hash(lsh, pl));
        pr += heap.back().pr;
        std::cout << pr << ", " << heap.back().pr << ":";
        heap.back().print(pl);

        Probe p = heap.back();

        {
            heap.back().extend();
            heap.back().setPr(pl);
            push_heap(heap.begin(), heap.end());
        }

        if (p.canShift()) {
            heap.push_back(p);
            heap.back().shift();
            heap.back().setPr(pl);
            push_heap(heap.begin(), heap.end());
        }

        if (p.canExpand()) {
            heap.push_back(p);
            heap.back().expand();
            heap.back().setPr(pl);
            push_heap(heap.begin(), heap.end());
        }
    }
    std::cout << "DONE" << std::endl;
}

}
