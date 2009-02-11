/* 
    Copyright (C) 2008 Wei Dong <wdong@princeton.edu>. All Rights Reserved.
  
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

/**
  * \file sketch-run.cpp
  * \brief Example of sketch construction with Gaussian LSH.
  * 
  * This program uses sketch filtering to accelerate K-NN search.
  * The idea is to first search against a dataset of sketches and keep the
  * top C*K points as candidates.  The candidates are than ranked using
  * the raw feature vectors.  The sketch database can be viewed as a index.
  * 
  * In this program, we use 2-stable LSH based sketch.  Each sketch is a bit-vector
  * of M bits, and each bit is produced by a independent hash function from the
  * family DeltaLSB<GaussianLsh>.
  *
  * The program reconstruct the sketches by default.  You can give the
  * --index option to make the program save the sketches.  The next
  * time you run the program with the same --index option, the program
  * will try to load the previously saved sketches.  When a saved sketch database is
  * used, you need to make sure that the dataset and other parameters match
  * the previous run.  However, the benchmark file, Q and K can be different.
\verbatim
Allowed options:
  -h [ --help ]          produce help message.
  -W [ -- ] arg (=1)
  -M [ -- ] arg (=1)     skech size / byte
  -C [ -- ] arg (=10)    # candidates = C x K
  -Q [ -- ] arg (=100)   # queries to use
  -K [ -- ] arg (=50)    K-NNs retrieved
  -D [ --data ] arg      data file
  -B [ --benchmark ] arg benchmark file
  --index arg            sketch file
  --asym                 Asymmetric distance estimation
\endverbatim
  */

#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/timer.hpp>
#include <lshkit.h>

using namespace std;
using namespace lshkit;
namespace po = boost::program_options; 

int main (int argc, char *argv[])
{
    string data_file;
    string index_file;
    string benchmark;

    float W;
    unsigned M;
    unsigned Q, K, C;
    bool do_benchmark = true;
    bool load_index = false; // load the index from a file
    bool save_index = false; // save the index to file
    bool do_asym = false;

    boost::timer timer;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message.")
        (",W", po::value(&W)->default_value(1.0), "")
        (",M", po::value(&M)->default_value(1), "skech size / byte")
        (",C", po::value(&C)->default_value(10), "# candidates = C x K")
        (",Q", po::value(&Q)->default_value(100), "# queries to use")
        (",K", po::value(&K)->default_value(50), "K-NNs retrieved")
        ("data,D", po::value(&data_file), "data file")
        ("benchmark,B", po::value(&benchmark), "benchmark file")
        ("index", po::value(&index_file), "sketch file")
        ("asym", "Asymmetric distance estimation")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm); 

    if (vm.count("help") || (vm.count("data") < 1))
    {
        cout << desc;
        return 0;
    }

    if ((Q == 0) || (vm.count("benchmark") == 0)) {
        do_benchmark = false;
    }

    if (vm.count("asym") == 1) {
        do_asym = true;
    }

    if (vm.count("index") == 1) {
        if (boost::filesystem::exists(index_file)) {
            load_index = true;
        } else {
            save_index = true;
        }
    }

    cout << "LOADING DATA..." << endl;
    timer.restart();
    FloatMatrix data(data_file);
    cout << boost::format("LOAD TIME: %1%s.") % timer.elapsed() << endl;

    typedef Sketch<DeltaLSB<GaussianLsh> > MySketch;

    MySketch sketcher;
    Matrix<unsigned char> index;

    if (load_index) {
        cout << "LOADING INDEX..." << endl;
        timer.restart();
        {
            ifstream is(index_file.c_str(), ios_base::binary);
            is.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);
            sketcher.load(is);
            index.load(is);
        }
        cout << boost::format("LOAD TIME: %1%s.") % timer.elapsed() << endl;
        verify(M <= sketcher.getChunks());
    }
    else {
        // We define a short name for the MPLSH index.
        MySketch::Parameter param;

        // Setup the parameters.  Note that L is not provided here.
        param.W = W;
        param.dim = data.getDim();
        DefaultRng rng;

        sketcher.reset(M, param, rng);
        // The accessor.

        // Initialize the index structure.  Note L is passed here.
        cout << "CONSTRUCTING INDEX..." << endl;

        index.reset(M, data.getSize());

        // CONSTRUCT SKETCHES FOR DATABASE
        timer.restart();
        {
            boost::progress_display progress(data.getSize());
            for (unsigned i = 0; i < data.getSize(); ++i)
            {
                sketcher.apply(data[i], index[i]);
                ++progress;
            }
        }
        cout << boost::format("CONSTRUCTION TIME: %1%s.") % timer.elapsed() << endl;

        if (save_index) {
            timer.restart();
            cout << "SAVING INDEX..." << endl;
            {
                ofstream os(index_file.c_str(), ios_base::binary);
                os.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);
                sketcher.save(os);
                index.save(os);
            }
            cout << boost::format("SAVING TIME: %1%s") % timer.elapsed() << endl;
        }
    }

    verify(index.getSize() == data.getSize());

    if (do_benchmark) {

        Benchmark<> bench;
        bench.resize(K, Q);
        cout << "LOADING BENCHMARK..." << endl;
        bench.load(benchmark);
        cout << "DONE." << endl;

        for (unsigned i = 0; i < Q; ++i)
        {
            for (unsigned j = 0; j < K; ++j)
            {
                assert(bench.getAnswer(i)[j].key < data.getSize());
            }
        }

        cout << "RUNNING QUERIES..." << endl;

        Stat recall;
        Topk<unsigned> candidate;
        Topk<unsigned> topk;
        unsigned char *query_sketch = new unsigned char[sketcher.getChunks()];
        boost::progress_display progress(Q);
        metric::l2<float> l2(data.getDim());

        timer.restart();
        if (do_asym)
        {
            float *asym = new float[sketcher.getBits()];
            WeightedHammingHelper<> helper(M);
            for (unsigned i = 0; i < Q; ++i)
            {
                // Query for one point.
                unsigned query = bench.getQuery(i);
                candidate.reset(C * K);
                topk.reset(K);

                sketcher.apply(data[query], query_sketch, asym);
                helper.update(query_sketch, asym);
                // Scan the sketches for candidates.
                for (unsigned j = 0; j < index.getSize(); ++j) {
                    if (j == query) continue;
                    Topk<unsigned>::Element e;
                    e.key = j;
                    e.dist = helper.distTo(index[j]);
                    candidate << e;
                }
                // Rank with raw feature.
                for (unsigned j = 0; j < candidate.size(); ++j) {
                    Topk<unsigned>::Element e;
                    e.key = candidate[j].key;
                    e.dist = l2(data[query], data[e.key]);
                    topk << e;
                }
                recall << bench.getAnswer(i).recall(topk);
                ++progress;
            }
            delete []asym;
        }
        else
        // symmetric sketch
        {
            metric::hamming<unsigned char> hamming(M);
            for (unsigned i = 0; i < Q; ++i)
            {
                // Query for one point.
                unsigned query = bench.getQuery(i);
                candidate.reset(C * K);
                topk.reset(K);

                sketcher.apply(data[query], query_sketch);
                for (unsigned j = 0; j < index.getSize(); ++j) {
                    if (j == query) continue;
                    Topk<unsigned>::Element e;
                    e.key = j;
                    e.dist = hamming(index[j], query_sketch);
                    candidate << e;
                }
                for (unsigned j = 0; j < candidate.size(); ++j) {
                    Topk<unsigned>::Element e;
                    e.key = candidate[j].key;
                    e.dist = l2(data[query], data[e.key]);
                    topk << e;
                }
                recall << bench.getAnswer(i).recall(topk);
                ++progress;
            }
        }

        delete []query_sketch;
        cout << boost::format("QUERY TIME: %1%s.") % timer.elapsed() << endl;
        cout << "[RECALL] " << recall.getAvg() << " +/- " << recall.getStd() << endl;

    }

    return 0;
}

