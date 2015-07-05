#include <map>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <unordered_map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "deftypes.hpp"
#include "jsoncons/json.hpp"
using namespace std;

#ifndef XFILE_COMPATIBLE
typedef cities* const cities_t;
#endif

typedef const jsoncons::json& config_t;
void load_cities(jsoncons::json&,
#ifndef XFILE_COMPATIBLE
    cities_t
#else
    int**&
#endif
    ) __nonnull();
void statistics(POPULATION* const, IPTR) __nonnull();
void init_report(
#ifdef XFILE_COMPATIBLE
                const char*,
#endif
                const POPULATION* const) __nonnull();
void init_pop(POPULATION* const,
#ifndef XFILE_COMPATIBLE
    const cities_t
#else
    const int** const
#endif
) __nonnull();
void raw_stat(FILE*, const POPULATION* const, IPTR) __nonnull();
void config_pop(config_t, POPULATION* const) __nonnull();
#ifndef XFILE_COMPATIBLE
double eval(const POPULATION* const, const cities* const, IPTR const) __nonnull();
#else
double eval(const POPULATION* const, const int** const, IPTR const) __nonnull();
#endif
void initialize(jsoncons::json& conf, POPULATION* const p,
#ifndef XFILE_COMPATIBLE
    cities_t c
#else
    int**& c
#endif
)
{
    load_cities(conf, c);
    printf("%s data loaded.\n", __PROMT_PASSED);
    config_pop(conf, p);
    printf("%s populaiton configured.\n", __PROMT_PASSED);
    init_pop(p,
#ifndef XFILE_COMPATIBLE
        c
#else
        const_cast<const int** const>(c)
#endif
        );
    printf("%s population initialized.\n", __PROMT_PASSED);
#ifdef STAT_ENABLED
    statistics(p, p->op);
    printf("%s statistical info. gathered from population.\n", __PROMT_PASSED);
    printf("\n-----------------------configs-----------------------\n\n");
    init_report(
#ifdef XFILE_COMPATIBLE
        conf["output-file"].as_string().c_str(),
#endif
        p);
#endif
}
void load_cities(jsoncons::json& conf,
#ifndef XFILE_COMPATIBLE
    cities_t c
#else
    int**& dist
#endif
    ) {
#ifndef XFILE_COMPATIBLE
    auto file = fopen(conf["input-file"].as_string().c_str(), "r");
    city tmp;
    int stat, continue_count = 0;
    while(true) {
        stat = fscanf(file, "%lf%*[ \t,]%lf", &tmp.x, &tmp.y);
        if(stat == EOF || continue_count > 100) break;
        if(stat != 2 && ++continue_count)   continue;
        c->push_back(tmp);
    }
    fclose(file);
    if(c->size() == 0) throw runtime_error("no data loaded");
    conf["chrom_size"] = c->size();
#else
    ifstream file(conf["input-file"].as_string());
    string line;
    // skip the configuration part
    while(true) {
        file >> line;
        if(line == conf["output-file"].as_string() && std::getline(file, line))
            break;
    }
    const size_t lchrom = conf["chrom_size"].as<size_t>();
    dist = new int*[lchrom];
    for(size_t i = 0; i < lchrom; i++) dist[i] = new int[lchrom];
    size_t lineNumber = -1;
    // has reached to the data part
    while(!file.eof() && ++lineNumber < lchrom) {
        // read the whole line and trim it
        std::getline(file, line), line = boost::trim_copy(line);
        // validate the input
        if(line.length() == 0) break;
        // split the line into its token
        vector<string> distances; boost::split(distances, line, boost::is_any_of(" "));
        // fail-safe check point
        assert(distances.size() == lchrom);
        // add the read distances into distance matrix
        for(size_t i = 0; i < distances.size(); i++)
        {
            try { dist[lineNumber][i] = boost::lexical_cast<int>(distances[i].c_str()); }
            catch(boost::bad_lexical_cast &) { throw runtime_error("Invalid input data"); }
        }
    }
#endif
}

void config_pop(config_t conf, POPULATION* const p)
{
    p->lchrom               = conf["chrom_size"].as<int>();
#ifdef XFILE_COMPATIBLE
    p->ofile                = const_cast<char*>(conf["output-file"].as_string().c_str());
#else
    p->ofile                = const_cast<char*>("log.out");
#endif
    p->pElit                = conf["elit_prop"].as_double();
    p->smallestEverGen      = std::numeric_limits<int>::max();
    p->smallestEverIndex    = std::numeric_limits<int>::max();
    p->popSize              = conf["population_size"].as_int();
    p->maxGen               = conf["max_generations"].as_int();
    p->pMut                 = conf["mutation_prob"].as_double();
    p->pCross               = conf["crossover_prob"].as_double();
    p->smallestEverFitness  = std::numeric_limits<double>::max();
    p->cityfile             = const_cast<char*>(conf["output-file"].as_string().c_str());
    if(p->popSize % 2 != 0) p->popSize++;
}
#ifndef DISABLE_INIT_CBLOCK
#ifndef XFILE_COMPATIBLE
map<size_t, vector<size_t>> make_cblocks(
    const cities_t c,
    double factor,
    double init_range,
    bool (*cmp)(const city&, const city&),
    bool (*ranger)(const double&, const city&))
{
    auto xc = c;
    std::sort(xc->begin(), xc->end(), cmp);
    double range = init_range;
    int j = 1;
    map<size_t, vector<size_t>> k;
    for(size_t i = 0; i < xc->size(); i++) {
        if(ranger(range + factor, xc->at(i))) {
            j++;
            range += (factor);
        }
        k[j].push_back(i);
    };
    return k;
}
#endif
#endif

const vector<int> pre_init_pop(POPULATION* const p,
#ifndef XFILE_COMPATIBLE
    const cities_t c
#else
    const int** const c
#endif
    , size_t* const index) {
    *index = 0;
#ifndef DISABLE_INIT_CBLOCK
#ifndef XFILE_COMPATIBLE
    double s = 0, _max = 0, _min = std::numeric_limits<double>::max();
    for(int i = 0; i < p->lchrom; i++) {
        auto o = (c->at(i).x + c->at(i).y); s += o;
        if(c->at(i).x > _max) _max = o;
        if(c->at(i).x < _min) _min = o;
    }
    s /= p->lchrom;
    double v = 0;
    for(int i = 0; i < p->lchrom; i++) {
        auto x = (c->at(i).x + c->at(i).y);
        v += pow((x - s), 2.0);
    }
    v /= (p->lchrom * _max + _min);
    auto k0 = make_cblocks(c, v, 0, [](const city& a, const city& b){
            return a.x + a.y < b.x + b.y;
    }, [](const double& x, const city& _c) { return x < _c.x + _c.y; });
    auto k1 = make_cblocks(c, v, 0, [](const city& a, const city& b){
        return a.y < b.y;
    }, [](const double& x, const city& _c) { return x < _c.y; });
    auto k2 = make_cblocks(c, v, 0, [](const city& a, const city& b){
        return a.x < b.x;
    }, [](const double& x, const city& _c) { return x < _c.x; });
    map<size_t, size_t> ll;
    map<size_t, vector<size_t>> ff;
    decltype(k1)* ka[] = {&k0, &k1, &k2};
    for(size_t i = 0; i < sizeof(ka) / sizeof(decltype(k1)*); i++) {
        auto k = *ka[i];
        for(auto& m : k) {
            for(auto& n: m.second) {
                ll[n] += m.first;
            }
        }
    }
    for(auto& b: ll) { ff[b.second].push_back({b.first}); }
    for(auto& b: ff) {
        if(b.second.size() < 2) {
            size_t cl, cc = 1;
            do {
                cl = b.first + cc++;
                if(ll[p->lchrom-1] == b.second.front()) cl = b.first - cc++;
            } while(!(ff.count(cl) && ff[cl].size()));
            ff[cl].push_back(b.second.front());
            ll[b.second.front()] = cl;
            b.second.clear();
        }
    }
    for (size_t j = 0; *index < sqrt(p->popSize); (*index)++, j = 0){
        IPTR po = &(p->op[*index]);
        po->chrom  = (decltype(po->chrom))malloc(p->lchrom * sizeof(int));
        for(auto& lff: ff) {
            auto clone = lff.second;
            if(!clone.size()) continue;
            std::random_shuffle(clone.begin(), clone.end());
            for(auto& iin : clone)
                po->chrom[j++] = iin;
        }
        po->fitness = eval(p, c, po);
    }
#if 0
    ofstream os("cblock.out");
    for(auto& b: ll) {
        auto _c = c->at(b.first);
        cout<<b.first<<" "<<_c.x<<" "<<_c.y<<" "<<b.second<<endl;
        os<<_c.x<<" "<<_c.y<<" "<<b.second<<endl;
    }
    os.close();
    exit(0);
#endif
#endif
#endif
    vector<int> genes; for(int gene = 0; gene < p->lchrom; gene++) genes.push_back(gene);
#ifdef XFILE_COMPATIBLE
    for (; *index < ceil(sqrt(p->popSize)); (*index)++) {
        IPTR po = &(p->op[*index]);
        po->chrom  = (decltype(po->chrom))malloc(p->lchrom * sizeof(int));
        for(int i = 0; i < 3; i++) std::random_shuffle(genes.begin(), genes.end());
        size_t index = genes[0];
        po->chrom[0] = index;
        genes.erase(genes.begin());
        for(int i = 1, counter = 0; genes.size(); counter++, i++) {
            po->chrom[i] = genes.front();
            auto xk = c[index][po->chrom[i]];
            genes.erase(genes.begin());
            if(xk == -1 && counter < 100) { genes.push_back(po->chrom[i]); i -= 1; continue; }
            index = po->chrom[i];
        }
        po->fitness = eval(p, c, po);
        for(int gene = 0; gene < p->lchrom; gene++) genes.push_back(gene);
    }
#endif
    return genes;
}

void init_pop(POPULATION* const p,
#ifndef XFILE_COMPATIBLE
    const cities_t c
#else
    const int** const c
#endif
)
{
    p->op = (IPTR) malloc (p->popSize * sizeof(INDIVIDUAL));
    p->np = (IPTR) malloc (p->popSize * sizeof(INDIVIDUAL));
    size_t i = 0;
    const vector<int> genes = pre_init_pop(p, c, &i);
    for (size_t j = 0; i < (size_t)p->popSize; i++, j = 0){
        IPTR po = &(p->op[i]);
        po->chrom  = (decltype(po->chrom))malloc(p->lchrom * sizeof(int));
        auto clone = genes;
        // shuffle the default solution, to generate new solution
        std::random_shuffle(clone.begin(), clone.end());
//        std::shuffle(clone.begin(), clone.end(), mt19937(dynamic_seed));
#ifdef CHECKSUM_IN_EFFECT
        assert((size_t)p->lchrom == clone.size());
#endif
        for(auto& i: clone)
            po->chrom[j++] = i;
        po->fitness = eval(p, c, po);
    }
#ifdef CHECKSUM_IN_EFFECT
    assert((size_t)p->lchrom == c->size());
    for(int i = 0; i<p->popSize; i++) {
        for(int j = 0; j<p->popSize; j++)
            assert(
                (size_t)p->op[i].chrom[j] < c->size());
    }
#endif
}
void init_report(
#ifdef XFILE_COMPATIBLE
    const char* ofile,
#endif
    const POPULATION* const p)
{
    FILE *fp;
    printf("Population Size             ( popsize )    %d\n", p->popSize);
    printf("Chromosome length           ( lchrom  )    %d\n", p->lchrom);
    printf("Crossover Probability       ( pcross  )    %lf\n", p->pCross);
    printf("Mutation Probability        ( pmut    )    %lf\n", p->pMut);
    printf("Maximum num of generations  ( maxgen  )    %d\n\n\n", p->maxGen);
    if((fp = fopen(
#ifdef XFILE_COMPATIBLE
        ofile
#else
        p->ofile
#endif
        , "w+")) == NULL)
    { fprintf(stderr, "error in opening file %s\n", p->ofile); exit(EXIT_FAILURE); }
#ifndef XFILE_COMPATIBLE
#   define XEF_HEADER " ", "SEG.", " ", "SEF.", " ", "SEI.", " "
#else
#   define XEF_HEADER " ", "LEG.", " ", "LEF.", " ", "LEI.", " "
#endif
#define header_format "%-10s%s%-7s%s%-8s%s%-8s%s%-10s%s%-11s%s%-13s%s%-12s\n", \
    " ", "gen#", " ", "max. fitness", " ", "avg. fitness", " ", "min. fitness", \
    XEF_HEADER
    // this header is pain in ass if want to plot log file's data
    /*fprintf(fp, header_format);*/
    raw_stat(fp, p, p->op); fclose(fp);
#ifdef OS_UNIX
    printf("\033[7m");
#endif
    fprintf(stdout, header_format);
#ifdef OS_UNIX
    printf("\033[m");
#endif
    raw_stat(stdout, p, p->op);
#undef header_format
#undef XEF_HEADER
}
