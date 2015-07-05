#ifndef GACONFIG_H
#define GACONFIG_H
#include "stdafx.hpp"
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
typedef double prob;
typedef double limit;
typedef prob  _ratio;
namespace TSPGA {
#define VALIDATE_PROB(x) if(x < 0 || x > 1) throw std::range_error("Expecting the probability be in range of [0, 1]")
    class gaConfig
    {
        typedef std::string mkey;
        typedef std::string mval;
        typedef std::unordered_map<mkey, mval> umap;
    private:
        umap    map;
        _ratio   elitism_ratio;
        prob    crossover_prob, mutation_prob;
        limit   genes_value_min, genes_value_max;
        size_t  population_size, genes_count_min, genes_count_max, gen_max_count;
    public:
        gaConfig                                (size_t population_size, prob crossover_prob, prob mutation_prob)
                                                                            { this->setPopulation_size(population_size); this->setCrossover_prob(crossover_prob); this->setMutation_prob(mutation_prob); }
        inline void     setGenes_ValueLimit     (limit min, limit max)      { if(max < min) throw std::invalid_argument("Min value cannot exceed fromo Min value"); this->genes_value_min = min; this->genes_value_max = max; }
        inline void     setGenes_CountLimit     (size_t min, size_t max)    { if(max < min) throw std::invalid_argument("Min value cannot exceed fromo Min value"); this->genes_count_min = min; this->genes_count_max = max; }
        inline void     setGenerationMaxCount   (size_t gen_max_count)      { this->gen_max_count = gen_max_count; }
        inline void     setPopulation_size      (size_t ps)                 { this->population_size = ps; }
        inline void     setCrossover_prob       (prob cp)                   { VALIDATE_PROB(cp); this->crossover_prob   = cp; }
        inline void     setMutation_prob        (prob mp)                   { VALIDATE_PROB(mp); this->mutation_prob    = mp; }
        inline void     setElitism_ratio        (_ratio er)                  { VALIDATE_PROB(er); this->elitism_ratio    = er; }
        inline size_t   getPopulation_Size      () const                    { return this->population_size;}
        inline prob     getCrossover_prob       () const                    { return this->crossover_prob; }
        inline prob     getMutation_prob        () const                    { return this->mutation_prob;  }
        inline _ratio   getElitism_ratio        () const                    { return this->elitism_ratio; }
        inline limit    getGenes_Value_MAX_Limit() const                    { return this->genes_value_max;  }
        inline limit    getGenes_Value_MIN_Limit() const                    { return this->genes_value_min;  }
        inline limit    getGenes_Count_MAX_Limit() const                    { return this->genes_count_max;}
        inline limit    getGenes_Count_MIN_Limit() const                    { return this->genes_count_min;}
        inline size_t   getGenerationMaxCount   () const                    { return this->gen_max_count;  }
        inline void     option                  (mkey k, mval v)            { this->map[k] = v; }
        inline mval     option                  (mkey k) const              { return this->map.at(k); }
        inline bool     option_exists           (mkey k) const              { return this->map.count(k); }
        inline size_t   options_count            ()       const              { return this->map.size(); }
        inline const umap& getOptions            ()       const              { return this->map; }
        inline friend   std::ostream& operator<<(std::ostream& stream, const gaConfig& gac) { stream<<gac.to_str(); return stream; }
        inline const    std::string to_str      () const {
            using namespace std;
            stringstream ss;
            ss<<"{"<<endl;
#           define _str(x, y) ss<<"    "<<x<<": "<<y<<endl
            _str("Population",
                endl<<"    {"<<endl
                    <<"        Size: "<<this->population_size<<endl
                    <<"        Elitism Ratio: "<<this->elitism_ratio<<endl
                    <<"    }"
                );
            _str("Genes",
                endl<<"    {"<<endl
                    <<"        Count Boundary: "<<"[ "<<this->genes_count_min<<" , "<<this->genes_count_max<<" ]"<<endl
                    <<"        Value Boundary: "<<"[ "<<this->genes_value_min<<" , "<<this->genes_value_max<<" ]"<<endl
                    <<"    }"
                );
            _str("Combinations",
                 endl<<"    {"<<endl
                     <<"        Crossover probability: "<<this->crossover_prob<<endl
                     <<"        Mutation  probability: "<<this->mutation_prob<<endl
                     <<"    }"
                );
            _str("Generation Max Limit", this->gen_max_count);
            _str("Options",endl<<"    {");
            for(auto it = this->map.begin(); it!= this->map.end(); it++)
                ss<<"        "<<it->first<<":    "<<it->second<<endl;
            ss<<"    }"<<endl;
            ss<<"}";
#           undef _str
            return ss.str();
        }
    };
#undef VALIDATE_PROB
}
#endif // GACONFIG_H
