#ifndef MAIN_HELPER_HPP
#define MAIN_HELPER_HPP
#include "inc/stdafx.hpp"
#include <string>
#include <vector>
#include <limits>
#include <unistd.h>
#include <iostream>
#include <stdexcept>
using namespace std;

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/algorithm/string/predicate.hpp>
namespace po = boost::program_options;

#pragma GCC diagnostic ignored "-Wunused-parameter"
#   include "jsoncons/json.hpp"
using jsoncons::json;
#pragma GCC diagnostic pop

#include "inc/population.hpp"
using namespace TSPGA;

#define __PROMT_FAILED "[\u2613]"
#define __PROMT_PASSED "[\u221A]"
/**
 * @brief validate_arg Validates the passed argument
 * @param argc The arguments count
 * @param argv The arguments
 * @return The input argument varibale map
 */
po::variables_map validate_arg(const int argc, char** const argv) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help",                                                                "Produce help message")
        ("data-file,d",     po::value<string>()->default_value("data.dat"),     "The data input file")
        ("config-file,c",   po::value<string>()->default_value("config.json"),  "A json formatted config file")
        ("cache-file",      po::value<string>(),                                "The cache file to cache loaded data")
    ;
    // the argument variable map
    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch(const std::exception& e) {
        cerr<<e.what()<<endl;
        cerr<<endl<<desc<<endl;
        exit(EXIT_FAILURE);
    }
    // on help argument
    if (vm.count("help")) { cout << desc << endl; exit(EXIT_SUCCESS); }
    // return the variable map
    return vm;
}
/**
 * @brief load json configurations
 * @param config_file The configuration file
 * @return A json instance
 */
json load_configuration(const char* config_file) {
    json config;
    try{
         config = json::parse_file(config_file);
    } catch(const jsoncons::json_exception& e) {
        cerr<<e.what()<<endl;
        exit(EXIT_FAILURE);
    }
    return config;
}
/**
 * @brief Builds GA configurations based on the input json-configuration and data files
 * @param jconfig   The json configurations instance
 * @param _data     The input data
 * @return gaConfig A GA configuration instance
 */
gaConfig build_gaConfig(json jconfig, const TSPGA::data* const _data) {
    // temp. var for building gaconfig
    std::string var_build_gaConfig; bool var_build_gaConfig_match = false, var_build_gaConfig_break = false;
    // config path container
    vector<string> path;
    /**
     * These definitions are bounded to the scope of current function.
     * These are defined to facilitate the loading gaConfig stuffs.
     */
    #   define _switch(x)  var_build_gaConfig = x; var_build_gaConfig_match = var_build_gaConfig_break = false;
    #   define _case(x, o)    if(x == var_build_gaConfig) { var_build_gaConfig_match = true; path.push_back(x); o  if(path.size()) path.pop_back(); }
    #   define _break      var_build_gaConfig_break = true; if(path.size()) path.pop_back(); continue
    #   define _default()
    #   define _value(x) x->value()
    #   define _name(x)  x->name()
    #   define _IS_NA(x)    (_value(x).as<string>() == jconfig["meta"]["nonspecific_flag"].as<string>())
    #   define _iter(it, x) for (auto it = x.begin_members(); it != x.end_members(); ++it)
    #   define _iter_branch(it, b, section, x) _case(b, { _iter(section, _value(it)) { _switch(_name(section)) x } })
    #   define _iter_branch_final(it, b, section, x) _iter_branch(it, b, section, { _case(_name(section), x) })
    #   define _print(x) cout<<x->name()<<": "<<x->value()<<endl
    #   define _path(x)     boost::algorithm::join(path, ".")
    #   define _skip(x) cout<<__PROMT_FAILED<<" config."<<_path(x)<<endl; _break;
    #   define _max(t)  std::numeric_limits<t>::max()
    #   define _min(t)  std::numeric_limits<t>::min()
    #   define _panic(x) if(!var_build_gaConfig_match && !var_build_gaConfig_break) throw std::runtime_error("Unexpected configuration `" + _path(x) +"."+ _name(x) + ":" + _value(x).as_string() +"`");
    #   define _SHOULD_NOT_BE_NA(x)    if(_IS_NA(x)) _panic(x)
    #   define _add_as_option(e)   if(!_IS_NA(e)) gaconfig.option(_path(e), _value(e).as_string())
    /* TO SEE THE INPUT CONFIRGURATION UN-COMMENT BELOW */
    // cout<<jsoncons::pretty_print(jconfig)<<endl;
    // create an gaconfig instance (dummy initilization for now)
    gaConfig gaconfig(0, 0, 0);
    // start iteration of json-configuration elems
    _iter(it, jconfig) {
        path.clear();
        _switch(_name(it)) {
            // HEADER <combinations>
            _iter_branch(it, "combinations", c, {
                // SUB-HEADER <crossover>
                _iter_branch(c, "crossover", cc, {
                    _case("prob", {
                        gaconfig.setCrossover_prob(_value(cc).as_double());
                        _break;
                    });
                    _case("method", {
                        _add_as_option(cc);
                    });
                    _panic(cc);
                });
                // SUB-HEADER <mutation>
                _iter_branch(c, "mutation", cc, {
                    _case("prob", {
                        gaconfig.setMutation_prob(_value(cc).as_double());
                        _break;
                    });
                    _case("method", {
                        _add_as_option(cc);
                    });
                    _panic(cc);
                });
                _panic(c);
            });
            // HEADER <generations>
            _iter_branch(it, "generations", g, {
                // SUB-HEADER <count_limit>
                _case("count_limit", {
                    _switch(_value(g).as_string()) {
                        // VALUE-TYPE <MAX>
                        _case("MAX", {
                            gaconfig.setGenerationMaxCount(_max(size_t));
                            _break;
                        });
                        // VALUE-TYPE <DEFAULT>
                        _default() {
                            gaconfig.setGenerationMaxCount(boost::lexical_cast<size_t>(_value(g).as_string()));
                            _break;
                        };
                    }
                    _break;
                });
                // SUB-HEADER <storage_type>
                _case("storage_type", {
                    _add_as_option(g);
                    _break;
                });
                // SUB-HEADER <storage_detail>
                _iter_branch(g, "storage_details", s, {
                    _case("mem", {
                        _add_as_option(s);
                        _break;
                    });
                    _iter_branch_final(s, "disk", d, {
                        _add_as_option(d);
                    });
                  _panic(s);
                });
                _panic(g);
            });
            // HEADER <population>
            _iter_branch(it, "population", p, {
                // SUB-HEADER <size>
                _case("size", {
                    gaconfig.setPopulation_size(_value(p).as_ulonglong());
                    _break;
                });
                // SUB-HEADER <selection>
                _iter_branch(p, "selection", s, {
                    _case("method_type", {
                        _add_as_option(s);
                        _break;
                    });
                    _iter_branch(s, "method_details", d, {
                        _iter_branch_final(d, "tournament", m, {
                            _add_as_option(m);
                        });
                        _panic(d);
                    });
                    _panic(s);
                });
                // SUB-HEADER <elitism>
                _iter_branch(p, "elitism", e, {
                    _case("method", {
                        _add_as_option(e);
                        _break;
                    });
                    _case("ratio", {
                        gaconfig.setElitism_ratio(_value(e).as_double());
                        _break;
                    });
                    _panic(e);
                });
                _panic(p);
             });
            // HEADER <solution>
            _iter_branch(it, "solution", s, {
                // SUB-HEADER <genes>
                _iter_branch(s, "genes", g, {
                    // SUB-SUB-HEADER <count_boundary>
                    _case("count_boundary", {
                        _switch(_value(g).as_string()) {
                            // VALUE-TYPE <AUTO>
                            _case("AUTO", {
                                gaconfig.setGenes_CountLimit(_data->size(), _data->size());
                                _break;
                            });
                            // VALUE-TYPE <MAX>
                            _case("MAX", {
                                gaconfig.setGenes_CountLimit(0, _max(size_t));
                                _break;
                            });
                            // VALUE-TYPE <ARRAY : 1*2>
                            if(_value(g).is_array()) {
                                auto gv = _value(g).as_vector<int>();
                                if(gv.size() != 2)
                                    throw std::logic_error("Expecting to have [min, max] format for `"+_path(g)+"`");
                                gaconfig.setGenes_CountLimit(gv[0], gv[1]);
                                _break;
                            }
                            // VALUE-TYPE <ANY>
                            _panic(g);
                        }
                        _break;
                    });
                    // SUB-SUB-HEADER <value_boundary>
                    _case("value_boundary", {
                         _switch(_value(g).as_string()) {
                            // VALUE-TYPE <MAX>
                            _case("MAX", {
                                gaconfig.setGenes_ValueLimit(-_max(limit), _max(limit));
                                _break;
                            });
                            // VALUE-TYPE <ARRAY : 1*2>
                            if(_value(g).is_array()) {
                                auto gv = _value(g).as_vector<string>();
                                if(gv.size() != 2)
                                    throw std::logic_error("Expecting to have [min, max] format for `"+_path(g)+"`");
                                if(boost::to_upper_copy(gv[1]) == "AUTO")
                                    gv[1] = boost::lexical_cast<std::string>(_data->size() - 1);
                                gaconfig.setGenes_ValueLimit(atof(gv[0].c_str()), atof(gv[1].c_str()));
                                _break;
                            }
                            // VALUE-TYPE <ANY>
                            _panic(g);
                        }
                        _break;
                    });
                });
                // SUB-SUB-HEADER <fitness>
                _iter_branch(s, "fitness", f, {
                    _case("method", {
                        _add_as_option(f);
                        _break;
                    });
                    _panic(f);
                });
             });
        }
    }
    // return the configuration
    return gaconfig;
    #   undef _switch
    #   undef _case
    #   undef _break
    #   undef _default
    #   undef _value
    #   undef _name
    #   undef _IS_NA
    #   undef _iter
    #   undef _iter_branch
    #   undef _iter_branch_final
    #   undef _print
    #   undef _skip
    #   undef _max
    #   undef _min
    #   undef _path
    #   undef _panic
    #   undef _SHOULD_NOT_BE_NA
    #   undef _add_as_option
}
/**
 * @brief loads data from file
 * @param data_file The target data file
 * @return A data instance
 */
TSPGA::data* load_data(const char* data_file) {
#   define _read_coord(ss, s, x) ss.clear (); ss << s; ss >> s; ss >> x.longitude >> x.latitude;
    // open the data file
    std::ifstream in(data_file);
    // check if it is open?
    if(!in) throw std::runtime_error("Couldn't open data file");
    // define variables
    stringstream ss;
    string line;
    TSPGA::coordinate location;
    coordinates* locations = new coordinates();
    TSPGA::data* _data = new TSPGA::data();
    // while reading
    while(in) {
        getline (in, line,'\n');
        if(boost::starts_with(line, "EDGE_WEIGHT_TYPE")) {
            const boost::regex e("EDGE_WEIGHT_TYPE\\s*:\\s*([A-Za-z0-9_]+).*");
            boost::match_results<std::string::const_iterator> matches;
            if (boost::regex_match(line, matches, e))
                _data->type(string(matches[1].first, matches[1].second).c_str());
            else
                throw std::runtime_error("Unrecognized `EDGE_WEIGHT_TYPE` format");
        }
        // prune the header text
        if(!isdigit(line[0])) continue;
        // using the stringstream read coordinate of current line
        // and store it to location instance
        _read_coord(ss, line, location);
        // add the location to crowd
        locations->push_back(location);
        printf("\r%lu coordinates has been read...", locations->size());
    }
    printf("\r%s %lu coordinates has been read.\n", __PROMT_PASSED, locations->size());
    // set the fetched locations
    _data->set(locations);
    // construct and return data
    return _data;
#   undef _read_coord
}
/**
 * @brief load data considering the cache file
 * @param vm The variable map from boost option lib.
 */
TSPGA::data* load_data_cached(const po::variables_map vm) {
    TSPGA::data* _data = NULL;
    auto use_cache  = vm.count("cache-file");
    auto cache_file = use_cache ? vm["cache-file"].as<string>().c_str() : "";
    auto data_file  = vm["data-file"].as<string>().c_str();
    // validate the existance of cache file
    if(!vm.count("cache-file") || !boost::filesystem::exists(cache_file)) {
        // load data
        _data = load_data(data_file);
        // if any cache file has been defined?
        if(use_cache) {
            cout<<"[.] Serializing retrieved data into `"<<cache_file<<"`\r"<<std::flush;
            // serialize data
            data::serialize(_data, cache_file);
            cout<<__PROMT_PASSED<<" Serialization stored "<<_data->size()<<" objects into `"<<cache_file<<'`'<<endl<<std::flush;
        }
    } else {
        if(!use_cache)
            throw std::runtime_error("Expecting to `--cache-file` be set.");
        cout<<"[.] Deserializing stored data from `"<<cache_file<<"`\r"<<std::flush;
        // deserialize data
        _data = data::deserialize(cache_file);
        cout<<__PROMT_PASSED<<" Deserlization retrieved "<<_data->size()<<" objects from `"<<cache_file<<'`'<<endl<<std::flush;
    }
    return _data;
}

#endif // MAIN_HELPER_HPP
