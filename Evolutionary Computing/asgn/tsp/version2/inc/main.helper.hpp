#ifndef MAIN_HELPER_HPP
#define MAIN_HELPER_HPP

#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#pragma GCC diagnostic ignored "-Wunused-parameter"
#   include "jsoncons/json.hpp"
using jsoncons::json;
#pragma GCC diagnostic pop

#define ELIT_SELECTION_RATIO 0.2

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
        ("input-file,i",    po::value<string>()->default_value("in.dat"),       "The data input file")
        ("output-file,o",   po::value<string>()->default_value("out.dat"),      "The result output file")
#ifndef XFILE_COMPATIBLE
        ("config-file,c",   po::value<string>()->default_value("config.json"),  "A json formatted config file")
#endif
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
    // validate the input file
    if (!boost::filesystem::exists(vm["input-file"].as<string>()))
    { cerr<<"input file `"<<vm["input-file"].as<string>()<<"` does not exist!"<<endl; exit(EXIT_FAILURE); }
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
#ifndef XFILE_COMPATIBLE
         config = json::parse_file(config_file);
#else
        ifstream infile(config_file);
        if(!infile.is_open()) throw jsoncons::json_exception_0("Unable to open input file");
        vector<string> settings = {
            "population_size", "i",
            "chrom_size", "i",
            "max_generations", "i",
            "crossover_prob", "d",
            "mutation_prob", "d",
            "output-file", "s"
        };
        for(size_t line = 0; line < settings.size(); line += 2) {
            switch(settings[line+1][0]) {
                case 'i': {
                    int set;
                    infile >> set;
                    config[settings[line]] = set;
                    break;
                }
                case 'd': {
                    double set;
                    infile >> set;
                    config[settings[line]] = set;
                    break;
                }
                case 's': {
                    string set;
                    infile >> set;
                    config[settings[line]] = set;
                    break;
                }
                default:
                    throw jsoncons::json_exception_0("Invalid data type");
            }
        }
        config["elit_prop"] = ELIT_SELECTION_RATIO;
#endif
    } catch(const jsoncons::json_exception& e) {
        cerr<<e.what()<<endl;
        exit(EXIT_FAILURE);
    }
    return config;
}

#endif // MAIN_HELPER_HPP
