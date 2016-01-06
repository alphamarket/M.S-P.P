#include "inc/stdafx.hpp"

#include <string>
#include <thread>
#include <fstream>
#include <signal.h>
#include <iostream>
#include <boost/filesystem.hpp>

#include "inc/nnet.hpp"
#include "inc/load_mat.hpp"
#include "inc/main.helper.hpp"


#define _start_of(s) cout<<PROMT_WAIT<<s;
#define _end_of(s) cout<<"\r"<<PROMT_PASSED<<s<<endl;

#define SECTION_BATCH_LEARNING      1
#define SECTION_MINI_BATCH_LEARNING 11
#define SECTION_ONLINE_LEARNING     2
#define SECTION_MOMENTUM            3
#define SECTION_REGULARIZATION      4

#define CONF_MAX_EPOC                   30
#define CONF_VALIDATION_SIZE            1e+4
#define CONF_HIDDEN_LAYERS              100

#define CONF_BATCH_LEARNING
#define CONF_LEARNING_RATE(x)           0.4
#define CONF_MOMENTUM_RATE(x)           0
#define CONF_REGULARIZATION_RATE(x)     0
#define CONF_MAX_ACCURACY_FAILURE_TRIAL 5

#define CONF_IMAGE_FILE                 string("images.mat")
#define CONF_LABEL_FILE                 string("labels.mat")

#define SECTION SECTION_MINI_BATCH_LEARNING

#if SECTION == SECTION_BATCH_LEARNING
#   undef   SECTION_MINI_BATCH_LEARNING
#elif SECTION == SECTION_ONLINE_LEARNING
#   undef   CONF_BATCH_LEARNING
#   undef   CONF_LEARNING_RATE
#   define  CONF_LEARNING_RATE(x)       0.3
#elif SECTION == SECTION_MOMENTUM
#   undef   CONF_MOMENTUM_RATE
#   define  CONF_MOMENTUM_RATE(x)       0.2
#elif SECTION == SECTION_REGULARIZATION
#   undef   CONF_REGULARIZATION_RATE
#   define  CONF_REGULARIZATION_RATE(x) 0.2
#elif SECTION == SECTION_MINI_BATCH_LEARNING
#   undef CONF_MINI_BATCH_SIZE
#   define CONF_MINI_BATCH_SIZE         10
#   undef CONF_IMAGE_FILE
#   define CONF_IMAGE_FILE              string("images_pca.mat")
#else
#   error UNDEFINED SECTION
#endif

#define INITIAL_LEARNING_RATE           CONF_LEARNING_RATE(0)
#define INITIAL_MOMENTUM_RATE           CONF_MOMENTUM_RATE(0)
#define INITIAL_REGULARIZATION_RATE     CONF_REGULARIZATION_RATE(0)

bool training = false;
bool stop_trainging = false;

string get_log_file_name() {
    string o = "";
#ifdef CONF_BATCH_LEARNING
#   ifdef CONF_MINI_BATCH_SIZE
    o += "MINI-";
#   endif
    o += "BATCH";
#else
    o += "ONLINE";
#endif
    o += "-LR-" + to_string(INITIAL_LEARNING_RATE).substr(0, 4);
    o += "-MR-" + to_string(INITIAL_MOMENTUM_RATE).substr(0, 4);
    o += "-RR-" + to_string(INITIAL_REGULARIZATION_RATE).substr(0, 4);
    o += "-DATE-" + get_date_time();
    return o;
}

scalar mean(const Eigen::VectorXf& v) {
    scalar o;
    LOOP_UP(i, 0, v.rows()) o += v(i);
    return o / v.rows();
}

int main (int, char**) {
#ifdef DEBUG
    chdir(string(const_cast<char*>(__FILE__)).replace(string(const_cast<char*>(__FILE__)).find(basename(__FILE__)), string(basename(__FILE__)).length(), "").c_str());
#endif
    std::cout.setf(std::ios_base::unitbuf);
    bool stop = false;
    // seed updator
    std::thread([&stop](){ while(!stop) { updateseed(); usleep(25); } }).detach();
    signal(SIGQUIT, [](int) { cout<<"\n[ABORT]"<<endl; exit(EXIT_FAILURE); });
    signal(SIGINT, [](int) {
        static int count = 0;
        if(!training || ++count > 1) raise(SIGQUIT);
        cout<<"\n[ Wating for learner to finish current epoch, PRESS CTRL+C AGAIN FOR HARD-KILL! ]"<<endl;
        stop_trainging = true;
    });
    // configs
    const size_t
        max_epoc = CONF_MAX_EPOC,
        size_of_validation = CONF_VALIDATION_SIZE;
    // loading trainset
    _start_of("Loading database")
        vector<Eigen::VectorXf*> images, labels, train_images, train_labels, test_images, test_labels, eval_images, eval_labels;
        images = load_mnist_mat("data/training/" + CONF_IMAGE_FILE, "images");
        labels = load_mnist_mat("data/training/" + CONF_LABEL_FILE, "labels");
        test_images = load_mnist_mat("data/test/" + CONF_IMAGE_FILE, "images");
        test_labels = load_mnist_mat("data/test/" + CONF_LABEL_FILE, "labels");
        train_images = vector<Eigen::VectorXf*>(vector<Eigen::VectorXf*>(images.begin(), images.begin() +  images.size() - size_of_validation));
        train_labels = vector<Eigen::VectorXf*>(vector<Eigen::VectorXf*>(labels.begin(), labels.begin() +  images.size() - size_of_validation));
        eval_images = vector<Eigen::VectorXf*>(vector<Eigen::VectorXf*>(images.end() - size_of_validation, images.end())),
        eval_labels = vector<Eigen::VectorXf*>(labels.end() - size_of_validation, labels.end());
        assert(train_images.size() == train_labels.size());
        assert(eval_images.size() == eval_labels.size());
        assert(test_images.size() == test_labels.size());
        images.clear(), labels.clear();
    _end_of("Loading database");
    // creating network
    _start_of("Creating a network");
        vector<uint> net_config = {(uint)images.front()->rows(), CONF_HIDDEN_LAYERS, (uint)labels.front()->rows()};
        nnet nn(net_config);
        nn.set_activations({[](scalar s) { return 1 / (1 + exp(-s)); }});
    _end_of("A network with configuration `"<<net_config<<"` created.");
    // logs
    const string log_dir = "logs";
    if(!boost::filesystem::exists(log_dir)) boost::filesystem::create_directory(log_dir);
    string log_file = log_dir + "/" + get_log_file_name();
    ofstream logs(log_file + ".log");
    logs.unsetf(ios_base::unitbuf);
    if(!logs.is_open()) throw runtime_error("Cannot open log file!");
    // promt
    cout
        <<"\n--------------------------\n"<<endl
        <<"Maximum epoc#: "<<max_epoc<<endl
        <<"Training set's size: "<<train_images.size()<<endl
        <<"Evaluation set's size: "<<eval_images.size()<<endl
        <<"Test set's size: "<<test_images.size()<<endl
        <<"Learning mode: "
#ifdef CONF_BATCH_LEARNING
#   ifdef CONF_MINI_BATCH_SIZE
        <<"MINI-BATCH"<<endl
        <<"Batch Size: "<<CONF_MINI_BATCH_SIZE<<endl
#   else
        <<"BATCH"<<endl
#   endif
#else
        <<"ONLINE"<<endl
#endif
        <<"Initial Learning rate: "<<INITIAL_LEARNING_RATE<<endl
        <<"Initial Momentum rate: "<<INITIAL_MOMENTUM_RATE<<endl
        <<"Initial Regularization rate: "<<INITIAL_REGULARIZATION_RATE<<endl
        <<"Maximum accuracy failure trail: "<<CONF_MAX_ACCURACY_FAILURE_TRIAL<<endl
        <<"\n--------------------------"<<endl;
    // previous accuracy
    scalar prev_accuracy = 0;
    // number of accuracy failure
    size_t accuracy_failure  = 0, accuracy_hold = 0;
    // learning configs
    learning_configs configs = { INITIAL_LEARNING_RATE, INITIAL_MOMENTUM_RATE, INITIAL_REGULARIZATION_RATE };
    // foreach epoc
    LOOP_UP(epoc, 0, max_epoc) {
        training = true;
        if(stop_trainging) break;
        cout<<endl<<configs<<endl;
        // training
        _start_of("Training epoc# "<<epoc)
#ifndef CONF_BATCH_LEARNING
            LOOP_UP(j, 0, train_images.size())
                nn.backprop({train_images[j]}, {train_labels[j]}, configs);
#else
#   ifndef CONF_MINI_BATCH_SIZE
            nn.backprop(train_images, train_labels, configs);
#   else
            LOOP_UP(j, 0, train_images.size() / CONF_MINI_BATCH_SIZE) {
                vector<Eigen::VectorXf*>
                    mini_batch_images(train_images.begin() + j * CONF_MINI_BATCH_SIZE, train_images.begin() + (j + 1) * CONF_MINI_BATCH_SIZE),
                    mini_batch_labels(train_labels.begin() + j * CONF_MINI_BATCH_SIZE, train_labels.begin() + (j + 1) * CONF_MINI_BATCH_SIZE);
                nn.backprop(mini_batch_images, mini_batch_labels, configs);
            }
#   endif
#endif
        _end_of("Epoc# "<<epoc<<" trained   ");
        // evaluating
        _start_of("Evaluating the epoc# "<<epoc)
            /**
             * acc./mse over trainset
             */
            vector<Eigen::VectorXf> results = nn.feedforward(train_images);
            Eigen::VectorXf train_mse = nn.get_MSE(results, train_labels);
            scalar accuracy_train = nn.get_accuracy(results, train_labels);
            /**
             * acc./mse over evalset
             */
            results = nn.feedforward(eval_images);
            Eigen::VectorXf eval_mse = nn.get_MSE(results, eval_labels);
            scalar accuracy_eval = nn.get_accuracy(results, eval_labels);
            /**
             * acc./mse over testset
             */
             results = nn.feedforward(test_images);
             Eigen::VectorXf test_mse = nn.get_MSE(results, test_labels);
             scalar accuracy_test = nn.get_accuracy(results, test_labels);
        _end_of("Evaluating the epoc# "<<epoc);
        scalar train_mse_mean = mean(train_mse), eval_mse_mean = mean(eval_mse), test_mse_mean = mean(test_mse);
        logs<<train_mse_mean<<" "<<eval_mse_mean<<" "<<test_mse_mean<<" "<<accuracy_train<<" "<<accuracy_eval<<" "<<accuracy_test<<endl;
        cout
            <<endl
            <<"mean(Train MSE): "<<train_mse_mean<<endl
            <<"mean(Eval. MSE): "<<eval_mse_mean<<endl
            <<"mean(Test  MSE): "<<test_mse_mean<<endl
            <<endl
            <<"Train Accuracy : "<<accuracy_train<<endl
            <<"Eval. Accuracy : "<<accuracy_eval<<endl
            <<"Test  Accuracy : "<<accuracy_test<<endl
            <<endl
            <<"--------------------------------"<<endl;
        if(eval_mse_mean < 1e-2) { cout<<"\nTraining ended at "<<epoc<<"iteration with MSE mean of: "<<eval_mse_mean<<endl; break; }
        if(accuracy_eval < prev_accuracy) {
            if(++accuracy_failure > CONF_MAX_ACCURACY_FAILURE_TRIAL) {
                cout<<"\n\nTraining stopped due to accuracy failure overflow("<<accuracy_failure<<" times)."<<endl;
                stop_trainging = true;
            }
            cout<<"Accuracy failure#"<<accuracy_failure<<" occurred."<<endl;
        }
        else if(accuracy_eval > prev_accuracy) accuracy_failure = 0;
        if(accuracy_eval == prev_accuracy) accuracy_hold++;
        else accuracy_hold = 0;
        prev_accuracy = accuracy_eval;
    }
    training = false;
    logs.close();
    _start_of("Releasing resouces");
        stop = true;
        sleep(1);
    _end_of("Releasing resouces");
    return EXIT_SUCCESS;
}
