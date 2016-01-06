#ifndef NNET_H
#define NNET_H

#include "inc/stdafx.hpp"
#include "inc/nnconf.hpp"

#include <vector>

#include "Eigen/Core"

#define LOOP_UP(i, start, end) for(decltype(start + end) i = start; i < end; i++)
#define LOOP_DOWN(i, start, end) for(decltype(start + end) i = start; i > end; i--)
#define FOREACH_LAYER_ASC(i) LOOP_UP(i, 0, this->get_layer_no())
#define FOREACH_LAYER_DESC(i) LOOP_DOWN(i, this->get_layer_no(), 0)

typedef Eigen::MatrixXf layer;

class nnet
{
protected:
    /**
     * @brief vector_activation_func_t signature of the activation func. for vector input
     */
    typedef Eigen::VectorXf (*vector_activation_func_t)(const Eigen::VectorXf&);
protected:
    /**
     * @brief _size the network's size
     */
    vector<uint> _size;
    /**
     * @brief batch_size hold the batch size used to train the network
     */
    size_t _batch_size;
    /**
     * @brief _layers the layers' weight matrix of network
     */
    vector<layer> _layers;
    /**
     * @brief _layers the network's layers' momentum container
     */
    vector<layer> _momentum_layers;
    /**
     * @brief _layers the network's layers' cumulative delta weights
     */
    vector<layer> _delta_layers;
    /**
     * @brief _layers the network's biases' cumulative delta weights
     */
    vector<Eigen::VectorXf> _delta_biases;
    /**
     * @brief _biases the biases' weight matrix of network
     */
    vector<Eigen::VectorXf> _biases;
    /**
     * @brief _momentum_biases the network's biases' momentum container
     */
    vector<Eigen::VectorXf> _momentum_biases;
    /**
     * @brief _outputs neurons deltas
     */
    vector<Eigen::VectorXf> _deltas;
    /**
     * @brief _outputs neurons outputs
     */
    vector<Eigen::VectorXf> _outputs;
    /**
     * @brief _scalar_activations scalar input activation functions
     */
    vector<activation_func_t> _scalar_activations;
    /**
     * @brief _activations neuron layers' activation functions
     */
    function<Eigen::VectorXf (const Eigen::VectorXf&, size_t)>* _vector_activations;
    /**
     * @brief feedforward feedforward's network
     * @param in the input array
     * @param out the output matrix
     */
    void feedforward(const Eigen::VectorXf& in, Eigen::VectorXf* const out) const __nonnull();
public:
    /**
     * @brief nnet construct a neural network shell
     */
    nnet();
    /**
    * @brief ~nnet dtor the network
    */
    ~nnet();
    /**
     * @brief nnet construct a neural network
     * @param size the neural size
     */
    nnet(vector<uint> size);
    /**
     * @brief set_activations set activation functions for each layer
     * @param vaf the activation functions' vector
     */
    void set_activations(vector<activation_func_t> vaf);
    /**
     * @brief feedforward feedforward's the network
     * @param in the input array
     * @return the output of network
     */
    vector<Eigen::VectorXf> feedforward(const vector<Eigen::VectorXf*>& in) const;
    /**
     * @brief backprop train the network with backpropagation
     * @param inputs the input vectors
     * @param outputs the target vectors
     * @param learning_configs the network's learning configs
     */
    void backprop(const vector<Eigen::VectorXf*>& inputs, const vector<Eigen::VectorXf*>& target, learning_configs configs);
    /**
     * @brief get_MSE calc. the MSE value between output and target values
     * @param outs the output value
     * @param target the target value
     * @return the MSE value foreach output neurons
     */
    Eigen::VectorXf get_MSE(vector<Eigen::VectorXf> outs, const vector<Eigen::VectorXf*>& target) const;
    /**
     * @brief get_accuracy get accuracy of output VS. targets
     * @param outs the output value
     * @param target the target value
     * @return the accuracy value
     */
    scalar get_accuracy(vector<Eigen::VectorXf> outs, const vector<Eigen::VectorXf*>& target) const;
    /**
     * @brief get_layer_no get number of layers
     * @return the number of layers
     */
    inline size_t get_layer_no() const { return this->_size.size() - 1; }
};

#endif // NNET_H
