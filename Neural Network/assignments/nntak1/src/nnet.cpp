#include "inc/nnet.hpp"

#include <cfloat>
#include <assert.h>
#include <algorithm>

#include <iostream>

nnet::nnet() { }

nnet::nnet(vector<uint> size)
    : _size(size), _batch_size(0)
{
    // alloc the layers' deltas
    this->_deltas = vector<Eigen::VectorXf>(size.size());
    // alloc the activations
    this->_vector_activations = new function<Eigen::VectorXf (const Eigen::VectorXf&, size_t)>[size.size()];
    // set activations to null
    std::fill_n(this->_vector_activations, size.size(), nullptr);
    // alloc the layers' outputs/biases & init with zero
    LOOP_UP(i, 0, size.size()) {
        // create a container for neurons outputs
        this->_outputs.push_back(Eigen::VectorXf(size[i]));
        // create a container for delta biases
        this->_delta_biases.push_back(Eigen::VectorXf(size[i]));
        // create a random valued biases weights
        this->_biases.push_back(Eigen::VectorXf::Random(size[i]));
        // create a zero valued momentum value
        this->_momentum_biases.push_back(Eigen::VectorXf::Zero(size[i]));
    }
	// a fail-safe, set input's layer's biases to zero[we won't use them!!]
    this->_biases.front().setZero();
    // foreach layer
    FOREACH_LAYER_ASC(i) {
        // create a randomly initialized layer
        this->_layers.push_back(layer::Random(size[i + 1], size[i]));
        // assign the zero values for delta weights in each layer
        this->_delta_layers.push_back(layer::Zero(size[i + 1], size[i]));
        // create a zero valued momentum value
        this->_momentum_layers.push_back(layer::Zero(size[i + 1], size[i]));
    }
    // validate the layers number
    assert(this->_layers.size() == this->get_layer_no());
}

nnet::~nnet() {
    // free-up every bit we have allocated in the ctor
    delete[] this->_vector_activations;
}

void nnet::set_activations(vector<activation_func_t> vaf) {
    // if only one activation func. passed, assuming that is for every neuron layers
    // or the input is omiting the input neurons' activation func. expand the input vector,
    //  the input neurons' acitvation func. will be altered later in this function
    if(vaf.size() == 1) while(vaf.size() != this->_size.size()) vaf.push_back(vaf.back());
    // make place for input neurons' activation functions
    if(vaf.size() == this->_size.size() - 1) vaf.insert(vaf.begin(), nullptr);
    // the size of activation func. vector should match with the neuron layers' size
    if(vaf.size() != this->_size.size()) e_raise("The activations' size should match the network's size!");
    // store the scalar activation functions
    this->_scalar_activations = vaf;
    // copy the activation call size
    LOOP_UP(i, 1/* exclude input neurons*/, this->_size.size()) {
        this->_vector_activations[i] = [&](const Eigen::VectorXf& v, size_t layer) {
            Eigen::VectorXf o(v);
            for(int i = 0; i < o.rows(); i++) o(i) = this->_scalar_activations[layer - 1](o(i));
            return o;
        };
    }
    // assign identical activation func. for input neurons
    this->_vector_activations[0] = [](const Eigen::VectorXf& v, size_t) { return v; };
}

void nnet::backprop(const vector<Eigen::VectorXf*>& inputs, const vector<Eigen::VectorXf*>& target, learning_configs configs) {
    // make sure we have paired input and target values correctly
    assert(inputs.size() == target.size());
    // 1 foreach input
    //  1.1 compute output errors
    //  1.2 foreach hidden layer, backpropagate the errors
    //      1.2.1 store the errors in `::_delta_layers` matrix for layer's weights and `::_delta_biases` for layer's biases
    // 2 update the weights with deltas stored in `::_delta_{biases|layers}` matrix
    LOOP_UP(index, 0, inputs.size()) {
        Eigen::VectorXf
            *x = inputs[index],
            *t = target[index];
        // input and output match with configuration
        assert(t->rows() == this->_size.back() && x->rows() == this->_size.front());
        // feedforward the network and capture the neurons' outputs
        this->feedforward(*x, this->_outputs.data());
        // compute the average error output's deltas
        this->_deltas.back() = *t - this->_outputs.back();
        // foreach layer
        FOREACH_LAYER_DESC(i) {
            // update the biases' deltas
            this->_delta_biases[i] += this->_deltas[i];
            // update the weights' deltas
            this->_delta_layers[i - 1] += this->_deltas[i] * this->_outputs[i - 1].transpose();
            // it don't matter if we compute delta for input neurons, because we will never loop into there >:)
            Eigen::MatrixXf
                // o * (1 - o)
                d_o = this->_outputs[i - 1].cwiseProduct(Eigen::MatrixXf(1 - this->_outputs[i - 1].array())),
                // W * delta
                d_w = this->_layers[i - 1].transpose() * this->_deltas[i];
            // o * (1 - o) .* (W * delta)
            this->_deltas[i - 1] = d_o.cwiseProduct(d_w);
        }
    }
    // update the weights
    FOREACH_LAYER_ASC(i) {
        // calc the delta weights
        layer delta_layer = configs._learning_rate * this->_delta_layers[i] / inputs.size();
        Eigen::VectorXf delta_bias = configs._learning_rate * this->_delta_biases[i] / inputs.size();
        // a work-around for times which the {layers|bias} as so small which change to NaN!
        // we replace them with a minimum value
        replace_nan(delta_bias, FLT_MIN);
        replace_nan(delta_layer, FLT_MIN);
        // update the weights : { weights = (1 - gamma) * weights + delta_weights - momentums}
        this->_biases[i] =
            this->_biases[i] +
            delta_bias  -
            configs._momentum_rate * this->_momentum_biases[i];
        this->_layers[i] =
            (1 - configs._learning_rate * configs._regularization_rate / inputs.size()) * this->_layers[i] +
            delta_layer -
            configs._momentum_rate * this->_momentum_layers[i];
        // store the momentum values
        this->_momentum_biases[i] = delta_bias;
        this->_momentum_layers[i] = delta_layer;
    }
    // clean up the network's buffers values
    for(auto& d : this->_deltas)        d.setZero();
    for(auto& d : this->_delta_layers)  d.setZero();
    for(auto& d : this->_delta_biases)  d.setZero();
}

Eigen::VectorXf nnet::get_MSE(vector<Eigen::VectorXf> outs, const vector<Eigen::VectorXf*>& target) const {
    // assert the sizes
    assert(outs.size() == target.size());
    // compute the (t - o)^2 for each input's output
    LOOP_UP(i, 0, outs.size()) outs[i] = pow(*(target[i]) - outs[i], 2);
    // create the error vector
    Eigen::VectorXf e(this->_size.back(), 1); e.setZero();
    // foreach input
    LOOP_UP(i, 0, outs.size()) {
        // foreach output neuron
        LOOP_UP(j, 0, outs[i].rows())
            // calc. the average over each column
            e(j) += outs[i](j) / outs.size();
    }
    // return the mse
    return e;
}

scalar nnet::get_accuracy(vector<Eigen::VectorXf> outs, const vector<Eigen::VectorXf*> &target) const {
    // assert the sizes
    assert(outs.size() == target.size());
    // return the pair of max_value and max_index
    auto max = [](const Eigen::VectorXf& v) -> pair<scalar, size_t> {
        if(!v.rows()) e_raise("Vector cannot be empty!");
        pair<scalar, size_t> p = { v[0], 0 };
        LOOP_UP(i, 1, v.rows()) if(p.first < v[i]) p = { v[i], i };
        return p;
    };
    scalar accuracy = 0;
    // count the corrects
    LOOP_UP(i, 0, outs.size()) if(max(outs[i]).second == max(*target[i]).second) accuracy++;
    // compute the accuracy
    return accuracy / outs.size();
}

vector<Eigen::VectorXf> nnet::feedforward(const vector<Eigen::VectorXf*>& ins) const {
    // alloc a virtual output matrix
    Eigen::VectorXf* o = new Eigen::VectorXf[this->_size.size()]; LOOP_UP(i, 0, this->_size.size()) o[i] = Eigen::VectorXf(this->_size[i]);
    // the outputs vector
    vector<Eigen::VectorXf> out;
    // foreach inputs
    for(auto in : ins) {
        // the input cannot be null
        assert(in);
        // feed forward the network with input array and output matrix
        this->feedforward(*in, o);
        // fetch the output vector for corresponding input
        out.push_back(o[this->_size.size() - 1]);
    }
    // free-up the allocated memory
    delete[] o;
    // returns the output vector
    return out;
}

void nnet::feedforward(const Eigen::VectorXf& in, Eigen::VectorXf* const out) const {
    // input's size should match with network's configuration
    assert(in.rows() == this->_size.front());
    // fetch the output value of input layer
    //  note the activation function for input layer is identical function,
    //  and this is for just code consistency.
    out[0] = this->_vector_activations[0](in, 0);
    LOOP_UP(i, 1/*exclude input layer*/, this->_size.size()) {
        // the layers should match with network's configuration
        assert(out[i].rows() == this->_size[i]);
        // current_layer's_output = activation(previous_layer_weights * previous_layer_output_as current's_input + current_layer's_biases)
        out[i] = this->_vector_activations[i](this->_layers[i - 1] * out[i - 1] + this->_biases[i], i);
    }
}

