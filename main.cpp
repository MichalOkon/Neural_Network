#include <cmath>
#include <initializer_list>
#include <iostream>
#include <list>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <cfloat>


template<typename T>
class Matrix {
private:

    std::vector<T> data;
    int rows;
    int cols;

public:

    template<typename U> friend
    class Matrix;

    Matrix() {};

    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        data = std::vector<T>(rows * cols, 0);
    }

    Matrix(int rows, int cols, const std::initializer_list<T> &list) : rows(rows), cols(cols) {
        if (rows * cols != (int) list.size()) {
            throw std::length_error(
                    "The number of columns and rows does not match the given list in the matrix construction");
        }
        // Check if that actually works
        data = std::vector(list);
    }

    Matrix(Matrix &other) : data(other.data), rows(other.rows), cols(other.cols) {}

    Matrix(Matrix &&other) noexcept: data(std::move(other.data)), rows(other.rows), cols(other.cols) {
        std::destroy(other.data.begin(), other.data.end());
        other.rows = 0;
        other.cols = 0;
    }

    ~Matrix() = default;

    // How to test those two?
    // Copy assignment operator
    Matrix &operator=(const Matrix &other) {
        if (this != &other) {
            data = other.data;
            rows = other.rows;
            cols = other.cols;
        }
        return *this;
    }

    // Move assignment operator
    Matrix &operator=(Matrix &&other) noexcept {
        if (this != &other) {
            data = other.data;
            rows = other.rows;
            cols = other.cols;
            other.rows = 0;
            other.cols = 0;
        }
        return *this;
    }

    T &operator[](const std::pair<int, int> &ij) {
        if (ij.first >= rows) {
            throw std::length_error("Number of rows exceeded");
        }
        if (ij.second >= cols) {
            throw std::length_error("Number of columns exceeded");
        }
        return data[ij.first * cols + ij.second];
    }

    // When is this one used???
    const T &operator[](const std::pair<int, int> &ij) const {
        if (ij.first >= rows) {
            throw std::length_error("Number of rows exceeded");
        }
        if (ij.second >= cols) {
            throw std::length_error("Number of columns exceeded");
        }
        return data[ij.first * cols + ij.second];
    }

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator*(U x) const {
        Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix(this->rows, this->cols);
        for (int i = 0; i < (int) this->data.size(); i++) {
            newMatrix.data[i] = this->data[i] * x;
        }
        return newMatrix;
    }

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator*(const Matrix<U> &B) const {
        if (this->cols != B.rows) {
            throw std::length_error("Matrices' sizes incompatible for multiplication");
        }
        Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix(this->rows, B.cols);
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < B.cols; j++) {
                typename std::common_type<T, U>::type sum = 0;
                for (int k = 0; k < this->cols; k++) {
                    sum += this->data[i * cols + k] * B.data[k * B.cols + j];
                }
                newMatrix.data[i * B.cols + j] = sum;
            }
        }
        return newMatrix;
    }

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator+(const Matrix<U> &B) const {
        if (B.rows != 1 && rows != 1 && (B.rows != this->rows || B.cols != this->cols)) {
            throw std::length_error("Matrices' sizes incompatible for addition");
        }
        Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix<typename std::common_type<T, U>::type>(
                std::max(this->rows, B.rows), std::max(this->cols, B.cols));
        if (B.rows == 1) {
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->cols; j++) {
                    newMatrix.data[i * cols + j] = this->data[i * cols + j] + B.data[j];
                }
            }
        } else if (rows == 1) {
            for (int i = 0; i < B.rows; i++) {
                for (int j = 0; j < B.cols; j++) {
                    newMatrix.data[i * cols + j] = B.data[i * cols + j] + this->data[j];
                }
            }
        } else {
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->cols; j++) {
                    newMatrix.data[i * cols + j] = this->data[i * cols + j] + B.data[i * cols + j];
                }
            }
        }
        return newMatrix;
    }

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator-(const Matrix<U> &B) const {
        if (B.rows != 1 && (B.rows != this->rows || B.cols != this->cols)) {
            throw std::length_error("Matrices' sizes incompatible for subtraction");
        }
        Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix<typename std::common_type<T, U>::type>(
                this->rows, this->cols);
        if (B.rows == 1) {
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->cols; j++) {
                    newMatrix.data[i * cols + j] = this->data[i * cols + j] - B.data[j];
                }
            }
        } else {
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->cols; j++) {
                    newMatrix.data[i * cols + j] = this->data[i * cols + j] - B.data[i * cols + j];
                }
            }
        }
        return newMatrix;
    }

    Matrix transpose() const {
        Matrix<T> transposedMatrix(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposedMatrix.data[j * rows + i] = this->data[i * cols + j];
            }
        }

        return transposedMatrix;
    }

    int getRows() const {
        return rows;
    }

    int getCols() const {
        return cols;
    }

    // For debugging
    void print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << this->data[i * cols + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::string toString() const {
        std::string result = "[";
        for (int i = 0; i < rows; i++) {
            result = result.append("[");
            for (int j = 0; j < cols; j++) {
                result = result.append(std::to_string(this->data[i * cols + j]));
                if (j != cols - 1) {
                    result = result.append(", ");
                }
            }
            result = result.append("]");
            if (i != rows - 1) {
                result = result.append(", ");
            }
        }
        result = result.append("]");
        return result;
    }
};

template<typename T>
class Layer {

    virtual Matrix<T> forward(const Matrix<T> &x) = 0;

    virtual Matrix<T> backward(const Matrix<T> &dy) = 0;

};

template<typename T>
class Linear : public Layer<T> {
private:
    Matrix<T> cache;
    Matrix<T> bias_grads;
    Matrix<T> weight_grads;
    Matrix<T> weights;
    Matrix<T> bias;
    int n_samples;
    int in_features;
    int out_features;

public:

    Linear(const int in_features,
           const int out_features,
           const int n_samples,
           const int seed) : n_samples(n_samples), in_features(in_features), out_features(out_features) {

        std::default_random_engine generator(seed);
        std::normal_distribution<T> distribution_normal(0.0, 1.0);
        std::uniform_real_distribution<T> distribution_uniform(0.0, 1.0);

        // create matrices
        this->weights = Matrix<T>(in_features, out_features);
        this->weight_grads = Matrix<T>(in_features, out_features);
        this->bias = Matrix<T>(1, out_features);
        this->bias_grads = Matrix<T>(1, out_features);
        this->cache = Matrix<T>(n_samples, in_features);

        // set initial weights and weight gradients
        for (int i = 0; i < in_features; i++) {
            for (int j = 0; j < out_features; j++) {
                this->weights[{i, j}] = distribution_normal(generator);
                this->weight_grads[{i, j}] = 0;
            }
        }

        // set initial bias and bias grads
        for (int i = 0; i < out_features; i++) {
            this->bias[{0, i}] = distribution_uniform(generator);
            this->bias_grads[{0, i}] = 0;
        }

        // set initial cache values
        for (int i = 0; i < n_samples; i++) {
            for (int j = 0; j < in_features; j++) {
                this->cache[{i, j}] = 0;
            }
        }

    }

    ~Linear() = default;

    virtual Matrix<T> forward(const Matrix<T> &x) override final {

        Matrix<T> out(1, this->out_features);

        for (int i = 0; i < this->out_features; i++) {

            T activation = 0;

            for (int w_i = 0; w_i < this->in_features; w_i++) {
                activation += this->weights[{w_i, i}] * x[{0,w_i}];
            }

            activation += this->bias[{0,i}];

            out[{0,i}] = activation;
        }

        return out;
    }

    virtual Matrix<T> backward(const Matrix<T> &dy) override final {
        Matrix<T> out(1, this->in_features);
        return out;
    }

};

template<typename T>
class ReLU : public Layer<T> {

private:
    Matrix<T> cache;
    Matrix<T> bias_grads;
    Matrix<T> weight_grads;
    Matrix<T> weights;
    Matrix<T> bias;
    int n_samples;
    int in_features;
    int out_features;

public:

    ReLU(const int in_features,
         const int out_features,
         const int n_samples,
         const int seed) : n_samples(n_samples), in_features(in_features), out_features(out_features) {

        std::default_random_engine generator(seed);
        std::normal_distribution<T> distribution_normal(0.0, 1.0);
        std::uniform_real_distribution<T> distribution_uniform(0.0, 1.0);

        // create matrices
        this->weights = Matrix<T>(in_features, out_features);
        this->weight_grads = Matrix<T>(in_features, out_features);
        this->bias = Matrix<T>(1, out_features);
        this->bias_grads = Matrix<T>(1, out_features);
        this->cache = Matrix<T>(n_samples, in_features);

        // set initial weights and weight gradients
        for (int i = 0; i < in_features; i++) {
            for (int j = 0; j < out_features; j++) {
                this->weights[{i, j}] = distribution_normal(generator);
                this->weight_grads[{i, j}] = 0;
            }
        }

        // set initial bias and bias grads
        for (int i = 0; i < out_features; i++) {
            this->bias[{0, i}] = distribution_uniform(generator);
            this->bias_grads[{0, i}] = 0;
        }

        // set initial cache values
        for (int i = 0; i < n_samples; i++) {
            for (int j = 0; j < in_features; j++) {
                this->cache[{i, j}] = 0;
            }
        }

    }

    ~ReLU() = default;

    virtual Matrix<T> forward(const Matrix<T> &x) override final {

        Matrix<T> out(1, this->out_features);

        for (int i = 0; i < this->out_features; i++) {

            // perform weighted sum
            T activation = 0;
            for (int w_i = 0; w_i < this->in_features; w_i++) {
                activation += this->weights[{w_i, i}] * x[{0,w_i}];
            }

            // add bias
            activation += this->bias[{0,i}];

            // apply ReLU nonlinearity
            activation = activation > 0 ? activation : 0;

            out[{0,i}] = activation;
        }

        return out;
    }

    virtual Matrix<T> backward(const Matrix<T> &dy) override final {
        Matrix<T> out(1, this->in_features);
        return out;
    }



};

template<typename T>
class Net {
    // Your implementation of the Net class starts here
};

// Function to calculate the loss
template<typename T>
T MSEloss(const Matrix<T> &y_true, const Matrix<T> &y_pred) {
    int n = y_true.getCols() * y_pred.getRows();
    double sum = 0;
    for (int i = 0; i < y_true.getRows(); i++) {
        for (int j = 0; j < y_true.getCols(); j++) {
            sum += std::pow((y_true[std::pair<int, int>(i, j)] - y_pred[std::pair<int, int>(i, j)]), 2.0);
        }
    }
    double loss = (1.0 / n) * sum;
    return loss;
}

// Function to calculate the gradients of the loss
template<typename T>
Matrix<T> MSEgrad(const Matrix<T> &y_true, const Matrix<T> &y_pred) {
    Matrix<double> gradient(y_true.getRows(), y_true.getCols());
    for (int i = 0; i < y_true.getRows(); i++) {
        for (int j = 0; j < y_true.getCols(); j++) {
            gradient[std::pair<int, int>(i, j)] =
                    2 * (y_true[std::pair<int, int>(i, j)] - y_pred[std::pair<int, int>(i, j)]);
        }
    }
    return gradient;
}

// Calculate the argmax
template<typename T>
Matrix<T> argmax(const Matrix<T> &y) {
    Matrix<T> maxMatrix = Matrix<T>(1, y.getRows());
    for (int i = 0; i < y.getRows(); i++) {
        int maxColIndex = 0;
        T maxColValue = y[{i, 0}];
        for (int j = 1; j < y.getCols(); j++) {
            T currValue = y[{i, j}];
            if (currValue > maxColValue) {
                maxColIndex = j;
                maxColValue = currValue;
            }
        }
        maxMatrix[std::pair<int, int>(0, i)] = maxColIndex;
    }
    return maxMatrix;
}

// Calculate the accuracy of the prediction, using the argmax
template<typename T>
T get_accuracy(const Matrix<T> &y_true, const Matrix<T> &y_pred) {
    auto predIndices = argmax(y_pred);
    auto trueIndices = argmax(y_true);
    T n = trueIndices.getCols();
    int accurateCount = 0;
    for (int i = 0; i < n; i++) {
        if (predIndices[{0, i}] == trueIndices[{0, i}]) {
            accurateCount++;
        }
    }
    T accuracy = (T) accurateCount / n;
    return accuracy;
}

template<typename T, typename U>
Matrix<typename std::common_type<T, U>::type> operator*(T x1, const Matrix<U> &x2) {
    Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix<U>(x2.getRows(), x2.getCols());
    for (int i = 0; i < x2.getRows(); i++) {
        for (int j = 0; j < x2.getCols(); j++) {
            newMatrix[{i, j}] = x1 * x2[{i, j}];
        }
    }
    return newMatrix;
}

bool test_matrix() {
    Matrix<int> matrix(2, 3);
    std::cout << "Vecotr with zeros" << std::endl;
    matrix.print();
    std::cout << "number of rows: " << matrix.getRows() << " number of cols: " << matrix.getCols() << std::endl;
    std::cout << "Transposed matrix:" << std::endl;
    matrix.transpose().print();

    Matrix<int> matrix2(2, 3, {1, 2, 3, 4, 5, 6});
    std::cout << "Initializing with list (row-major):" << std::endl;
    matrix2.print();

    std::cout << "Copy assignment matrix3 to matrix2:" << std::endl;
    Matrix<int> copied_matrix = matrix2;
    copied_matrix.print();

    Matrix<int> movedMatrix = Matrix<int>(2, 3, {1, 2, 3, 4, 5, 6});
    std::cout << "Move assignment matrix3 to matrix2:" << std::endl;
    Matrix<int> moved_matrix = std::move(movedMatrix);
    moved_matrix.print();

    std::cout << "Setting values 1 and 2" << std::endl;
    matrix[std::pair<int, int>(0, 0)] = 1;
    matrix[std::pair<int, int>(1, 2)] = 2;
    matrix.print();


    std::cout << "Multiplying scalar 2 by matrix [[1, 2, 3], [4, 5, 6]]: " << std::endl;
    (matrix2 * 2).print();
    (2 * matrix2).print();
    Matrix<int> vector(3, 1, {1, 2, 3});
    std::cout << "Multiplying vector [1, 2, 3] by matrix [[1, 2, 3], [4, 5, 6]]: " << std::endl;
    (matrix2 * vector).print();

    Matrix<int> matrix3(3, 2, {1, 2, 3, 4, 5, 6});
    std::cout << "Multiplying matrix by matrix [[1, 2, 3], [4, 5, 6]] * [[1, 2],  [3, 4], [5, 6]]: " << std::endl;
    (matrix2 * matrix3).print();

    std::cout << "Adding matrices [[1, 2, 3], [4, 5, 6]] and [[1, 0, 0], [0, 0, 2]]:" << std::endl;
    (matrix + matrix2).print();

    Matrix<int> biasMatrix(1, 3, {1, 2, 3});
    std::cout << "Adding matrix [[1], [2], [3]] to matrix[[1, 2, 3], [4, 5, 6]]: " << std::endl;
    (biasMatrix + matrix2).print();
    std::cout << "Reverse order of addition" << std::endl;
    (matrix2 + biasMatrix).print();

    std::cout << "Subtracting matrices [[1, 2, 3], [4, 5, 6]] and [[1, 0, 0], [0, 0, 2]]:" << std::endl;
    (matrix2 - matrix).print();

    std::cout << "Subtracting vector [[1], [2], [3]] from matrix [[1, 2, 3], [4, 5, 6]]:" << std::endl;
    (matrix2 - biasMatrix).print();

    std::cout << "Adding elements (1, 2) of matrices [[1, 2, 3], [4, 5, 6]] * [[1, 0, 0], [0, 0, 2]]: " << std::endl;
    std::cout << (matrix2[std::pair<int, int>(1, 2)] + matrix[std::pair<int, int>(1, 2)]) << std::endl;


    std::cout << "Constant indices test (1,2) of  [[1, 2, 3], [4, 5, 6]]" << std::endl;
    const Matrix<int> constant_matrix(2, 3, {1, 2, 3, 4, 5, 6});
    std::cout << constant_matrix[std::pair<int, int>(1, 2)] << std::endl;

    std::cout << "Move and copy operators" << std::endl;
    matrix3 = matrix2;
    std::cout << "Copying matrix2 into matrix3. Matrix2:" << matrix2.toString() << "Matrix3: " << matrix3.toString()
              << std::endl;
    matrix = std::move(matrix3);
    std::cout << "Moving matrix3 into matrix. Matrix: " << matrix.toString() << "Matrix3: " << matrix3.toString()
              << std::endl;

    std::cout << "Empty constructor:" << std::endl;
    Matrix<int> emptyMatrix;

    std::cout << "Add matrices of different types (should be 2.7):" << std::endl;
    Matrix<float> floatMatrix = Matrix<float>(1, 1, {1.5});
    Matrix<double> doubleMatrix = Matrix<double>(1, 1, {1.2});
    (floatMatrix + doubleMatrix).print();
    return true;
}

bool test_linear() {
//    Matrix<int> matrix(2, 3);

    // make layer
    Linear<double> l(3, 2, 1, 0);

    std::cout << "weights" << std::endl;
    // l.weights.print();

    std::cout << "biases" << std::endl;
    // l.bias.print();

    // input
    Matrix<double> x(1, 3);
    x[{0, 0}] = 1;
    x[{0, 1}] = 2;
    x[{0, 2}] = 3;

    std::cout << "input" << std::endl;
    x.print();

    Matrix<double> out = l.forward(x);

    std::cout << "output" << std::endl;
    out.print();

    return true;
}



bool test_MSE() {
    Matrix<double> y_pred = Matrix<double>(1, 4, {0.7, 1.0, 0.3, 0.0});
    Matrix<double> y_true = Matrix<double>(1, 4, {1.0, 1.0, 1.0, 0.0});

    std::cout << "Calculate first loss(should be 0.58/4 = 0.145) " << std::to_string(MSEloss(y_true, y_pred))
              << std::endl;
    std::cout << "Calculate gradient(should be [0.6, 0, 1.4, 0] ) " << (MSEgrad(y_true, y_pred)).toString()
              << std::endl;
    return true;
}

bool test_accuracy() {
    Matrix<double> y_pred = Matrix<double>(2, 2, {0.0, 1.0, 1.0, 0.0});
    Matrix<double> y_true = Matrix<double>(2, 2, {1.0, 0.0, 1.0, 0.0});

    std::cout << "Calculate the argmax (should be [1, 0]) " << argmax(y_pred).toString() << std::endl;
    std::cout << "Calculate accuracy (should be 0.5): " << std::to_string(get_accuracy(y_true, y_pred))
              << std::endl;
    return true;
}

int main(int argc, char *argv[]) {
    // Your training and testing of the Net class starts here
    test_linear();
    test_matrix();
    test_MSE();
    test_accuracy();
    return 0;
}
