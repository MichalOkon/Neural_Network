#include <cmath>
#include <initializer_list>
#include <iostream>
#include <list>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>


template<typename T>
class Matrix {
private:

    // Data held in row-major notation
    std::vector<T> data;
    int rows;
    int cols;

public:

    template<typename U> friend
    class Matrix;

    Matrix() {};

    // Initializes new matrix with 0s
    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        data = std::vector<T>(rows * cols, 0);
    }

    // Initializes new matrix with values from the initializer_list
    Matrix(int rows, int cols, const std::initializer_list<T> &list) : rows(rows), cols(cols) {
        if (rows * cols != (int) list.size()) {
            throw std::length_error(
                    "The number of columns and rows does not match the given list in the matrix construction");
        }

        data = std::vector(list);
    }

    // Copy constructor
    Matrix(const Matrix &other) : data(other.data), rows(other.rows), cols(other.cols) {}

    // Move constructor
    Matrix(Matrix &&other) noexcept: data(std::move(other.data)), rows(other.rows), cols(other.cols) {
        std::destroy(other.data.begin(), other.data.end());
        other.rows = 0;
        other.cols = 0;
    }

    ~Matrix() = default;

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

    // Returns value at the position (i, j)
    T &operator[](const std::pair<int, int> &ij) {
        if (ij.first >= rows) {
            throw std::length_error("Number of rows exceeded");
        }
        if (ij.second >= cols) {
            throw std::length_error("Number of columns exceeded");
        }
        return data[ij.first * cols + ij.second];
    }

    // Returns value at the position (i, j) as a constant
    const T &operator[](const std::pair<int, int> &ij) const {
        if (ij.first >= rows) {
            throw std::length_error("Number of rows exceeded");
        }
        if (ij.second >= cols) {
            throw std::length_error("Number of columns exceeded");
        }
        return data[ij.first * cols + ij.second];
    }

    // Matrix-scalar multiplication
    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator*(U x) const {
        Matrix<typename std::common_type<T, U>::type> newMatrix(this->rows, this->cols);
        for (int i = 0; i < (int) this->data.size(); i++) {
            newMatrix.data[i] = this->data[i] * x;
        }
        return newMatrix;
    }

    // Matrix-matrix multiplication
    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator*(const Matrix<U> &B) const {
        if (this->cols != B.getRows()) {
            throw std::length_error("Matrices' sizes incompatible for multiplication");
        }
        Matrix<typename std::common_type<T, U>::type> newMatrix(this->rows, B.getCols());
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < B.getCols(); j++) {
                typename std::common_type<T, U>::type sum = 0;
                for (int k = 0; k < this->cols; k++) {
                    sum += this->data[i * cols + k] * B.data[k * B.cols + j];
                }
                newMatrix.data[i * B.cols + j] = sum;
            }
        }
        return newMatrix;
    }

    // Matrix-matrix addition
    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator+(const Matrix<U> &B) const {
        if (B.rows != 1 && rows != 1 && (B.rows != this->rows || B.cols != this->cols)) {
            throw std::length_error("Matrices' sizes incompatible for addition");
        }
        Matrix<typename std::common_type<T, U>::type> newMatrix(std::max(this->rows, B.rows),
                                                                std::max(this->cols, B.cols));
        // Consider special cases where one vector is a bias vector and needs to be broadcast
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

    // Matrix-matrix subtraction
    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator-(const Matrix<U> &B) const {
        if (B.rows != 1 && (B.rows != this->rows || B.cols != this->cols)) {
            throw std::length_error("Matrices' sizes incompatible for subtraction");
        }
        Matrix<typename std::common_type<T, U>::type> newMatrix(this->rows, this->cols);
        // Consider the special where the second vector is a bias vector and needs to be broadcast
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

    // Transpose the matrix
    Matrix transpose() const {
        Matrix<T> transposedMatrix(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposedMatrix.data[j * rows + i] = this->data[i * cols + j];
            }
        }

        return transposedMatrix;
    }

    // Return the number of rows in the matrix
    int getRows() const {
        return rows;
    }

    // Return the number of colums in the matrix
    int getCols() const {
        return cols;
    }

    // For debugging, prints the matrix
    void print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << this->data[i * cols + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Returns the matrix as a string
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

        // cache input value for backpropagation
        this->cache = x;

        return x * this->weights + this->bias;
    }

    virtual Matrix<T> backward(const Matrix<T> &dy) override final {
        // update the weights gradient, element-wise multiplication, multiply each gradient
        this->weight_grads = this->cache.transpose() * dy;
        // update the bias gradient
        for (int i = 0; i < (int) dy.getCols(); i++) {
            T sum = 0;
            for (int j = 0; j < (int) dy.getRows(); j++) {
                sum += dy[{j, i}];
            }
            this->bias_grads[{0, i}] = sum;
        }
        Matrix<T> out = dy * this->weights.transpose();
        // calculate the downstream gradient and return it
        return out;
    }

    void optimize(T learning_rate) {
        this->weights = this->weights - weight_grads * learning_rate;
        this->bias = this->bias - bias_grads * learning_rate;
    }

    Matrix<T> weights;
    Matrix<T> bias;
};

template<typename T>
class ReLU : public Layer<T> {

private:
    Matrix<T> cache;
    int n_samples;
    int in_features;
    int out_features;

public:

    ReLU(const int in_features,
         const int out_features,
         const int n_samples) : n_samples(n_samples), in_features(in_features), out_features(out_features) {

        // create matrices
        this->cache = Matrix<T>(n_samples, in_features);

        // set initial cache values
        for (int i = 0; i < n_samples; i++) {
            for (int j = 0; j < in_features; j++) {
                this->cache[{i, j}] = 0;
            }
        }

    }

    ~ReLU() = default;

    virtual Matrix<T> forward(const Matrix<T> &x) override final {
        // cache input value for backpropagation
        this->cache = x;

        // apply reLu nonlinearity
        Matrix<T> out = x;
        for (int row = 0; row < out.getRows(); row++) {
            for (int col = 0; col < out.getCols(); col++) {
                if (x[{row, col}] < 0) {
                    out[{row, col}] = 0;
                }
            }
        }

        return out;
    }

    virtual Matrix<T> backward(const Matrix<T> &dy) override final {

        Matrix<T> out = dy;

        for (int row = 0; row < out.getRows(); row++) {
            for (int col = 0; col < out.getCols(); col++) {
                if (dy[{row, col}] < 0) {
                    out[{row, col}] = 0;
                } else {
                    out[{row, col}] *= 1;
                }
            }
        }

        return out;
    }
};

template<typename T>
class Net {

public:

    Linear<T> in_layer;
    ReLU<T> hidden_layer;
    Linear<T> out_layer;

    Net(int in_features, int hidden_dim, int out_features, int n_samples, int seed)
            : in_layer(Linear<T>(in_features, hidden_dim, n_samples, seed)),
              hidden_layer(ReLU<T>(hidden_dim, hidden_dim, n_samples)),
              out_layer(Linear<T>(hidden_dim, out_features, n_samples, seed)) {}

    ~Net() = default;

    Matrix<T> forward(const Matrix<T> &x) {
        Matrix<T> x_cp(x);
        x_cp = in_layer.forward(x_cp);
        x_cp = hidden_layer.forward(x_cp);
        x_cp = out_layer.forward(x_cp);
        return x_cp;
    }

    Matrix<T> backward(const Matrix<T> &dy) {
        Matrix<T> dy_cp(dy);
        dy_cp = out_layer.backward(dy_cp);
        dy_cp = hidden_layer.backward(dy_cp);
        dy_cp = in_layer.backward(dy_cp);
        return dy_cp;
    }

    void optimize(T learning_rate) {
        in_layer.optimize(learning_rate);
        out_layer.optimize(learning_rate);
    }
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
            gradient[{i, j}] = 2 * (y_true[{i, j}] - y_pred[{i, j}]);
        }
    }
    return gradient;
}

// Calculate the one-row matrix containing the index of the column with the highest value in each row
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
    int n = trueIndices.getCols();
    int accurateCount = 0;
    for (int i = 0; i < n; i++) {
        if (predIndices[{0, i}] == trueIndices[{0, i}]) {
            accurateCount++;
        }
    }
    T accuracy = (T) accurateCount / n;
    return accuracy;
}

// Scalar-matrix multiplication
template<typename T, typename U>
Matrix<typename std::common_type<T, U>::type> operator*(T x1, const Matrix<U> &x2) {
    Matrix<typename std::common_type<T, U>::type> newMatrix(x2.getRows(), x2.getCols());
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

    std::cout << "Add matrices of different types (should be 4.5):" << std::endl;
    Matrix<float> floatMatrix = Matrix<float>(1, 1, {2.0});
    Matrix<double> doubleMatrix = Matrix<double>(1, 1, {2.5});
    (floatMatrix + doubleMatrix).print();

    std::cout << "Multiplying matrices of different types (should be 5.0):" << std::endl;
    (floatMatrix * doubleMatrix).print();

    std::cout << "Subtracting matrices of different types (should be -0.5):" << std::endl;
    (floatMatrix - doubleMatrix).print();
    return true;
}

bool test_linear() {
//    Matrix<int> matrix(2, 3);

    // make layer
    Linear<double> l(3, 2, 1, 0);

    std::cout << "weights" << std::endl;
    l.weights.print();

    std::cout << "biases" << std::endl;
    l.bias.print();

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
    // test_linear();
    // test_matrix();
    // test_MSE();
    // test_accuracy();

    double learning_rate = 0.0005;
    int optimizer_steps = 100;
    int seed = 1;

    int in_features = 2;
    int hidden_dim = 100;
    int out_features = 2;
    Net<double> net(in_features, hidden_dim, out_features, 4, seed);

    Matrix<double> x_xor(4, 2);
    x_xor[{1, 1}] = 1;
    x_xor[{2, 0}] = 1;
    x_xor[{3, 0}] = 1;
    x_xor[{3, 1}] = 1;
    std::cout << "x_xor" << std::endl;
    x_xor.print();

    Matrix<double> y_true(4, 2);
    y_true[{0, 0}] = 1;
    y_true[{1, 1}] = 1;
    y_true[{2, 1}] = 1;
    y_true[{3, 0}] = 1;
    std::cout << "y_true" << std::endl;
    y_true.print();

    std::cout << std::endl;

    // optimize
    for (int i=0; i<optimizer_steps; i++) {

        Matrix<double> y_pred = net.forward(x_xor);
        Matrix<double> mse_grad = MSEgrad(y_true, y_pred);

        double accuracy = get_accuracy(y_true, y_pred);
        std::cout << accuracy << std::endl;

        net.backward(mse_grad);
        net.optimize(learning_rate);
    }

    std::cout << "prediction: " << std::endl;
    Matrix<double> out = net.forward(x_xor);

    out.print();
    argmax(out).print();

    return 0;
}
