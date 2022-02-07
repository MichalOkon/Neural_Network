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

    std::vector<T> data;
    int rows;
    int cols;

public:

    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        data = std::vector<T>(rows * cols, 0);
    }

    Matrix(int rows, int cols, const std::initializer_list<T> &list) : rows(rows), cols(cols) {
        if (rows * cols != (int)list.size()) {
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
        for (int i = 0; i < (int)this->data.size(); i++) {
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
    Matrix<float> cache;
    Matrix<float> bias;
    Matrix<float> weights;
    Matrix<float> bias_grads;
    Matrix<float> weight_grads;
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
        this->weights(in_features, out_features);
        this->weight_grads(in_features, out_features);
        this->bias(1, out_features);
        this->bias_grads(1, out_features);
        this->cache(n_samples, in_features);

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

        for (int i = 0; i < x.getRows(); i++) {
            for (int j = 0; j < x.getCols(); j++) {
                // cache the input value
                this->cache[{i, j}] = x[{i, j}];

                // TODO not sure if this should be bias[0,j] or bias[0,i]
                x[{i, j}] = x[{i, j}] * this->weights[{i, j}] + this->bias[{0, j}];
            }
        }


    }

    virtual Matrix<T> backward(const Matrix<T> &dy) override final {

    }
};

template<typename T>
class ReLU : public Layer<T> {
    // Your implementation of the ReLU class starts here
};

template<typename T>
class Net {
    // Your implementation of the Net class starts here
};

// Function to calculate the loss
template<typename T>
T MSEloss(const Matrix<T> &y_true, const Matrix<T> &y_pred) {
    // Your implementation of the MSEloss function starts here
}

// Function to calculate the gradients of the loss
template<typename T>
Matrix<T> MSEgrad(const Matrix<T> &y_true, const Matrix<T> &y_pred) {
    // Your implementation of the MSEgrad function starts here
}

// Calculate the argmax
template<typename T>
Matrix<T> argmax(const Matrix<T> &y) {
    // Your implementation of the argmax function starts here
}

// Calculate the accuracy of the prediction, using the argmax
template<typename T>
T get_accuracy(const Matrix<T> &y_true, const Matrix<T> &y_pred) {
    // Your implementation of the get_accuracy starts here
}

template<typename T, typename U>
Matrix<typename std::common_type<T, U>::type> operator*(T x1, const Matrix<U> &x2) {
    Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix<U>(x2.getRows(), x2.getCols());
    std::initializer_list<U> list;
    for (int i = 0; i < x2.getRows(); i++) {
        for (int j = 0; j < x2.getCols(); j++) {
            newMatrix[std::pair<int, int>(i, j)] = x1 * x2[std::pair<int, int>(i, j)];
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

    std::cout << "Move assignment matrix3 to matrix2:" << std::endl;
    Matrix<int> moved_matrix = std::move(Matrix<int>(2, 3, {1, 2, 3, 4, 5, 6}));
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

    return true;
}

int main(int argc, char *argv[]) {
    // Your training and testing of the Net class starts here
    // Testing the matrix
    test_matrix();
    return 0;
}
