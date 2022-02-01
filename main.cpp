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
        rows = rows;
        cols = cols;
    }

    Matrix(int rows, int cols, const std::initializer_list<T> &list) : rows(rows), cols(cols) {
        if (rows * cols != list.size()) {
            throw std::length_error(
                    "The number of columns and rows does not match the given list in the matrix construction");
        }
        // Check if that actually works
        data = std::vector(list);
        rows = rows;
        cols = cols;
    }

    Matrix(Matrix &other) : data(other.data), rows(other.rows), cols(other.cols) {}

    Matrix(Matrix &&other) noexcept: data(other.data), rows(other.rows), cols(other.cols) {
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
    Matrix &&operator=(const Matrix &&other) {
        if (this != &other) {
            data = other.data;
            rows = other.rows;
            cols = other.cols;
            other.data = nullptr;
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
        Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix(this->rows, this->columns);
        for (int i = 0; i < this->data.size; i++) {
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
                    sum += this->data[i* cols + k] * B.data[k * B.cols + j];
                }
                newMatrix.data[i * B.cols + j] = sum;
            }
        }
        return newMatrix;
    }

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator+(const Matrix<U> &B) const {
        if (B.rows != 1 && (B.rows != this->rows || B.cols != this->cols)) {
            throw std::length_error("Matrices' sizes incompatible for addition");
        }
        Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix<typename std::common_type<T, U>::type>(
                this->rows, this->cols);
        if (B.rows == 1) {
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->cols; j++) {
                    newMatrix.data[i * cols + j] = this->data[i * cols + j] + B.data[j];
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
        Matrix<std::common_type<T, U>> newMatrix = Matrix<std::common_type<T, U>>(this->rows, this->cols);
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
            result += "[";
            for (int j = 0; j < cols; j++) {
                result += this->data[i * cols + j];
                if (j != cols - 1) {
                    result += ", ";
                }
            }
            result += "]";
            if (i != rows - 1) {
                result += ", ";
            }
        }
        result += "]";
        return result;
    }
};

template<typename T>
class Layer {
    // Your implementation of the Layer clas starts here
};

template<typename T>
class Linear : public Layer<T> {
    // Your implementation of the Linear class starts here
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

    std::cout << "Setting values" << std::endl;
    matrix[std::pair<int, int>(0, 0)] = 9;
    matrix[std::pair<int, int>(1, 2)] = 7;
    matrix.print();

    std::cout << "Adding vectors:" << std::endl;
    (matrix + matrix2).print();

    Matrix<int> vector(3, 1, {1, 2, 3});
    std::cout << "Multiplying vector [1, 2, 3] by matrix [[1, 2, 3], [4, 5, 6]]: " <<std::endl;
    (matrix2 * vector).print();

    Matrix<int> matrix3(3, 2, {1, 2, 3, 4, 5, 6});
    std::cout << "Multiplying matrix by matrix [[1, 2, 3], [4, 5, 6]] * [[1, 2],  [3, 4], [5, 6]]: " <<std::endl;
    (matrix2 * matrix3).print();
    return true;
}

int main(int argc, char *argv[]) {
    // Your training and testing of the Net class starts here
    // Testing the matrix
    test_matrix();
    return 0;
}
