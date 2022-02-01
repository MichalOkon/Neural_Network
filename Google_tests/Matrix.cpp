//
// Created by micha on 2/1/2022.
//

#include <stdexcept>
#include <vector>
#include <iostream>
#include "Matrix.h"

template<typename T>
Matrix<T>::Matrix(int rows, int cols) : rows(rows), cols(cols) {
    data = std::vector<T>(rows * cols, 0);
}

template<typename T>
Matrix<T>::Matrix(int rows, int cols, const std::initializer_list <T> &list) : rows(rows), cols(cols) {
    if (rows * cols != list.size()) {
        throw std::length_error(
                "The number of columns and rows does not match the given list in the matrix construction");
    }
    // Check if that actually works
    data = std::vector(list);
}

template<typename T>
Matrix<T>::Matrix(Matrix &&other) noexcept: data(other.data), rows(other.rows), cols(other.cols) {
    std::destroy(other.data.begin(), other.data.end());
    other.rows = 0;
    other.cols = 0;
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix &other) {
    if (this != &other) {
        data = other.data;
        rows = other.rows;
        cols = other.cols;
    }
    return *this;
}

template<typename T>
Matrix<T> &&Matrix<T>::operator=(const Matrix &&other) noexcept {
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

template<typename T>
T &Matrix<T>::operator[](const std::pair<int, int> &ij) {
    if (ij.first >= rows) {
        throw std::length_error("Number of rows exceeded");
    }
    if (ij.second >= cols) {
        throw std::length_error("Number of columns exceeded");
    }
    return data[ij.first * cols + ij.second];
}

template<typename T>
const T &Matrix<T>::operator[](const std::pair<int, int> &ij) const {
    if (ij.first >= rows) {
        throw std::length_error("Number of rows exceeded");
    }
    if (ij.second >= cols) {
        throw std::length_error("Number of columns exceeded");
    }
    return data[ij.first * cols + ij.second];
}

template<typename T>
template<typename U>
Matrix<typename std::common_type<T,U>::type> Matrix<T>::operator*(U x) const {
    Matrix<typename std::common_type<T, U>::type> newMatrix = Matrix(this->rows, this->columns);
    for (int i = 0; i < this->data.size; i++) {
        newMatrix.data[i] = this->data[i] * x;
    }
    return newMatrix;
}

template<typename T>
template<typename U>
Matrix<typename std::common_type<T, U>::type> Matrix<T>::operator+(const Matrix<U> &B) const {
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

template<typename T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> transposedMatrix(cols, rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            transposedMatrix.data[j * rows + i] = this->data[i * cols + j];
        }
    }

    return transposedMatrix;
}

template<typename T>
template<typename U>
Matrix<typename std::common_type<T, U>::type> Matrix<T>::operator-(const Matrix<U> &B) const {
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

template<typename T>
std::string Matrix<T>::toString() const {
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

template<typename T>
int Matrix<T>::getCols() const {
    return cols;
}

template<typename T>
template<typename U>
Matrix<typename std::common_type<T, U>::type> Matrix<T>::operator*(const Matrix<U> &B) const {
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

template<typename T>
int Matrix<T>::getRows() const {
    return rows;
}

template<typename T>
void Matrix<T>::print() const {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << this->data[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
}