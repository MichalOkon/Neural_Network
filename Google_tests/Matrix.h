//
// Created by micha on 2/1/2022.
//

#ifndef NEURAL_NETWORK_MATRIX_H
#define NEURAL_NETWORK_MATRIX_H


class type;

template<typename T>
class Matrix {
private:

    std::vector<T> data{};
    int rows;
    int cols;

public:

    Matrix(int rows, int cols);

    Matrix(int rows, int cols, const std::initializer_list<T> &list);

    Matrix(Matrix &other) : data(other.data), rows(other.rows), cols(other.cols) {}

    Matrix(Matrix &&other) noexcept;

    ~Matrix() = default;

    // Copy assignment operator
    Matrix<T> &operator=(const Matrix &other);

    // Move assignment operator
    Matrix<T> &&operator=(const Matrix &&other) noexcept ;

    T &operator[](const std::pair<int, int> &ij);

    const T &operator[](const std::pair<int, int> &ij) const;

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator*(U x) const;

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator*(const Matrix<U> &B) const;

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator+(const Matrix<U> &B) const;

    template<typename U>
    Matrix<typename std::common_type<T, U>::type> operator-(const Matrix<U> &B) const;

    Matrix transpose() const;

    int getRows() const;

    int getCols() const;

    // For debugging
    void print() const;

    std::string toString() const;
};


#endif //NEURAL_NETWORK_MATRIX_H
