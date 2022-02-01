#include "gtest/gtest.h"
#include "Matrix.h"

TEST(DummyTest, TEST_1_IS_1) {
    EXPECT_EQ(1, 1);
}

TEST(Test, TEST_ZERO_VECTOR) {
    Matrix<int> matrix1 (2, 3);
    EXPECT_EQ(matrix1.toString(), "[[0, 0, 0], [0, 0, 0]]");
}