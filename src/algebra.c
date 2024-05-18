#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if((a.rows!=b.rows)||(a.cols!=b.cols)){
        // 如果a和b的行数或列数不相等，会给出错误提示"Error: Matrix a and b must have the same rows and cols.\n"并返回一个空矩阵
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }else{
        Matrix result = create_matrix(a.rows, a.cols);
        // 执行矩阵加法
        for(int i=0;i<a.rows;i++){
            for(int j=0;j<a.cols;j++){
                result.data[i][j] = a.data[i][j] + b.data[i][j];
            }
        }
        return result;
    }
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if((a.rows!=b.rows)||(a.cols!=b.cols)){
        // 如果a和b的行数或列数不相等，会给出错误提示"Error: Matrix a and b must have the same rows and cols.\n"并返回一个空矩阵
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }else{
        Matrix result = create_matrix(a.rows, a.cols);
        // 执行矩阵减法
        for(int i=0;i<a.rows;i++){
            for(int j=0;j<a.cols;j++){
                result.data[i][j] = a.data[i][j] - b.data[i][j];
            }
        }
        return result;
    }
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if(a.cols!=b.rows){
        //如果a的列数不等于b的行数，会给出错误提示"Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n"并返回一个空矩阵
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }else{
        Matrix result = create_matrix(a.rows, b.cols);
        //初始化输出矩阵
        for(int i=0;i<a.rows;i++){
            for(int j=0;j<b.cols;j++){
                result.data[i][j] = 0;
            }
        }
        //执行矩阵乘法
        double temp;
        for(int i=0;i<a.rows;i++){
            for(int j=0;j<b.cols;j++){
                for(int k=0;k<a.cols;k++){
                    temp = a.data[i][k]*b.data[k][j];
                    result.data[i][j] += temp;
                }
            }
        }
        return result;
    }
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result = create_matrix(a.rows, a.cols);
    // 执行矩阵的数乘
    for(int i=0;i<a.rows;i++){
        for(int j=0;j<a.cols;j++){
            result.data[i][j] = a.data[i][j]*k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result = create_matrix(a.rows, a.cols);
    // 执行矩阵的转置
    for(int i=0;i<a.rows;i++){
        for(int j=0;j<a.cols;j++){
            result.data[i][j] = a.data[j][i];
        }
    }
    return result;
}

double det_matrix(Matrix a) 
{
    if (a.rows != a.cols) {
        //如果a不是方阵，会给出错误提示"Error: The matrix must be a square matrix.\n"并返回0
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    if (a.rows == 1) {
        return a.data[0][0];
    }
    double det = 0.0;
    for (int i = 0; i < a.cols; i++) {
        Matrix subMatrix = create_matrix(a.rows - 1, a.cols - 1);
        for (int j = 1; j < a.rows; j++) {
            for (int k = 0, subCol = 0; k < a.cols; k++) {
                if (k == i) continue; // 跳过当前列
                subMatrix.data[j-1][subCol] = a.data[j][k];
                subCol++;
            }
        }

        // 计算行列式的当前项（沿第一列展开）
        double sign = (i % 2 == 0) ? 1.0 : -1.0; // 交替符号因子
        det += sign * a.data[0][i] * det_matrix(subMatrix);
    }
    return det;
}

Matrix inv_matrix(Matrix a)
{
    // 声明一个矩阵用于存放逆矩阵
    Matrix inverse = create_matrix(a.rows, a.cols);
    double det = det_matrix(a);
    if (a.rows != a.cols) {
        //如果a不是方阵，会给出错误提示"Error: The matrix must be a square matrix.\n"并返回0
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    if(det_matrix(a)==0){
        //如果a的逆矩阵不存在，会给出错误提示"Error: The matrix is singular.\n"并返回一个空矩阵。
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    // 计算伴随矩阵
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            Matrix subMatrix = create_matrix(a.rows - 1, a.cols - 1);
            int row_mark = 0;
            for (int k = 0; k < a.rows; k++) {
                if (k == i) continue; // 跳过当前行
                int col_mark = 0;
                for (int l = 0; l < a.cols; l++) {
                    if (l == j) continue; // 跳过当前列
                    subMatrix.data[row_mark][col_mark] = a.data[k][l];
                    col_mark++;
                }
                row_mark++;
            }
            // 计算代数余子式并将其赋值给伴随矩阵转置位置
            double cofactor = ((i + j) % 2 == 0 ? 1 : -1) * det_matrix(subMatrix);
            inverse.data[j][i] = cofactor / det;
        }
    }
    return inverse;
}

int rank_matrix(Matrix a) {
    int rank = 0;
    for (int col = 0; col < a.cols; ++col) {
        int pivot = rank;
        for (int row = rank; row < a.rows; ++row) {
            if (fabs(a.data[row][col]) > fabs(a.data[pivot][col])) {
                pivot = row;
            }
        }
        if (fabs(a.data[pivot][col]) <= 1e-6) {
            // 此列是零列，跳过
            continue;
        }

        // 将最大绝对值行交换到当前行
        for (int k = 0; k < a.cols; ++k) {
            double temp = a.data[rank][k];
            a.data[rank][k] = a.data[pivot][k];
            a.data[pivot][k] = temp;
        }

        // 使用当前行消去下面所有行的此列元素
        for (int row = rank + 1; row < a.rows; ++row) {
            double factor = a.data[row][col] / a.data[rank][col];
            for (int k = col; k < a.cols; ++k) {
                a.data[row][k] -= a.data[rank][k] * factor;
            }
        }

        // 增加秩
        ++rank;

        // 如果达到了行或列的最大数，就结束
        if (rank == a.rows || rank == a.cols) {
            break;
        }
    }
    return rank;
}

double trace_matrix(Matrix a)
{
    if(a.rows!=a.cols){
        //如果a不是方阵，会给出错误提示"Error: The matrix must be a square matrix.\n"并返回0
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }else{
        int result=0;
        for(int i=0;i<a.rows;i++){
            result += a.data[i][i];
        }
        return result;
    }
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}