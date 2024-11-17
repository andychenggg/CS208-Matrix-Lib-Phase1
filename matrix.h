#include <stdio.h>
#ifndef __C__MATRIX__H__
#define __C__MATRIX__H__


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part1:
 定义错误编号，
 在使用库的过程中出现错误，会有相应的错误编号提示，可以根据下面对应的数字找到对应的错误类型

///////////////////////////////////////////////////////////////////////////////////////////////////*/

#define ERROR_IN_ALLOCATING_MEMORY 0x00                                //内存申请分配失败
#define ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO 0x01                   //矩阵边长,float指针的size为负数
#define ERROR_CAUSE_BY_NULL_POINTER 0x02                               //空指针异常
#define ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix 0x03                 //矩阵并非使用库中的方法创建
#define ERROR_INDEX_OUT_OF_BOUND 0x04                                  //访问越界
#define ERROR_MODIFICATION_IN_MATRIX 0x05                              //用户可能自己将Matrix的性质设置为非法值
#define ERROR_NOT_THE_SAME_SIZE 0x06                                   //复制，加减法要求连个矩阵是相同的规模
#define ERROR_INVALID_COLS_AND_ROWS_IN_MULTIPLICATION 0x07             //矩阵乘法要求第一个矩阵列数等于第二个矩阵的行数
#define ERROR_INVALID_MATRIX_ELEMENT_FROM_KEYBOARD 0x08                //从键盘输入的矩阵元素不合法
#define ERROR_INVALID_NUMBER_OF_LINES_FROM_KEYBOARD 0x09               //从键盘输入的行数不合法
#define ERROR_FAIL_TO_OPEN_FILE 0x0a                                   //文件打开失败
#define ERROR_INVALID_MATRIX_ELEMENT_FROM_FILE 0x0b                    //从文件输入的矩阵元素不合法
#define ERROR_UNSUPPORTED_MODE 0x0c                                    //不支持该种模式的文件操作
#define ERROR_NOT_A_SQUARE_MATRIX 0x0d                                 //有一些计算如求行列式只能对方阵进行运算
#define ERROR_MATRIX_IRREVERSIBLE 0x0e                                 //矩阵不可逆
typedef unsigned ERROR_ID;

/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part1:
 声明矩阵的结构体，

///////////////////////////////////////////////////////////////////////////////////////////////////*/
typedef struct mat Matrix;




/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part4:
 声明矩阵的创建方法，初始化方法，设置方法，删除方法，复制方法，打印方法等等

///////////////////////////////////////////////////////////////////////////////////////////////////*/

Matrix* createMatrix(long long rows, long long cols);
Matrix* createMatrix_with_pFloat( long long rows, long long cols, const float *pFloat, long long arr_size);
Matrix* createMatrix_from_keyboard(long long rows, long long cols);
Matrix* createMatrix_with_a_file(long long rows, long long cols, const char* filename_with_directory);
Matrix* create_an_identity_matrix(long long rows);
Matrix* create_an_zero_matrix(long long rows, long long cols);
Matrix* create_an_one_matrix(long long rows, long long cols);
void set_Matrix_same_value(Matrix * matrix, float same_value);
void set_Matrix_with_an_array(Matrix *matrix, const float *pFloat, long long arr_size);
void set_Matrix_from_keyboard(Matrix *pMat);
void set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory);
void set_Matrix(Matrix *matrix, long long row, long long col, float value);
void set_whole_Matrix_self_define(Matrix *matrix, float (*self_define_function)(float float_element));
void deleteMatrix(Matrix **p_pMat);
void delete_all_Matrix();
void copyMatrix(const Matrix* pMat, Matrix* des);
Matrix *copyMatrix_and_return(const Matrix* pMat);
void printMatrix(const Matrix*);
void printMatrix_to_a_file(const Matrix* pMat, const char* filename_with_directory, const char* mode);
void print_matrices_list();


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part5:
 求矩阵的一些特性，不会对矩阵做出任何的改动

///////////////////////////////////////////////////////////////////////////////////////////////////*/

float minimum_in_a_matrix(const Matrix *);
float maximum_in_a_matrix(const Matrix *);
float determinant(const Matrix *);
float trace_matrix(const Matrix *);
size_t rank(const Matrix *);


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part6:
 进行一些矩阵的运算，大部分运算都提供了两种方法：
 一种是动态申请一个新的矩阵储存结果并返回（上半部分）
 一种是对原有的某一个参数中的矩阵进行修改（下半部分）

///////////////////////////////////////////////////////////////////////////////////////////////////*/
Matrix* addConstMatrices(const Matrix *, ...);
Matrix* subtractConstMatrices(const Matrix *, ...);
Matrix* MultiplyConstMatrices(const Matrix *, ...);
Matrix *add_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar);
Matrix * subtract_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar);
Matrix * multiply_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar);
Matrix *InvertConstMatrix(const Matrix *);
Matrix *transposeConstMatrix(const Matrix*);
Matrix *adj_const_matrix(const Matrix*);
Matrix * swapRows_const_matrix(const Matrix *, size_t rows1, size_t rows2);
Matrix * addRows_const_matrix(const Matrix *, size_t rows1, size_t rows2, float coefficient);
Matrix * multiplyRows_const_matrix(const Matrix *, size_t rows1, float coefficient);
Matrix *solve_equation_and_return(const Matrix *);


void addMatrix(Matrix* pMat, ...);
void subtractMatrix(Matrix* pMat, ...);
void MultiplyMatrix(Matrix* pMat, ...);
void add_a_scalar(Matrix *, float );
void subtract_a_scalar(Matrix *, float );
void multiply_a_scalar(Matrix *, float );
void invertMatrix(Matrix*);
void transposeMatrix(Matrix*);
void adj_matrix(Matrix*);
void swapRows(Matrix *, size_t rows1, size_t rows2);
void addRows(Matrix *, size_t rows1, size_t rows2, float coefficient);
void multiplyRows(Matrix *, size_t rows1, float coefficient);
void solve_equation(Matrix *);
#endif
