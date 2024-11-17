#include "matrix.h"
#include <malloc.h>
#include <string.h>
#include <stdbool.h>
#include <minmax.h>
#include <stdarg.h>

#define SIZE_OF_POINTER (sizeof(NULL))

/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part1:
 将mat的定义放在这里；
 为监测用户行为创建的数据结构:矩阵链表;
 只允许在本文件中访问，不允许和外部建立链接,所以调用能够忽略很多非法的问题。
 更好的进行内存释放和检查;

///////////////////////////////////////////////////////////////////////////////////////////////////*/
struct mat {
    size_t rows_numbers_in_Matrix;
    size_t cols_numbers_in_Matrix;
    float *pointer_of_data;
};
typedef struct mat_node {
    Matrix *item;
    struct mat_node *before;
    struct mat_node *next;
} Matrix_Node;
typedef struct mat_linked_list {
    size_t Matrix_counts;
    Matrix_Node *head;
    Matrix_Node *tail;
} Matrix_List;


//只允许本文件访问
static Matrix_List matrices_list = {0, 0, 0};

static void add_node_to_list(Matrix_Node *pNode);

static Matrix_Node *create_node(Matrix *pMat);

static void delete_node(Matrix_Node *pNode);

static void deleteAll();

static Matrix_Node *find_node_by_matrix(const Matrix *pMat);

static bool has_before(const Matrix_Node *pNode);

static bool has_next(const Matrix_Node *pNode);


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part2:
 核心算法的声明：
 static方法，仅供内部调用，无需检查数据合法性
 能够更加简洁地看清楚不同矩阵功能的具体实现

///////////////////////////////////////////////////////////////////////////////////////////////////*/





static void add_two_mat(Matrix* des, const Matrix *pMat);

static void subtract_two_mat(Matrix* des, const Matrix *pMat);

static void multiply_two_mat(Matrix* des, const Matrix *pMat);

static float det(const float * pFloat, size_t size);

static void cofactor(const float * pFloat,float * des, size_t size, size_t index_i, size_t index_j);

static void invert(const float *pFloat, float * des, size_t size);

static void transpose(const float *pFloat, float *des, size_t rows, size_t cols);

static void cofactor_Matrix(const float *pFloat,  float *des, size_t sides);

static void adjoint_Matrix(const float *pFloat, float *des, size_t sides);

static float trace(const float *pFloat, size_t sides);

static void Gause_elimination_to_upper_triangle_mat(float *augmented_pFloat, size_t rows, size_t cols);

static size_t triangle_mat_rank_for_aug(const float *simplified_augmented_pFloat, size_t rows, size_t cols);

static size_t triangle_mat_rank_for_cof(const float *simplified_augmented_pFloat, size_t rows, size_t cols);

static void Gause_elimination_to_solve_equations(float *simplified_augmented_pFloat, size_t cols, size_t rank);

static void swap_rows(float *pFloat, size_t cols, size_t i1, size_t i2);

static void add_rows(float  *pFloat, size_t cols, size_t i1, size_t i2, float coefficient);

static void multiply_rows(float *pFloat, size_t cols, size_t i1, float coefficient);


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part3:
 异常处理的声明


///////////////////////////////////////////////////////////////////////////////////////////////////*/

static bool check_error(const char *function_name, const char *description1, const void *errors_key_numbers, ERROR_ID id);


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part1:
 将mat的定义放在这里；
 为监测用户行为创建的数据结构:矩阵链表;
 只允许在本文件中访问，不允许和外部建立链接,所以调用能够忽略很多非法的问题。
 更好的进行内存释放和检查;

///////////////////////////////////////////////////////////////////////////////////////////////////*/

static void add_node_to_list(Matrix_Node *pNode) {
    if (!matrices_list.head) {
        //head == null
        matrices_list.head = (Matrix_Node *) calloc(1, 3 * SIZE_OF_POINTER);
        matrices_list.tail = matrices_list.head;
    }
    matrices_list.tail->next = (Matrix_Node *) pNode;
    pNode->before = matrices_list.tail;
    matrices_list.tail = pNode;
    matrices_list.Matrix_counts++;
}

static Matrix_Node *create_node(Matrix *pMat) {
    Matrix_Node *pMatNode = (Matrix_Node *) calloc(1, 3 * SIZE_OF_POINTER);
    pMatNode->item = pMat;
    return pMatNode;
}

static void delete_node(Matrix_Node *pNode) {
    free(pNode->item->pointer_of_data);
    pNode->item->pointer_of_data = 0;
    free(pNode->item);
    pNode->item = 0;
    pNode->before->next = pNode->next;
    if (has_next(pNode)) {
        pNode->next->before = pNode->before;
    }
    free(pNode);
    pNode = 0;
    matrices_list.Matrix_counts--;
}

static void deleteAll(){
    for(int i = 0; i<matrices_list.Matrix_counts;){
        Matrix_Node *pNode = matrices_list.head->next;
        delete_node(pNode);
    }
}

//如果没找到返回NULL
//这个方法能够判断矩阵是不是由createMatrix创建的，如果不是会返回NULL
static Matrix_Node *find_node_by_matrix(const Matrix *pMat) {
    Matrix_Node *m = matrices_list.head;
    for (int i = 0; i < matrices_list.Matrix_counts; i++) {
        m = m->next;
        if (m->item == pMat) {
            return m;
        }
    }
    return 0;
}

static bool has_before(const Matrix_Node *pNode) {
    return pNode->before != 0;
}

static bool has_next(const Matrix_Node *pNode) {
    return pNode->next != 0;
}

static void reset(Matrix* pMat, size_t rows, size_t cols, float* pFloat){
    pMat->cols_numbers_in_Matrix = cols;
    pMat->rows_numbers_in_Matrix = rows;
    float* monitor = pMat->pointer_of_data;
    pMat->pointer_of_data = realloc(pMat->pointer_of_data, rows*cols*sizeof (float ));
    if(check_error("reset(Matrix* pMat, size_t rows, size_t cols, float* pFloat)", "Failure in resetting", pMat->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        pMat->pointer_of_data = monitor;
    }
    for(int i = 0; i<rows*cols; i++){
        pMat->pointer_of_data[i] = pFloat[i];
    }
}






/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part2:
 核心算法的具体实现：
 static方法，仅供内部调用，无需检查数据合法性
 能够更加简洁地看清楚不同矩阵功能的具体实现

///////////////////////////////////////////////////////////////////////////////////////////////////*/



//求行列式
static float det(const float * pFloat, size_t size){
    size_t rows = size;
    size_t cols = size;

    float cof[(size-1)*(size-1)];

    if(rows == 1){
        return pFloat[0];
    }
    else if(rows == 2){
        return pFloat[0]*pFloat[3] - pFloat[1]*pFloat[2];
    }
    else {
        float result = 0;
        for(size_t i = 0; i<size; i++){
            cofactor(pFloat, cof, size, 0, i);
            if(i%2 == 0){
                result += pFloat[i] * det(cof, size - 1);
            }
            else{
                result -= pFloat[i] * det(cof, size - 1);
            }
        }
        return result;
    }
}
//求余子式
static void cofactor(const float * pFloat,float * des, size_t size, size_t index_i, size_t index_j){
    //size是pFloat的边长
    for(int i = 0; i<size; i++){
        for(int j = 0; j<size; j++){
            if(i<index_i && j <index_j){
                des[i*(size-1)+j] = pFloat[i*size + j];
            }
            if(i<index_i && j > index_j){
                des[i*(size-1)+j-1] = pFloat[i*size + j];
            }
            if(i>index_i && j<index_j){
                des[(i-1)*(size-1)+j] = pFloat[i*size + j];
            }
            if(i>index_i && j>index_j){
                des[(i-1)*(size-1)+j-1] = pFloat[i*size + j];
            }
        }
    }
}

static void invert(const float *pFloat, float * des, size_t size){
    float deter = det(pFloat, size);
    adjoint_Matrix(pFloat, des, size);
    for(int i = 0; i<size; i++){
        for(int j = 0; j<size; j++){
            des[i*size+j] /= deter;
        }
    }
}

static void transpose(const float *pFloat, float *des, size_t rows, size_t cols){
    for(int i = 0; i<rows; i++){
        for(int j =0; j<cols; j++){
            des[j*rows + i]= pFloat[i*cols +j];
        }
    }
}

static void cofactor_Matrix(const float *pFloat, float *des, size_t sides){
    for(int i = 0; i<sides; i++){
        for(int j = 0; j<sides; j++){
            float mat[(sides-1)*(sides-1)];
            cofactor(pFloat, mat, sides, i, j);
            des[i*sides+j] = det(mat, sides-1)*(float)((i+j)%2 == 0?1:-1);
        }
    }
}

static void adjoint_Matrix(const float *pFloat, float *des, size_t sides){
    float mid[sides*sides];
    cofactor_Matrix(pFloat, mid, sides);
    transpose(mid, des,sides,sides);
}

static float trace(const float *pFloat, size_t sides){
    float result = 0;
    for(int i =0; i<sides; i++){
        result+= pFloat[i*sides+i];
    }
    return result;
}

static void Gause_elimination_to_upper_triangle_mat(float *augmented_pFloat, size_t rows, size_t cols){
    int col_index =0;
    for(int i = 0; i<rows; i++){
        //先判断等不等于0
        while(augmented_pFloat[i*cols+col_index] == 0){
            //先找到这个位置上不是0的行
            for(int k = i+1; k<rows;k++){
                if(augmented_pFloat[k*cols+col_index] != 0){
                    swap_rows(augmented_pFloat, cols, i, k);
                    goto eliminate;
                }
            }
            //如果找不到
            col_index++;
        }
        eliminate:
        for(int j = i+1; j<rows; j++){
            add_rows(augmented_pFloat, cols, i, j, -(augmented_pFloat[j*cols+col_index])/augmented_pFloat[i*cols+col_index]);
        }
        col_index++;
    }
}

static size_t triangle_mat_rank_for_aug(const float *simplified_augmented_pFloat, size_t rows, size_t cols){
    for(int i =rows-1; i>=0; i--){
        if(simplified_augmented_pFloat[i*cols+cols-1] != 0 || simplified_augmented_pFloat[i*cols+cols-2] != 0){
            return i+1;
        }
    }
    return rows;
}

static size_t triangle_mat_rank_for_cof(const float *simplified_augmented_pFloat, size_t rows, size_t cols){
    for(int i =rows-1; i>=0; i--){
        if(simplified_augmented_pFloat[i*cols+cols-2] != 0){
            return i+1;
        }
    }
    return rows;
}

static void printPF(const float*pF, size_t cols, size_t rows){
    for(int i = 0; i<rows; i++){
        for(int j = 0; j<cols; j++){
            printf("%f ",pF[i*cols+j]);
        }
        printf("\n");
    }
}
static void Gause_elimination_to_solve_equations(float *simplified_augmented_pFloat, size_t cols, size_t rank){
    for(int i = rank - 1; i>=0; i--){
        multiply_rows(simplified_augmented_pFloat, cols, i, 1.0f/simplified_augmented_pFloat[i*cols+i]);
        for(int j = i-1; j>=0; j--){
            add_rows(simplified_augmented_pFloat,cols,i,j,-(simplified_augmented_pFloat[j*cols+i])/simplified_augmented_pFloat[i*cols+i]);
        }
    }
}

///矩阵的初等行变换：\n
///1.行交换\n
///2.将i1行的coefficient倍加到第i2行\n
///3.行的倍乘\n
static void swap_rows(float *pFloat, size_t cols, size_t i1, size_t i2){
    for(int j = 0; j<cols; j++){
        float swap = pFloat[i1*cols+j];
        pFloat[i1*cols+j] = pFloat[i2*cols+j];
        pFloat[i2*cols+j] = swap;
    }
}

static void add_rows(float  *pFloat, size_t cols, size_t i1, size_t i2, float coefficient){
    for(int j = 0; j<cols; j++){
        pFloat[i2*cols+j] += coefficient*pFloat[i1*cols+j];
    }
}

static void multiply_rows(float *pFloat, size_t cols, size_t i1, float coefficient){
    for(int j = 0; j<cols; j++){
        pFloat[i1*cols+j] *= coefficient;
    }
}


static void add_two_mat(Matrix* des, const Matrix *pMat){
    for(int i = 0; i<des->rows_numbers_in_Matrix * des->cols_numbers_in_Matrix; i++){
        des->pointer_of_data[i] += pMat ->pointer_of_data[i];
    }
}

static void subtract_two_mat(Matrix* des, const Matrix *pMat){
    for(int i = 0; i<des->rows_numbers_in_Matrix * des->cols_numbers_in_Matrix; i++){
        des->pointer_of_data[i] -= pMat ->pointer_of_data[i];
    }
}


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part3:
 异常处理的定义


///////////////////////////////////////////////////////////////////////////////////////////////////*/


static bool check_error(const char *function_name, const char *description1, const void *errors_key_numbers, ERROR_ID id) {
    switch (id) {
        case ERROR_IN_ALLOCATING_MEMORY: {
            if(!errors_key_numbers){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s.\nError ID: %d\n",
                        description1, __FILE__, function_name, id);
                return false;
            }
            else {
                return true;
            }
        }
        case ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO: {
            if (*((long long *) errors_key_numbers) <= 0) {
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong Number: %lld\nError ID: %d\n",
                        description1, __FILE__, function_name, *((long long *) errors_key_numbers), id);
                return false;
            } else {
                return true;
            }
        }
        case ERROR_CAUSE_BY_NULL_POINTER: {
            if (errors_key_numbers == NULL) {
                fprintf(stdout, "%s \n",
                        description1);
                fprintf(stdout, "Error first occurs in: %s, in %s. Wrong address: %p\n", __FILE__, function_name,
                        errors_key_numbers);
                fprintf(stdout, "Error ID: %d\n", id);
                return false;
            } else {
                return true;
            }
        }
        case ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix: {
            Matrix_Node *matrixNode = find_node_by_matrix((Matrix *) errors_key_numbers);
            if (!matrixNode) {
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong address: %p\nError ID: %d\n",
                        description1, __FILE__, function_name, errors_key_numbers, id);
                return false;
            } else {
                return true;
            }
        }
        case ERROR_INDEX_OUT_OF_BOUND: {
            long long checked = ((long long*)errors_key_numbers)[0];
            long long standard = ((long long*)errors_key_numbers)[1];
            if(checked > standard || checked == 0) {
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong numbers: index = %lld\nError ID: %d\n",
                        description1, __FILE__, function_name, *(long long*)errors_key_numbers, id);
                return false;
            }
            else{
                return true;
            }
        }
        case ERROR_MODIFICATION_IN_MATRIX:{
            if(((Matrix*)errors_key_numbers)->rows_numbers_in_Matrix == 0 ){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong numbers: rows_numbers_in_Matrix = 0\nError ID: %d\n",
                        description1, __FILE__, function_name, id);
                return false;
            }
            else if(((Matrix*)errors_key_numbers)->cols_numbers_in_Matrix == 0){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong numbers: cols_numbers_in_Matrix = 0\nError ID: %d\n",
                        description1, __FILE__, function_name, id);
                return false;
            }
            else if(((Matrix*)errors_key_numbers)->pointer_of_data == 0){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong numbers: pointer_of_data = NULL\nError ID: %d\n",
                        description1, __FILE__, function_name, id);
                return false;
            }
            else {
                return true;
            }
        }
        case ERROR_NOT_THE_SAME_SIZE: {
            Matrix * m1 = ((Matrix**)errors_key_numbers)[0];
            Matrix * m2 = ((Matrix**)errors_key_numbers)[1];
            if(m1 -> rows_numbers_in_Matrix != m2->rows_numbers_in_Matrix){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong numbers: rows of Matrix1 = %zu, rows of Matrix2 = %zu\nError ID: %d\n",
                        description1, __FILE__, function_name, m1 -> rows_numbers_in_Matrix, m2->rows_numbers_in_Matrix, id);
                return false;
            }
            else if(m1 -> cols_numbers_in_Matrix != m2->cols_numbers_in_Matrix){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong numbers: cols of Matrix1 = %zu, cols of Matrix2 = %zu\nError ID: %d\n",
                        description1, __FILE__, function_name, m1 -> cols_numbers_in_Matrix, m2->cols_numbers_in_Matrix, id);
                return false;
            }
            else {
                return true;
            }
        }
        case ERROR_INVALID_COLS_AND_ROWS_IN_MULTIPLICATION:{
            Matrix * m1 = ((Matrix**)errors_key_numbers)[0];
            Matrix * m2 = ((Matrix**)errors_key_numbers)[1];
            if(m1 -> cols_numbers_in_Matrix != m2->rows_numbers_in_Matrix){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong numbers: cols of Matrix1 = %zu, rows of Matrix2 = %zu\nError ID: %d\n",
                        description1, __FILE__, function_name, m1 -> cols_numbers_in_Matrix, m2->rows_numbers_in_Matrix, id);
                return false;
            }
            else{
                return true;
            }
        }
        case ERROR_INVALID_MATRIX_ELEMENT_FROM_KEYBOARD: {
            //只负责检查不释放内存
            char* string = (char*)errors_key_numbers;
            char *p = string;
            while(*p){
                float test_float = strtof(string, &p);
                if(string == p && *string != '\0'){
                    while(*p == ' ' || *p == '\n'){
                        p++;
                    }
                    if(*p == '\0'){
                        return true;
                    }
                    //printf("here\n");
                    fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong sites: %c\nError ID: %d\n",
                            description1, __FILE__, function_name, (int)*p, id);
                    return false;
                }
                else if(string == p && *string == '\0')
                {
                    return true;
                }
                else{
                    string = p;
                }
            }
            return true;
        }
        case ERROR_INVALID_NUMBER_OF_LINES_FROM_KEYBOARD: {
            //只有一个数且这个数大于0
            char* string = (char*)errors_key_numbers;
            char *p = string;
            size_t str_size = strlen(string);
            long test_long = strtol(string, &p,10);
            char* terminal = &string[str_size];  //这个位置是\0
            while(p != terminal){
                //检查第一个数后面有没有非空格的字符
                if(*p != ' '){
                    fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong sites: %c\nError ID: %d\n",
                            description1, __FILE__, function_name, *p, id);
                    return false;
                }
                p++;
            }
            if(test_long<=0){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong sites: %ld\nError ID: %d\n",
                        description1, __FILE__, function_name, test_long, id);
                return false;
            }
            return true;
        }
        case ERROR_FAIL_TO_OPEN_FILE: {
            if (errors_key_numbers == NULL) {
                fprintf(stdout, "%s \n",
                        description1);
                fprintf(stdout, "Error first occurs in: %s, in %s. Wrong file pointer address: %p\n", __FILE__, function_name,
                        errors_key_numbers);
                fprintf(stdout, "Error ID: %d\n", id);
                return false;
            } else {
                return true;
            }
        }
        case ERROR_INVALID_MATRIX_ELEMENT_FROM_FILE:{
            char* string = (char*)errors_key_numbers;
            char *p = string;
            while(*p){
                float test_float = strtof(string, &p);
                if(string == p && *string != '\0'){
                    while(*p == ' ' || *p == '\n'){
                        p++;
                    }
                    if(*p == '\0'){
                        return true;
                    }
                    fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong sites: %c\nError ID: %d\n",
                            description1, __FILE__, function_name, *p, id);
                    return false;
                }
                else if(string == p && *string == '\0')
                {
                    return true;
                }
                else{
                    string = p;
                }
            }
            return true;
        }
        case ERROR_UNSUPPORTED_MODE: {
            if((!strcmp((char*)errors_key_numbers, "w")) || (!strcmp((char*)errors_key_numbers, "a"))){
                return true;
            }
            else{
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong sites: %s\nError ID: %d\n",
                        description1, __FILE__, function_name, (char*)errors_key_numbers, id);
                return false;
            }
        }
        case ERROR_NOT_A_SQUARE_MATRIX: {
            Matrix * mat = (Matrix*)errors_key_numbers;
            if(mat->rows_numbers_in_Matrix != mat->cols_numbers_in_Matrix){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. Wrong rows: %zu, cols: %zu\nError ID: %d\n",
                        description1, __FILE__, function_name, mat->rows_numbers_in_Matrix, mat->cols_numbers_in_Matrix, id);
                return false;
            }
            else{
                return true;
            }
        }
        case ERROR_MATRIX_IRREVERSIBLE: {
            Matrix * mat = (Matrix*)errors_key_numbers;
            float de = det(mat->pointer_of_data, mat->rows_numbers_in_Matrix);
            if(de == 0){
                fprintf(stdout, "%s \nError first occurs in: %s, in %s. irreversible matrix address: %p\nError ID: %d\n",
                        description1, __FILE__, function_name, mat, id);
                return false;
            }
            else{
                return true;
            }
        }
        default: {
            printf("unexpected error occur!\n");
            return false;
        }
    }
}





/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part4:
 声明矩阵的创建方法，初始化方法，设置方法，删除方法，复制方法，打印方法等等

///////////////////////////////////////////////////////////////////////////////////////////////////*/



//在这里使用long long作为参数类型，是为了在用户输入负数时提醒用户，而不是建一个规模很大的矩阵
// 默认初始化，每个位置的值随机
Matrix *createMatrix( long long rows, long long cols) {
    bool valid = check_error("createMatrix(long long rows, long long cols)", "Failure in creation!",
                             &rows,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO) &&
                 check_error("createMatrix(long long rows, long long cols)", "Failure in creation!",
                             &cols,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO);
    if (!valid) {
        return 0;
    }
    Matrix *m = (Matrix *) malloc(sizeof(size_t) + sizeof(size_t)+ SIZE_OF_POINTER); //避免不同操作系统的字节数不一样
    m->rows_numbers_in_Matrix = rows;
    m->cols_numbers_in_Matrix = cols;
    m->pointer_of_data = (float *) malloc(rows * cols * sizeof(float));
    if(!check_error("createMatrix(long long rows, long long cols)","failure in allocating memory",m->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        return NULL;
    }
    add_node_to_list(create_node(m));
    return m;
}

Matrix *createMatrix_with_pFloat( long long rows, long long cols, const float *pFloat, long long arr_size){
    bool valid = check_error("createMatrix_with_pFloat( long long rows, long long cols, const float *pFloat, long long arr_size)", "Failure in creation!",
                             &rows,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO) &&
                 check_error("createMatrix_with_pFloat( long long rows, long long cols, const float *pFloat, long long arr_size)", "Failure in creation!",
                             &cols,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO);
    if (!valid) {
        return 0;
    }
    Matrix *m = (Matrix *) malloc(sizeof(size_t) + sizeof(size_t)+ SIZE_OF_POINTER); //避免不同操作系统的字节数不一样
    m->rows_numbers_in_Matrix = rows;
    m->cols_numbers_in_Matrix = cols;
    m->pointer_of_data = (float *) malloc(rows * cols * sizeof(float));
    if(!check_error("createMatrix_with_pFloat( long long rows, long long cols, const float *pFloat, long long arr_size)","failure in allocating memory",
                    m->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        return NULL;
    }
    add_node_to_list(create_node(m));
    int i = 0;
    for (; i < min((m->rows_numbers_in_Matrix) * (m->cols_numbers_in_Matrix), arr_size); i++) {
        m->pointer_of_data[i] = pFloat[i];
    }
    while(i < max((m->rows_numbers_in_Matrix) * (m->cols_numbers_in_Matrix), arr_size)){
        m->pointer_of_data[i] = 0;
        i++;
    }
    return m;
}

Matrix *createMatrix_from_keyboard(long long rows, long long cols){
    bool valid = check_error("createMatrix(long long rows, long long cols)", "Failure in creation!",
                             &rows,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO) &&
                 check_error("createMatrix(long long rows, long long cols)", "Failure in creation!",
                             &cols,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO);
    if (!valid) {
        return 0;
    }
    Matrix *pMat = (Matrix *) malloc(sizeof(size_t) + sizeof(size_t)+ SIZE_OF_POINTER); //避免不同操作系统的字节数不一样
    pMat->rows_numbers_in_Matrix = rows;
    pMat->cols_numbers_in_Matrix = cols;
    pMat->pointer_of_data = (float *) malloc(rows * cols * sizeof(float));
    if(!check_error("createMatrix(long long rows, long long cols)","failure in allocating memory",pMat->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        return NULL;
    }
    add_node_to_list(create_node(pMat));
    set_Matrix_from_keyboard(pMat);
    return pMat;
}

Matrix* createMatrix_with_a_file(long long rows, long long cols, const char* filename_with_directory){
    bool valid = check_error("createMatrix(long long rows, long long cols)", "Failure in creation!",
                             &rows,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO) &&
                 check_error("createMatrix(long long rows, long long cols)", "Failure in creation!",
                             &cols,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO);
    if (!valid) {
        return 0;
    }
    Matrix *pMat = (Matrix *) malloc(sizeof(size_t) + sizeof(size_t)+ SIZE_OF_POINTER); //避免不同操作系统的字节数不一样
    pMat->rows_numbers_in_Matrix = rows;
    pMat->cols_numbers_in_Matrix = cols;
    pMat->pointer_of_data = (float *) malloc(rows * cols * sizeof(float));
    if(!check_error("createMatrix(long long rows, long long cols)","failure in allocating memory",pMat->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        return NULL;
    }
    add_node_to_list(create_node(pMat));
    set_Matrix_from_a_file(pMat,filename_with_directory);
    return pMat;
}

//返回一个单位矩阵的指针
Matrix *create_an_identity_matrix(long long rows) {
    bool valid = check_error("create_an_identity_matrix(long long rows)", "Failure in creation!", &rows,
                             ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO);
    if (!valid) {
        return 0;
    }
    Matrix *m = (Matrix *) malloc(sizeof(size_t) + sizeof(size_t)+ SIZE_OF_POINTER);
    m->rows_numbers_in_Matrix = rows;
    m->cols_numbers_in_Matrix = rows;
    m->pointer_of_data = (float *) malloc(rows * rows * sizeof(float));
    if(!check_error("create_an_identity_matrix(long long rows)","failure in allocating memory",m->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        return NULL;
    }
    memset(m->pointer_of_data, 0, rows * rows * sizeof(float));
    for (int i = 0; i < rows; i++) {
        m->pointer_of_data[i * rows + i] = 1.0f;
    }
    add_node_to_list(create_node(m));
    return m;
}

//返回一个零矩阵的指针
Matrix *create_an_zero_matrix(long long rows, long long cols) {
    bool valid = check_error("create_an_zero_matrix(long long rows, long long cols)", "Failure in creation!", &rows,
                             ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO) &&
                 check_error("create_an_zero_matrix(long long rows, long long cols)", "Failure in creation!", &cols,
                             ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO);
    if (!valid) {
        return 0;
    }
    Matrix *m = (Matrix *) malloc(sizeof(size_t) + sizeof(size_t)+ SIZE_OF_POINTER);
    m->rows_numbers_in_Matrix = rows;
    m->cols_numbers_in_Matrix = cols;
    m->pointer_of_data = (float *) malloc(rows * cols * sizeof(float));
    if(!check_error("create_an_zero_matrix(long long rows, long long cols) ","failure in allocating memory",m->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        return NULL;
    }
    memset(m->pointer_of_data, 0, rows * cols * sizeof(float));
    add_node_to_list(create_node(m));
    return m;
}

Matrix* create_an_one_matrix(long long rows, long long cols){
    bool valid = check_error("create_an_zero_matrix(long long rows, long long cols)", "Failure in creation!", &rows,
                             ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO) &&
                 check_error("create_an_zero_matrix(long long rows, long long cols)", "Failure in creation!", &cols,
                             ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO);
    if (!valid) {
        return 0;
    }
    Matrix *m = (Matrix *) malloc(sizeof(size_t) + sizeof(size_t)+ SIZE_OF_POINTER);
    m->rows_numbers_in_Matrix = rows;
    m->cols_numbers_in_Matrix = cols;
    m->pointer_of_data = (float *) malloc(rows * cols * sizeof(float));
    if(!check_error("create_an_zero_matrix(long long rows, long long cols) ","failure in allocating memory",m->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
        return NULL;
    }
    memset(m->pointer_of_data, 0, rows * cols * sizeof(float));
    for(int i = 0; i<m->cols_numbers_in_Matrix*m->rows_numbers_in_Matrix; i++){
        m->pointer_of_data[i] = 1;
    }
    add_node_to_list(create_node(m));
    return m;
}

//初始化矩阵
void set_Matrix_same_value(Matrix *matrix, float same_value) {
    bool valid = check_error("init_Matrix_same_value(Matrix *matrix, float same_value)", "Failure in initialization!", matrix,
                             ERROR_CAUSE_BY_NULL_POINTER);
    if (!valid) {
        return;
    }
    bool valid2 = check_error("init_Matrix_same_value(Matrix *matrix, float same_value)", "Failure in initialization!", matrix,
                              ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix);
    if (!valid2) {
        return;
    }
    for (int i = 0; i < (matrix->rows_numbers_in_Matrix) * (matrix->cols_numbers_in_Matrix); i++) {
        matrix->pointer_of_data[i] = same_value;
    }
}

void set_Matrix_with_an_array(Matrix *matrix, const float *pFloat, long long arr_size) {
    bool valid = check_error("init_Matrix_with_an_array(Matrix *matrix, const float *pFloat, long long arr_size) ", "Failure in initialization!",
                             &arr_size,ERROR_MATRIX_SIZE_LESS_OR_EQUAL_TO_ZERO) &&
                 check_error("init_Matrix_with_an_array(Matrix *matrix, const float *pFloat, long long arr_size) ", "Failure in initialization!",
                             matrix,ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("init_Matrix_with_an_array(Matrix *matrix, const float *pFloat, long long arr_size) ", "Failure in initialization!",
                             pFloat,ERROR_CAUSE_BY_NULL_POINTER);
    if (!valid) {
        return;
    }
    bool valid2 = check_error("init_Matrix_with_an_array()", "Failure in initialization!", matrix,
                              ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix);
    if (!valid2) {
        return;
    }
    for (int i = 0; i < min((matrix->rows_numbers_in_Matrix) * (matrix->cols_numbers_in_Matrix), arr_size); i++) {
        matrix->pointer_of_data[i] = pFloat[i];
    }
}

//用户从键盘输入，如果输入的数超过矩阵规模，后面的数会被舍弃，如果小于矩阵规模，矩阵剩下的数会自动补0；
void set_Matrix_from_keyboard(Matrix *pMat){
    bool valid = check_error("set_Matrix_from_keyboard(Matrix *pMat)", "Failure in setting!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("set_Matrix_from_keyboard(Matrix *pMat)", "Failure in setting!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("set_Matrix_from_keyboard(Matrix *pMat)", "Failure in setting!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    //先看看用户想输入多少行：
    printf("How many lines you want to input?\n");
    size_t lines_size = 0;
    char * lines_num = NULL;
    while(true){
        char c= getchar();
        if(c == '\n'){
            break;
        }
        else{
            if(lines_num == NULL){
                lines_num = malloc((++lines_size) * sizeof(char));
                if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","failure in allocating memory",
                                lines_num,ERROR_IN_ALLOCATING_MEMORY)){
                    return;
                }
                lines_num[lines_size-1] = c;
            }
            else{
                char *monitor_pointer = lines_num;
                lines_num = realloc(lines_num, ++lines_size);
                if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","failure in allocating memory",
                                lines_num,ERROR_IN_ALLOCATING_MEMORY)){
                    free(monitor_pointer); //避免因为realloc失败导致内存泄漏；
                    return;
                }
                lines_num[lines_size-1] = c;
            }
        }
    }
    //添上终止符：
    char *monitor_pointer = lines_num;
    lines_num = realloc(lines_num, (lines_size+=2));
    if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","failure in allocating memory",lines_num,ERROR_IN_ALLOCATING_MEMORY)){
        free(monitor_pointer); //避免因为realloc失败导致内存泄漏；
        return;
    }
    lines_num[lines_size-2] = '\0';
    lines_num[lines_size-1] = 'A';
    if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","Invalid input!",
                    lines_num,ERROR_INVALID_NUMBER_OF_LINES_FROM_KEYBOARD)){
        free(lines_num);
        return;
    }
    char *arg_p;
    long lines = strtol(lines_num, &arg_p, 10);
    free(lines_num);


    //接下来让用户输入字符串
    printf("Please input %ld lines:\n",lines);
    float * pFloat = pMat->pointer_of_data;
    size_t signal = 0;
    for(int i = 0; i<lines; i++){
        size_t str_size = 0;
        char *string = NULL;
        while(true){
            char c= getchar();
            if(c == '\n'){
                break;
            }
            else{
                if(string == NULL){
                    string = malloc((++str_size) * sizeof(char));
                    if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","failure in allocating memory",
                                    string,ERROR_IN_ALLOCATING_MEMORY)){
                        return;
                    }
                    string[str_size-1] = c;
                }
                else{
                    char *monitor_pointer2 = string;
                    string = realloc(string, (++str_size)* sizeof(char));
                    if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","failure in allocating memory",
                                    string,ERROR_IN_ALLOCATING_MEMORY)){
                        free(monitor_pointer2); //避免因为realloc失败导致内存泄漏；
                        return;
                    }
                    string[str_size-1] = c;
                }
            }
        }
        //添上终止符：
        char *monitor_pointer2 = string;
        string = realloc(string, (str_size+=2)* sizeof(char));
        if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","failure in allocating memory",
                        string,ERROR_IN_ALLOCATING_MEMORY)){
            free(monitor_pointer2); //避免因为realloc失败导致内存泄漏；
            return;
        }
        string[str_size-2] = '\0';
        string[str_size-1] = 'A';
        if(!check_error("set_Matrix_from_keyboard(Matrix *pMat)","Invalid input!",
                        string,ERROR_INVALID_MATRIX_ELEMENT_FROM_KEYBOARD)){
            free(string);
            return;
        }
        char *p_front = string;
        char *p_back = string;
        while(*p_front){
            if(signal >= pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix){
                free(string);
                return;
            }
            float test_float = strtof(p_front, &p_back);
            if(p_front == p_back){
                break;
            }
            else{
                p_front = p_back;
            }
            *pFloat = test_float;
            signal++;
            pFloat++;
        }
        free(string);
    }
    while(signal++<pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix){
        *pFloat = 0.f;
        pFloat++;
    }
}

void set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory){
    bool valid = check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)", "Failure in setting!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)", "Failure in setting!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)", "Failure in setting!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    FILE * file = fopen(filename_with_directory, "r");
    bool valid2 = check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)", "Failure in setting!",
                              file, ERROR_FAIL_TO_OPEN_FILE);
    if(!valid2){
        fclose(file);
        return;
    }
    float * pFloat = pMat->pointer_of_data;
    size_t signal = 0;
    size_t str_size = 0;
    char *string = NULL;
    char c;
    while( (c= fgetc(file)) != EOF) {
        //!feof(file)
        if (string == NULL) {
            string = malloc((++str_size) * sizeof(char));
            if (!check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)", "failure in allocating memory",
                             string, ERROR_IN_ALLOCATING_MEMORY)) {
                fclose(file);
                return;
            }
            string[str_size - 1] = c;
        } else {
            char *monitor_pointer2 = string;
            string = realloc(string, (++str_size)* sizeof(char));
            if (!check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)", "failure in allocating memory",
                             string, ERROR_IN_ALLOCATING_MEMORY)) {
                free(monitor_pointer2); //避免因为realloc失败导致内存泄漏；
                fclose(file);
                return;
            }
            string[str_size - 1] = c;
        }
    }
    //添上终止符：
    char *monitor_pointer2 = string;
    string = realloc(string, (str_size+=2) * sizeof(char));
    if(!check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)","failure in allocating memory",string,ERROR_IN_ALLOCATING_MEMORY)){
        free(monitor_pointer2); //避免因为realloc失败导致内存泄漏；
        fclose(file);
        return;
    }
    string[str_size-2] = '\0';
    string[str_size-1] = 'A';
    if(!check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)","Invalid character in file!",
                    string,ERROR_INVALID_MATRIX_ELEMENT_FROM_FILE)){
        free(string);
        fclose(file);
        return;
    }
    char *p_front = string;
    char *p_back = string;
    while(*p_front){
        if(signal >= pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix){
            free(string);
            fclose(file);
            return;
        }
        float test_float = strtof(p_front, &p_back);
        if(p_front == p_back){
            break;
        }
        else{
            p_front = p_back;
        }
        *pFloat = test_float;
        signal++;
        pFloat++;
    }
    free(string);

    while(signal++<pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix){
        *pFloat = 0.f;
        pFloat++;
    }
    fclose(file);
}

//将第row行第col列的数设为相应的数值
void set_Matrix(Matrix *matrix, long long row, long long col, float value){
    bool valid = check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                             matrix, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                             matrix, ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix);
    if(!valid){
        return;
    }
    bool valid2 = check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                              matrix, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid2){
        return;
    }
    long long rows [2] = {row, matrix->rows_numbers_in_Matrix};
    long long cols [2] = {col, matrix->cols_numbers_in_Matrix};
    bool valid3 = check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                              rows, ERROR_INDEX_OUT_OF_BOUND) &&
                  check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                              cols, ERROR_INDEX_OUT_OF_BOUND);
    if(!valid3){
        return;
    }
    matrix->pointer_of_data[(row - 1) * matrix->cols_numbers_in_Matrix + (col - 1)] = value;
}

void set_whole_Matrix_self_define(Matrix *matrix, float (*self_define_function)(float float_element)){
    bool valid = check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                             matrix, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                             matrix, ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("set_Matrix(Matrix *matrix, int row, int col, float value)", "Failure in setting!",
                             matrix, ERROR_MODIFICATION_IN_MATRIX);;
    if(!valid){
        return;
    }
    for(int i =0; i<matrix->rows_numbers_in_Matrix*matrix->cols_numbers_in_Matrix; i++){
        matrix->pointer_of_data[i] = self_define_function(matrix->pointer_of_data[i]);
    }
}

void deleteMatrix(Matrix **p_pMat) {
    bool valid = check_error("deleteMatrix(Matrix **p_pMat)", "Failure in destruction!", p_pMat, ERROR_CAUSE_BY_NULL_POINTER)&&
                 check_error("deleteMatrix(Matrix **p_pMat)", "Failure in destruction!", *p_pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("deleteMatrix(Matrix **p_pMat)", "Failure in destruction!", *p_pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("deleteMatrix(Matrix **p_pMat)", "Failure in destruction!", *p_pMat, ERROR_MODIFICATION_IN_MATRIX);
    if (!valid) {
        return;
    }
    Matrix_Node *matrixNode = find_node_by_matrix(*p_pMat);
    delete_node(matrixNode);
    *p_pMat = 0;
}

void delete_all_Matrix(){
    deleteAll();
}

void copyMatrix(const Matrix* pMat, Matrix* des){
    bool valid = check_error("copyMatrix(const Matrix* pMat, Matrix* des)", "Failure in copy!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("copyMatrix(const Matrix* pMat, Matrix* des)", "Failure in copy!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("copyMatrix(const Matrix* pMat, Matrix* des)", "Failure in copy!", pMat, ERROR_MODIFICATION_IN_MATRIX) &&
                 check_error("copyMatrix(const Matrix* pMat, Matrix* des)", "Failure in copy!", des, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("copyMatrix(const Matrix* pMat, Matrix* des)", "Failure in copy!", des,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("copyMatrix(const Matrix* pMat, Matrix* des)", "Failure in copy!", des, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }

    des->rows_numbers_in_Matrix = pMat->rows_numbers_in_Matrix;
    des->cols_numbers_in_Matrix = pMat->cols_numbers_in_Matrix;
    float *monitor_pointer2 = des->pointer_of_data;
    des->pointer_of_data = realloc(des->pointer_of_data, des->cols_numbers_in_Matrix*des->rows_numbers_in_Matrix*sizeof(float ));
    if(!check_error("set_Matrix_from_a_file(Matrix *pMat, const char* filename_with_directory)","failure in allocating memory",des->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)) {
        free(monitor_pointer2); //避免因为realloc失败导致内存泄漏；
        return;
    }
    for(int i  =0; i<pMat->rows_numbers_in_Matrix* pMat->cols_numbers_in_Matrix; i++){
        des->pointer_of_data[i] = pMat ->pointer_of_data[i];
    }
}

Matrix *copyMatrix_and_return(const Matrix* pMat){
    bool valid = check_error("copyMatrix_and_return(const Matrix* pMat)", "Failure in copy!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("copyMatrix_and_return(const Matrix* pMat)", "Failure in copy!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("copyMatrix_and_return(const Matrix* pMat)", "Failure in copy!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix *des = createMatrix(pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    for(int i  =0; i<pMat->rows_numbers_in_Matrix* pMat->cols_numbers_in_Matrix; i++){
        des->pointer_of_data[i] = pMat ->pointer_of_data[i];
    }
    return des;
}

void printMatrix(const Matrix* pMat){
    bool valid = check_error("printMatrix(const Matrix* pMat)", "Failure in printing!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("printMatrix(const Matrix* pMat)", "Failure in printing!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("printMatrix(const Matrix* pMat)", "Failure in printing!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    for(int i = 0; i<pMat->rows_numbers_in_Matrix; i++){
        printf("[  ");
        for(int j = 0; j<pMat->cols_numbers_in_Matrix; j++){
            printf("%f  ", pMat->pointer_of_data[i*pMat->cols_numbers_in_Matrix+j]);
        }
        printf("]\n");
    }
}

void printMatrix_to_a_file(const Matrix* pMat, const char* filename_with_directory, const char* mode){
    bool valid = check_error("printMatrix_to_a_file(const Matrix* pMat, const char* filename_with_directory, const char* mode)", "Failure in printing!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("printMatrix_to_a_file(const Matrix* pMat, const char* filename_with_directory, const char* mode)", "Failure in printing!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("printMatrix_to_a_file(const Matrix* pMat, const char* filename_with_directory, const char* mode)", "Failure in printing!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    bool valid1 = check_error("printMatrix_to_a_file(const Matrix* pMat, const char* filename_with_directory, const char* mode)", "Failure in printing!",
                              mode, ERROR_UNSUPPORTED_MODE);
    if(!valid1){
        return;
    }
    FILE *file = fopen(filename_with_directory, mode);
    bool valid2 = check_error("printMatrix_to_a_file(const Matrix* pMat, const char* filename_with_directory, const char* mode)", "Failure in printing!",
                              file, ERROR_FAIL_TO_OPEN_FILE);
    if(!valid2){
        return;
    }
    for(int i = 0; i<pMat->rows_numbers_in_Matrix; i++){
        fprintf(file,"[  ");
        for(int j = 0; j<pMat->cols_numbers_in_Matrix; j++){
            fprintf(file,"%f  ", pMat->pointer_of_data[i*pMat->cols_numbers_in_Matrix+j]);
        }
        fprintf(file,"]\n");
    }
}

void print_matrices_list() {
    Matrix_Node *m = matrices_list.head;
    if(matrices_list.Matrix_counts == 0){
        printf("No matrix left.\n");
        return;
    }
    for (int i = 0; i < matrices_list.Matrix_counts; i++) {
        m = m->next;
        printf("Matrix %d: %zu rows, %zu columns\n", i, m->item->rows_numbers_in_Matrix,
               m->item->cols_numbers_in_Matrix);
    }
}

/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part5:
 求矩阵的一些特性，不会对矩阵做出任何的改动

///////////////////////////////////////////////////////////////////////////////////////////////////*/

inline float minimum_in_a_matrix(const Matrix *pMat) {
    bool valid = check_error("deleteMatrix()", "Failure in destruction!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("deleteMatrix()", "Failure in destruction!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("deleteMatrix()", "Failure in destruction!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0;
    }
    float min = pMat->pointer_of_data[0];
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix -1; i++){
        min = min(pMat->pointer_of_data[i], min);
    }
    return min;
}

inline float maximum_in_a_matrix(const Matrix *pMat) {
    bool valid = check_error("deleteMatrix()", "Failure in destruction!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("deleteMatrix()", "Failure in destruction!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("deleteMatrix()", "Failure in destruction!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0;
    }
    float max = pMat->pointer_of_data[0];
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix -1; i++){
        max = max(pMat->pointer_of_data[i], max);
    }
    return max;
}

float determinant(const Matrix* pMat){
    bool valid = check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0.f;
    }
    bool valid1 = check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                              pMat, ERROR_NOT_A_SQUARE_MATRIX);
    if(!valid1){
        return 0.f;
    }
    float result = det(pMat->pointer_of_data, pMat->rows_numbers_in_Matrix);
    return result;
}

float trace_matrix(const Matrix *pMat){
    bool valid = check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0.f;
    }
    bool valid1 = check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                              pMat, ERROR_NOT_A_SQUARE_MATRIX);
    if(!valid1){
        return 0.f;
    }
    float result = trace(pMat->pointer_of_data, pMat->rows_numbers_in_Matrix);
    return result;
}
size_t rank(const Matrix *pMat){
    bool valid = check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("determinant(const Matrix* pMat)", "Failure in solving determinant!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0;
    }
    float pF[pMat->rows_numbers_in_Matrix*pMat->cols_numbers_in_Matrix];
    for(int i = 0; i<pMat->rows_numbers_in_Matrix*pMat->cols_numbers_in_Matrix; i++){
        pF[i] = pMat->pointer_of_data[i];
    }
    Gause_elimination_to_upper_triangle_mat(pF, pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    size_t result = triangle_mat_rank_for_aug(pF,pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    return result;
}


/*///////////////////////////////////////////////////////////////////////////////////////////////////
 part6:
 进行一些矩阵的运算，大部分运算都提供了两种方法：
 一种是动态申请一个新的矩阵储存结果并返回（上半部分）
 一种是对原有的某一个参数中的矩阵进行修改（下半部分）

///////////////////////////////////////////////////////////////////////////////////////////////////*/



Matrix* addConstMatrices(const Matrix* pMat, ...) {
    bool valid = check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                           check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                           check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->rows_numbers_in_Matrix,pMat->cols_numbers_in_Matrix);
    copyMatrix(pMat, result);
    va_list ap;
    va_start(ap, pMat);
    while(true){
        const Matrix * p = va_arg(ap, Matrix*);
        if(!find_node_by_matrix(p)){
            va_end(ap);
            return result;
        }
        const Matrix * mats[2] = {result, p};
        bool continue_or_not = check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!",
                                           p, ERROR_CAUSE_BY_NULL_POINTER) &&
                                            check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!",
                                                        p, ERROR_MODIFICATION_IN_MATRIX) &&
                                            check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!",
                                                        mats, ERROR_NOT_THE_SAME_SIZE) ;
        if(!continue_or_not ){
            va_end(ap);
            return result;
        }
        else {
            add_two_mat(result, p);
        }
    }
    va_end(ap);
    return result;
}

Matrix* subtractConstMatrices(const Matrix* pMat, ...) {
    bool valid = check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->rows_numbers_in_Matrix,pMat->cols_numbers_in_Matrix);
    copyMatrix(pMat, result);
    va_list ap;
    va_start(ap, pMat);
    while(true){
        const Matrix * p = va_arg(ap, Matrix*);
        if(!find_node_by_matrix(p)){
            va_end(ap);
            return result;
        }
        const Matrix * mats[2] = {result, p};
        bool continue_or_not = check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                                           p, ERROR_CAUSE_BY_NULL_POINTER) &&
                               check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                                           p, ERROR_MODIFICATION_IN_MATRIX) &&
                               check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                                           mats, ERROR_NOT_THE_SAME_SIZE) ;
        if(!continue_or_not ){
            va_end(ap);
            return result;
        }
        else {
            subtract_two_mat(result, p);
        }
    }
    va_end(ap);
    return result;
}

Matrix* MultiplyConstMatrices(const Matrix* pMat, ...) {
    bool valid = check_error("MultiplyMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in multiplication!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("MultiplyMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in multiplication!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("MultiplyMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in multiplication!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix * first_mat = createMatrix(pMat->rows_numbers_in_Matrix,pMat->cols_numbers_in_Matrix);
    copyMatrix(pMat, first_mat);
    Matrix * mid = first_mat;
    va_list ap;
    va_start(ap, pMat);
    while(true){
        const Matrix * p = va_arg(ap, Matrix*);
        if(!find_node_by_matrix(p)){
            float *pFloat = mid->pointer_of_data;
            size_t final_rows = mid->rows_numbers_in_Matrix;
            size_t final_cols = mid->cols_numbers_in_Matrix;
            mid = createMatrix(final_rows, final_rows);
            set_Matrix_with_an_array(mid, pFloat, final_rows*final_cols);
            deleteMatrix(&first_mat);
            va_end(ap);
            return mid;
        }
        const Matrix * mats[2] = {mid, p};
        bool continue_or_not = check_error("MultiplyMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in multiplication!",
                                           p, ERROR_CAUSE_BY_NULL_POINTER) &&
                               check_error("MultiplyMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in multiplication!",
                                           p, ERROR_MODIFICATION_IN_MATRIX) &&
                               check_error("MultiplyMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in multiplication!",
                                           mats, ERROR_INVALID_COLS_AND_ROWS_IN_MULTIPLICATION) ;
        if(!continue_or_not){
            float *pFloat = mid->pointer_of_data;
            size_t final_rows = mid->rows_numbers_in_Matrix;
            size_t final_cols = mid->cols_numbers_in_Matrix;
            mid = createMatrix(final_rows, final_rows);
            set_Matrix_with_an_array(mid, pFloat, final_rows*final_cols);
            deleteMatrix(&first_mat);
            va_end(ap);
            return mid;
        }
        else {
            size_t times = mid->cols_numbers_in_Matrix;
            float ans[mid->rows_numbers_in_Matrix*p->cols_numbers_in_Matrix];
            for(int i = 0; i<mid->rows_numbers_in_Matrix; i++){
                for(int j = 0; j<p->cols_numbers_in_Matrix; j++){
                    ans[i*p->cols_numbers_in_Matrix+j] = 0;
                    for(int k = 0; k<times; k++){
                        ans[i*p->cols_numbers_in_Matrix+j] += mid->pointer_of_data[i*mid->cols_numbers_in_Matrix+k] *
                                p->pointer_of_data[k*mid->cols_numbers_in_Matrix+j];
                        //printf("%f \n",ans[i*p->cols_numbers_in_Matrix+j]);
                    }
                }
            }
            free(mid->pointer_of_data);
            mid->pointer_of_data = (float *) calloc(mid->rows_numbers_in_Matrix*p->cols_numbers_in_Matrix, sizeof(float));
            if(!check_error("MultiplyMatrices(int MaxSize, const Matrix* pMat, ...)","failure in allocating memory",mid->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
                return NULL;
            }
            for(int i = 0; i<mid->rows_numbers_in_Matrix*p->cols_numbers_in_Matrix;i++){
                mid->pointer_of_data[i] = ans[i];
            }
            mid->cols_numbers_in_Matrix = p->cols_numbers_in_Matrix;
        }
    }
    va_end(ap);
    return first_mat;
}

Matrix *add_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar) {
    bool valid = check_error("add_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in addition!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("add_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in addition!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("add_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in addition!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0;
    }
    Matrix * new_mat = createMatrix(pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix -1; i++){
        new_mat->pointer_of_data[i] = pMat->pointer_of_data[i] + scalar;
    }
    return new_mat;
}

Matrix * subtract_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar) {
    bool valid = check_error("subtract_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in subtraction!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("subtract_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in subtraction!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("subtract_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in subtraction!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0;
    }
    Matrix * new_mat = createMatrix(pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix -1; i++){
        new_mat->pointer_of_data[i] = pMat->pointer_of_data[i] - scalar;
    }
    return new_mat;
}

Matrix * multiply_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar) {
    bool valid = check_error("multiply_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in multiplication!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("multiply_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in multiplication!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("multiply_a_scalar_to_a_const_Matrix(const Matrix *pMat, float scalar)", "Failure in multiplication!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return 0;
    }
    Matrix * new_mat = createMatrix(pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix -1; i++){
        new_mat->pointer_of_data[i] = pMat->pointer_of_data[i] * scalar;
    }
    return new_mat;
}

Matrix *InvertConstMatrix(const Matrix *pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    bool valid1 = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_NOT_A_SQUARE_MATRIX) &&
                  check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_MATRIX_IRREVERSIBLE);
    if(!valid1){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    invert(pMat->pointer_of_data, result->pointer_of_data,pMat->rows_numbers_in_Matrix);
    return result;
}

Matrix *transposeConstMatrix(const Matrix* pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->cols_numbers_in_Matrix, pMat->rows_numbers_in_Matrix);
    transpose(pMat->pointer_of_data, result->pointer_of_data,pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    return result;
}

Matrix *adj_const_matrix(const Matrix* pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    bool valid1 = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_NOT_A_SQUARE_MATRIX) &&
                  check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_MATRIX_IRREVERSIBLE);
    if(!valid1){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->cols_numbers_in_Matrix, pMat->rows_numbers_in_Matrix);
    adjoint_Matrix(pMat->pointer_of_data, result->pointer_of_data,pMat->rows_numbers_in_Matrix);
    return result;
}

Matrix * swapRows_const_matrix(const Matrix *pMat, size_t rows1, size_t rows2){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->cols_numbers_in_Matrix, pMat->rows_numbers_in_Matrix);
    copyMatrix(pMat, result);
    swap_rows(result->pointer_of_data,result->cols_numbers_in_Matrix, rows1-1, rows2-1);
    return result;
}
Matrix * addRows_const_matrix(const Matrix *pMat, size_t rows1, size_t rows2, float coefficient)
{
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->cols_numbers_in_Matrix, pMat->rows_numbers_in_Matrix);
    copyMatrix(pMat, result);
    add_rows(result->pointer_of_data,result->cols_numbers_in_Matrix, rows1-1, rows2-1, coefficient);
    return result;
}
Matrix * multiplyRows_const_matrix(const Matrix *pMat, size_t rows1, float coefficient){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return NULL;
    }
    Matrix * result = createMatrix(pMat->cols_numbers_in_Matrix, pMat->rows_numbers_in_Matrix);
    copyMatrix(pMat, result);
    multiply_rows(result->pointer_of_data,result->cols_numbers_in_Matrix, rows1-1, coefficient);
    return result;
}

Matrix *solve_equation_and_return(const Matrix *pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid) {
        return NULL;
    }
    Matrix * result = createMatrix(pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    copyMatrix(pMat, result);
    Gause_elimination_to_upper_triangle_mat(result->pointer_of_data, result->rows_numbers_in_Matrix, result->cols_numbers_in_Matrix);
    size_t rank_cof = triangle_mat_rank_for_cof(result->pointer_of_data, result->rows_numbers_in_Matrix, result->cols_numbers_in_Matrix);
    size_t rank_aug = triangle_mat_rank_for_aug(result->pointer_of_data, result->rows_numbers_in_Matrix, result->cols_numbers_in_Matrix);

    if(rank_cof<rank_aug){
        printf("Infinity solutions!\n");
    }
    else if(rank_cof<result->cols_numbers_in_Matrix-1){
        printf("No solutions!\n");
    }
    else if(rank_cof == result->cols_numbers_in_Matrix-1){
        Gause_elimination_to_solve_equations(result->pointer_of_data, result->cols_numbers_in_Matrix, rank_cof);
    }
    return NULL;
}

//--------------------------------------------------------------------------------------------//


void addMatrix(Matrix* pMat, ...){
    bool valid = check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    va_list ap;
    va_start(ap, pMat);
    while(true){

        const Matrix * p = va_arg(ap, Matrix*);
        if(!find_node_by_matrix(p)){
            va_end(ap);
            return;
        }
        const Matrix * mats[2] = {pMat, p};
        bool continue_or_not = check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!",
                                           p, ERROR_CAUSE_BY_NULL_POINTER) &&
                               check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!",
                                           p, ERROR_MODIFICATION_IN_MATRIX) &&
                               check_error("addMatrices(const Matrix* pMat, ...)", "Failure in addition!",
                                           mats, ERROR_NOT_THE_SAME_SIZE) ;
        if(!continue_or_not ){
            va_end(ap);
            return;
        }
        else {
            add_two_mat(pMat, p);
        }
    }
    va_end(ap);
}

void subtractMatrix(Matrix* pMat, ...){
    bool valid = check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }

    va_list ap;
    va_start(ap, pMat);
    while(true){
        const Matrix * p = va_arg(ap, Matrix*);
        if(!find_node_by_matrix(p)){
            va_end(ap);
            return;
        }
        const Matrix * mats[2] = {pMat, p};
        bool continue_or_not = check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                                           p, ERROR_CAUSE_BY_NULL_POINTER) &&
                               check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                                           p, ERROR_MODIFICATION_IN_MATRIX) &&
                               check_error("subtractMatrices(int MaxSize, const Matrix* pMat, ...)", "Failure in subtraction!",
                                           mats, ERROR_NOT_THE_SAME_SIZE) ;
        if(!continue_or_not){
            va_end(ap);
            return;
        }
        else {
            subtract_two_mat(pMat, p);
        }
    }
    va_end(ap);
}

void MultiplyMatrix(Matrix* pMat, ...){
    bool valid = check_error("MultiplyMatrix(Matrix* pMat, ...)", "Failure in multiplication!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("MultiplyMatrix(Matrix* pMat, ...)", "Failure in multiplication!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("MultiplyMatrix(Matrix* pMat, ...)", "Failure in multiplication!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    Matrix * first_mat = createMatrix(pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    copyMatrix(pMat, first_mat);
    Matrix * mid = first_mat;
    va_list ap;
    va_start(ap, pMat);
    while(true){
        const Matrix * p = va_arg(ap, Matrix*);
        if(!find_node_by_matrix(p)){
            float *pFloat = mid->pointer_of_data;
            size_t final_rows = mid->rows_numbers_in_Matrix;
            size_t final_cols = mid->cols_numbers_in_Matrix;
            reset(pMat,final_rows, final_cols, pFloat);
            deleteMatrix(&first_mat);
            va_end(ap);
            return;
        }
        const Matrix * mats[2] = {mid, p};
        bool continue_or_not = check_error("MultiplyMatrix(Matrix* pMat, ...)", "Failure in multiplication!",
                                           p, ERROR_CAUSE_BY_NULL_POINTER) &&
                               check_error("MultiplyMatrix(Matrix* pMat, ...)", "Failure in multiplication!",
                                           p, ERROR_MODIFICATION_IN_MATRIX) &&
                               check_error("MultiplyMatrix(Matrix* pMat, ...)", "Failure in multiplication!",
                                           mats, ERROR_INVALID_COLS_AND_ROWS_IN_MULTIPLICATION) ;
        if(!continue_or_not ){
            float *pFloat = mid->pointer_of_data;
            size_t final_rows = mid->rows_numbers_in_Matrix;
            size_t final_cols = mid->cols_numbers_in_Matrix;
            mid = createMatrix(final_rows, final_rows);
            set_Matrix_with_an_array(mid, pFloat, final_rows*final_cols);
            deleteMatrix(&first_mat);
            va_end(ap);
            return;
        }
        else {
            size_t times = mid->cols_numbers_in_Matrix;
            float ans[mid->rows_numbers_in_Matrix*p->cols_numbers_in_Matrix];
            for(int i = 0; i<mid->rows_numbers_in_Matrix; i++){
                for(int j = 0; j<p->cols_numbers_in_Matrix; j++){
                    ans[i*p->cols_numbers_in_Matrix+j] = 0;
                    for(int k = 0; k<times; k++){
                        ans[i*p->cols_numbers_in_Matrix+j] += mid->pointer_of_data[i*mid->cols_numbers_in_Matrix+k] *
                                                              p->pointer_of_data[k*mid->cols_numbers_in_Matrix+j];
                        //printf("%f \n",ans[i*p->cols_numbers_in_Matrix+j]);
                    }
                }
            }
            free(mid->pointer_of_data);
            mid->pointer_of_data = (float *) calloc(mid->rows_numbers_in_Matrix*p->cols_numbers_in_Matrix, sizeof(float));
            if(!check_error("MultiplyMatrix(Matrix* pMat, ...)","failure in allocating memory",mid->pointer_of_data,ERROR_IN_ALLOCATING_MEMORY)){
                return;
            }
            for(int i = 0; i<mid->rows_numbers_in_Matrix*p->cols_numbers_in_Matrix;i++){
                mid->pointer_of_data[i] = ans[i];
            }
            mid->cols_numbers_in_Matrix = p->cols_numbers_in_Matrix;
        }
    }
    va_end(ap);
}

void add_a_scalar(Matrix *pMat, float scalar) {
    bool valid = check_error("add_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("add_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("add_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix; i++){
        pMat->pointer_of_data[i] += scalar;
    }
}

void subtract_a_scalar(Matrix *pMat, float scalar) {
    bool valid = check_error("subtract_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("subtract_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("subtract_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix; i++){
        pMat->pointer_of_data[i] -= scalar;
    }
}

void multiply_a_scalar(Matrix *pMat, float scalar) {
    bool valid = check_error("multiply_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("multiply_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("multiply_a_scalar(Matrix *pMat, float scalar)", "Failure in destruction!", pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    for(int i = 0; i < pMat->rows_numbers_in_Matrix * pMat->cols_numbers_in_Matrix; i++){
        pMat->pointer_of_data[i] *= scalar;
    }
}

void invertMatrix( Matrix *pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    bool valid1 = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_NOT_A_SQUARE_MATRIX) &&
                  check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_MATRIX_IRREVERSIBLE);
    if(!valid1){
        return;
    }
    float data[pMat->cols_numbers_in_Matrix *pMat->rows_numbers_in_Matrix];
    invert(pMat->pointer_of_data, data,pMat->rows_numbers_in_Matrix);
    for(int i = 0; i<pMat->cols_numbers_in_Matrix*pMat->rows_numbers_in_Matrix; i++){
        pMat->pointer_of_data[i] = data[i];
    }
}

void transposeMatrix(Matrix* pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    float data[pMat->cols_numbers_in_Matrix *pMat->rows_numbers_in_Matrix];
    transpose(pMat->pointer_of_data, data,pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    for(int i = 0; i<pMat->rows_numbers_in_Matrix*pMat->rows_numbers_in_Matrix; i++){
        pMat->pointer_of_data[i] = data[i];
    }
    size_t c = pMat->rows_numbers_in_Matrix;
    pMat->rows_numbers_in_Matrix= pMat->cols_numbers_in_Matrix;
    pMat->cols_numbers_in_Matrix = c;
}


void adj_matrix(Matrix* pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    bool valid1 = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_NOT_A_SQUARE_MATRIX) &&
                  check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                              pMat, ERROR_MATRIX_IRREVERSIBLE);
    if(!valid1){
        return;
    }
    float data[pMat->cols_numbers_in_Matrix *pMat->rows_numbers_in_Matrix];
    adjoint_Matrix(pMat->pointer_of_data, data, pMat->rows_numbers_in_Matrix);
    for(int i = 0; i<pMat->rows_numbers_in_Matrix*pMat->rows_numbers_in_Matrix; i++){
        pMat->pointer_of_data[i] = data[i];
    }
}
void solve_equation(Matrix *pMat){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid) {
        return;
    }

    Gause_elimination_to_upper_triangle_mat(pMat->pointer_of_data, pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    size_t rank_cof = triangle_mat_rank_for_cof(pMat->pointer_of_data, pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);
    size_t rank_aug = triangle_mat_rank_for_aug(pMat->pointer_of_data, pMat->rows_numbers_in_Matrix, pMat->cols_numbers_in_Matrix);

    if(rank_cof<rank_aug){
        printf("Infinity solutions!\n");
    }
    else if(rank_cof<pMat->cols_numbers_in_Matrix-1){
        printf("No solutions!\n");
    }
    else if(rank_cof == pMat->cols_numbers_in_Matrix-1){
        Gause_elimination_to_solve_equations(pMat->pointer_of_data, pMat->cols_numbers_in_Matrix, rank_cof);
    }

}
void swapRows(Matrix *pMat, size_t rows1, size_t rows2){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    swap_rows(pMat->pointer_of_data,pMat->cols_numbers_in_Matrix, rows1-1, rows2-1);
}
void addRows(Matrix *pMat, size_t rows1, size_t rows2, float coefficient){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    add_rows(pMat->pointer_of_data,pMat->cols_numbers_in_Matrix, rows1-1, rows2-1,coefficient);
}
void multiplyRows(Matrix *pMat, size_t rows1, float coefficient){
    bool valid = check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_CAUSE_BY_NULL_POINTER) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat,ERROR_NOT_A_MATRIX_CREATE_BY_createMatrix) &&
                 check_error("invert(const Matrix *pMat)", "Failure in inverting matrix!",
                             pMat, ERROR_MODIFICATION_IN_MATRIX);
    if(!valid){
        return;
    }
    multiply_rows(pMat->pointer_of_data,pMat->cols_numbers_in_Matrix, rows1-1, coefficient);
}

