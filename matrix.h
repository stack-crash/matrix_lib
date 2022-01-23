//matrix.h is a general purpose matrix library for
//data manipulation and data parsing in C. A C++
//implementation will follow to permit the use of
//member functions and paralleleized processes

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <time.h>

#define TRUE 1
#define FALSE 0

#define min(x, y) ((x > y) ? y : x)
#define max(x, y) ((x < y) ? y : x)

#define MAX_MAT_DISPLAY_DIM 8

#define BUFFLEN 0xFF

//Add more specialty formats later
//Note: this is a performance optimization
//To efficiently calculate a determinant through LU decomp,
//the important half (U) is simply the row echelon form of
//a matrix with row swapping disallowed. The determinant is
//just the product of the entries in its diagonal. ech_matrix(),
//has *one* line of code that immediately orders the rows by pivots
//to simplify and expedite the following operations. Disabling
//this line through a "format" property of the matrix
//structure seems more appropriate than passing a boolean 
//(sort vs. no sort) every time the function is called, and
//copying the whole function minus this one line seems redundant.

//One obvious problem here is the fact that a matrix cannot be in
//row echelon form AND have row swapping disallowed simultaneously
//by the current format definition, but specific cases in which this
//might be a serious issue aren't obvious, so it's being left for now.
enum {NONE, IN_ROW_ECH, NO_ROW_SWAP};

//The matrix structure used everywhere below
typedef struct {
	unsigned short int rows;
	unsigned short int cols;
	unsigned char format;
	double ** entries;
} MAT_t;


//define a minimum exponent for which 2^-n is just a zero, fixes numbers blowing up when divided by ~0
//To find exponent, double subtracts 0x3FF from the exponent bits (11 bits right after sign bit)
//Thus, for a min coeff of 2^-32, or exp: -0x20, use 0x3FF - 0x020 = 0x3DF
//min coeff of 2^-3 has exp 0x3 and 0x3FF - 0x003 = 3FC
//Benefit? Called frequently in row operations so must be a rapid bitwise operation
#define IEEE754_MIN_EXP 0x3F0

typedef struct {
	long int mant : 52;
	int expn : 11;
	int sign : 1;
} ieee754_double;

typedef union {
	double fmt_double;
	ieee754_double bits;
} bitmap;

double bound_double(double input){
	//11-bit exponent in a double
	bitmap map;

	map.fmt_double = input;

	if ((map.bits.expn & 0x7FF) <= IEEE754_MIN_EXP) return 0;

	return input;
}

//Print dimension-formatted matrix type (can also be vector or singular element)
void print_matrix(MAT_t * mat){
	
	if (mat == NULL){
		puts("attempted to print from uninit. matrix. Halting...");
		exit(-1);
	}

	unsigned char row_bound = min(MAX_MAT_DISPLAY_DIM, mat->rows);
	unsigned char col_bound = min(MAX_MAT_DISPLAY_DIM, mat->cols);

	for (int i = 0, j; i < row_bound; i++){
		
		for (j = 0; j < col_bound; j++){
			if (mat->entries[i][j]){
				printf("%+.2f\t", mat->entries[i][j]);
			} else {
				printf("+0\t");
			}
		}

		puts((col_bound < mat->cols) ? "[ ... ]" : "");

	}

	if (row_bound < mat->rows || col_bound < mat->cols) {

		for (int i = 0; i < col_bound; i++) printf("[ ... ]\t");

		printf("\n\ntruncated with actual size of [%d x %d]\n\n", mat->rows, mat->cols);
	}
}

void print_vector(MAT_t * vector){
	if (vector->rows == 1) /*row vector*/ {

		for (int col = 0; col < vector->cols; col++) printf("%+0.6g,  ", vector->entries[0][col]);
	
	} else if (vector->cols == 1) /*col vector*/ {
		
		for (int row = 0; row < vector->rows; row++) printf("%+0.6g\n", vector->entries[row][0]);

	} else {
		puts("Error. Invalid dimensions for print_vector() call.");
		exit(-1);
	}

}

//Memory allocation functions
MAT_t * alloc_matrix(int rows, int cols){
	
	MAT_t * mat = (MAT_t *)calloc(sizeof(MAT_t), 1);

	double ** matrix = (double **)calloc(rows, sizeof(double *));

	for (int i = 0; i < rows; i++){
		matrix[i] = (double *)calloc(cols, sizeof(double));
	}

	mat->rows = rows;
	mat->cols = cols;
	mat->format = NONE;
	mat->entries = matrix;
	
	return mat;
}

MAT_t * get_col(MAT_t * mat, int col){
	MAT_t * vector = alloc_matrix(mat->rows, 1);

	vector->format = NONE;

	for (int row = 0; row < mat->rows; row++){
		vector->entries[row][0] = mat->entries[row][col];
	}

	return vector;
}

MAT_t * get_row(MAT_t * mat, int row){
	MAT_t * vector = alloc_matrix(1, mat->cols);

	vector->format = NONE;

	for (int col = 0; col < mat->cols; col++){
		vector->entries[0][col] = mat->entries[row][col];
	}

	return vector;
}

void delete_matrix(MAT_t * mat){

	for (int row = 0; row < mat->rows; row++){	
		free(mat->entries[row]);
		mat->entries[row] = NULL;
	}

	mat->cols = 0;
	mat->rows = 0;

	free(mat->entries);
	mat->entries = NULL;

	free(mat);
	mat = NULL;
}

//Save / load / read / make operations
//(This group is effectively any function that alloc's space)
MAT_t * load_matrix(char * filename){

	FILE * fp = fopen(filename, "r");

	if (fp == NULL){
		printf("Error: Tried to read from '%s' but found no matches.\n", filename);
		exit(-1);
	}

	char * buff = (char *)calloc(sizeof(char), BUFFLEN);

	fgets(buff, BUFFLEN, fp);

	int rows = 1, cols = 1;

	//Gets row / col count of csv file (trust me, reading the file twice is,
	//in fact, the most efficient way of doing a static matrix allocation
	//(a growing linked list would be very relatively inefficient for operations)
	while (fgets(buff, BUFFLEN, fp)){	

		for (int i = 0; (buff[i]); i++){

			if ((rows == 1) && (buff[i] == ',')){
				cols += 1;
			}
		}
		rows++;
	}

	rewind(fp);

	MAT_t * mat = alloc_matrix(rows, cols);

	printf("loading [%dx%d] matrix to %p\n\n", rows, cols, mat);

	int i, row = 0, col = 0;

	while (fgets(buff, BUFFLEN, fp) && (row < rows)){
		
		i = 0;

		while ((buff[i] != '\n') && (i < BUFFLEN) && (col < cols)){
			
			mat->entries[row][col] = atof((char *)(buff + i));

			while ((buff[i] != ',') && (buff[i] != '\n')) i++;

			i += 1;
			
			col++;
		}

		col = 0;
		row++;
	}	

	free(buff);
	buff = NULL;

	fclose(fp);
	fp = NULL;
	
	mat->format = NONE;

	return mat;	
}

//Does same as above, but instead of reading a file, gets contents of a csv through stdin
MAT_t * gets_matrix(char * meme){

	puts("Error. Fix me!");

	exit(-1);

	return NULL;
}

MAT_t * identity(int dim){

	MAT_t * mat = alloc_matrix(dim, dim);

	for (int i = 0; i < dim; i++){
		mat->entries[i][i] = 1;
	}

	return mat;
}

MAT_t * mat_rand_int(int rows, int cols, int min, int max, int seed){

	if (min >= max){
		puts("Error. Invalid min and max assignments in mat_rand_int() call.");
		exit(-1);
	}

	srand(time(0) + seed);

	MAT_t * new = alloc_matrix(rows, cols);

	for (int row = 0; row < rows; row++){
		for (int col = 0; col < cols; col++){
			new->entries[row][col] = (rand() % (max - min)) + min;
		}
	}

	return new;
}

//Matrix operations and helper functions
int get_pivot(MAT_t * mat, int row){
	int i = 0;
	while ((i < mat->cols) && (mat->entries[row][i]==0)) i++;

	return i;
}

void swap_rows(MAT_t * mat, int r1, int r2){
	
	double * buff = mat->entries[r1];

	mat->entries[r1] = mat->entries[r2];

	mat->entries[r2] = buff;

	buff = NULL;
}

void order_pivots(MAT_t * mat){

	int piv_pos_a = 0, 
		piv_pos_b = 0;

	//Index to next row. Scan every row below it. 
	//If a row below it has a more significant pivot, swap the rows.
	for (int row = 0; row < mat->rows - 1; row++){

		for (int i = row + 1; i < mat->rows; i++){

			//find pivot position of slower indexing row
			for(piv_pos_a = 0; 	(piv_pos_a < mat->cols) && 
							(mat->entries[row][piv_pos_a] == 0); piv_pos_a++);

			//if pivot_pos(mat->entries[i]) < pivot_pos(mat->entries[row]) swap rows i and row
			for(piv_pos_b = 0; 	(piv_pos_b < mat->cols) && 
								(mat->entries[i][piv_pos_b] == 0); piv_pos_b++);
			
			if (piv_pos_b < piv_pos_a) swap_rows(mat, i, row);
		}
	}

	mat->format = IN_ROW_ECH;
}

void ech_matrix(MAT_t * mat){

	//reference for manual algorithm:
	//http://www.math.odu.edu/~bogacki/cgi-bin/lat.cgi?c=ref

	if (mat->format == IN_ROW_ECH) return;

	//not having this leads to weird echelon formatting errors like once in 15 random sets.
	//see "read.csv" for an example of a set that fails the picto ordering process
	if (mat->format != NO_ROW_SWAP) order_pivots(mat);

	double scaling_factor;

	double normalizing_constant;

	for (int row = 0, subt_row, pivot; row < mat->rows; row++){

		pivot = get_pivot(mat, row);

		normalizing_constant = mat->entries[row][pivot];

		//subtract row([pivot in subt_row]*[row]) from [subt_row]
		for (subt_row = row + 1; subt_row < mat->rows; subt_row++){

			scaling_factor = mat->entries[subt_row][pivot];

			//traverse rows horizontally to complete row subtraction
			for (int i = 0; i < mat->cols; i++){
				mat->entries[subt_row][i] -= 
				((scaling_factor/normalizing_constant) * mat->entries[row][i]);

				mat->entries[subt_row][i] = bound_double(mat->entries[subt_row][i]);
			}
		}
	}

	mat->format=IN_ROW_ECH;
}

void ref_matrix(MAT_t * mat){

	if (mat->format != IN_ROW_ECH) ech_matrix(mat);

	int pivot;

	double normalizing_constant;

	for (int row = 0; row < mat->rows; row++){

		//double_bound entries in row to find pivots more accurately
		for (int i = 0; i < mat->cols; i++) mat->entries[row][i] = bound_double(mat->entries[row][i]);

		pivot = get_pivot(mat, row);


		normalizing_constant = mat->entries[row][pivot];
		//normalize row by dividing pivot and every row member 
		//after pivot by "normalizing constant" (itself)
		for (int col = pivot; col < mat->cols; col++){
			mat->entries[row][col] /= normalizing_constant;
		}
	}
}

void rref_matrix(MAT_t * mat){
	
	//http://www.math.odu.edu/~bogacki/cgi-bin/lat.cgi

	ref_matrix(mat);

	//find lowest row with pivot, and use it to clear entries in rows above. Repeat.

	int pivot = 0, row;

	double scaling_constant;

	for (row = mat->rows - 1; (row > 0); row--){

		pivot = get_pivot(mat, row);

		for (int subt_row = row - 1; subt_row >= 0; subt_row--){

			scaling_constant = mat->entries[subt_row][pivot];

			for (int subt_col = pivot; subt_col < mat->cols; subt_col++){
				mat->entries[subt_row][subt_col] 
					-= (scaling_constant * mat->entries[row][subt_col]);
			}
		}
	}
}

void mat_scale(MAT_t * mat, double scale){
	for (int row = 0; row < mat->rows; row++){
		for (int col = 0; col < mat->cols; col++){
			mat->entries[row][col] *= scale;
		}
	}
	return;
}

MAT_t * mat_mul(MAT_t * mat1, MAT_t * mat2){

	if (mat1->cols != mat2->rows){
		printf("Error. Attempted to multiply [%02dx%02d] and [%02dx%02d] matrices.\n", mat1->rows, mat1->cols, mat2->rows, mat2->cols);
		exit(-1);
	}

	MAT_t * mat = alloc_matrix(mat1->rows, mat2->cols);

	int entry_summation_bound = mat1->cols;

	double buff = 0;

	for (int row = 0; row < mat->rows; row++){
		for (int col = 0; col < mat->cols; buff = 0, col++){
			buff = 0;
			for (int i = 0; i < entry_summation_bound; i++){
				buff += mat1->entries[row][i] * mat2->entries[i][col];
			}


			mat->entries[row][col] = buff;
		}
	}

	return mat;
}

void mat_add(MAT_t * mat1, MAT_t * mat2){
	if ((mat1->cols != mat2->cols) || (mat1->rows != mat2->rows)){
		puts("[!] dimension mismatch. Halting...");
		exit(-1);
	}

	for (int row = 0; row < mat1->rows; row++){
		for (int col = 0; col < mat1->cols; col++){
			mat1->entries[row][col]+= mat2->entries[row][col];
		}
	}
	
	//if both are in row echelon form, the sum will be as well
	if (mat1->format == IN_ROW_ECH &&
		mat2->format == IN_ROW_ECH) {
			mat1->format = (mat1->format && mat2->format);
		} else {
			mat1->format = NONE;
		}

	return;
}

void mat_sub(MAT_t * mat1, MAT_t * mat2){
	if ((mat1->cols != mat2->cols) || (mat1->rows != mat2->rows)){
		puts("[!] dimension mismatch. Halting...");
		exit(-1);
	}

	for (int row = 0; row < mat1->rows; row++){
		for (int col = 0; col < mat1->cols; col++){
			mat1->entries[row][col]-= mat2->entries[row][col];
		}
	}
	
	//if both are in row echelon form, the sum will be as well
	if (mat1->format == IN_ROW_ECH &&
		mat2->format == IN_ROW_ECH) {
			mat1->format = (mat1->format && mat2->format);
		} else {
			mat1->format = NONE;
		}

	return;
}

MAT_t * mat_copy(MAT_t * mat){

	MAT_t * duplicate = alloc_matrix(mat->rows, mat->cols);

	duplicate->format = mat->format;

	for (int row = 0; row < mat->rows; row++){
		for (int col = 0; col < mat->cols; col++){
			duplicate->entries[row][col] = mat->entries[row][col];
		}
	}

	return duplicate;
}

//Matrix transpose function
MAT_t * mat_tps(MAT_t * mat){
	MAT_t * mat_tp = alloc_matrix(mat->cols, mat->rows);

	for (int row = 0; row < mat->rows; row++){
		for (int col = 0; col < mat->cols; col++){
			mat_tp->entries[col][row] = mat->entries[row][col];
		}
	}

	return mat_tp;
}

double mat_det(MAT_t * mat){
	//This is an implementation using LU decomposition
	//Only the matrix U is needed in finding det(mat)

	if (mat->rows != mat->cols){
		puts("Attempted to find determinant of non-square matrix. Halting...");
		exit(-1);
	}

	MAT_t * U = mat_copy(mat);

	U->format = NO_ROW_SWAP;

	ech_matrix(U);

	double determinant = 1;

	for (int i = 0; i < U->rows; i++){
		determinant *= U->entries[i][i];
	}

	delete_matrix(U);

	return determinant;
}

MAT_t * mat_drop_col(MAT_t * mat, int index){
	if (index >= mat->cols){
		puts("Error. Attempted to drop invalid index on call to mat_drop_col().");
		exit(-1);
	} else if (mat->cols == 1){
		puts("Error. Attempted to drop column from singular column matrix (vector) on call to mat_drop_col().");
		exit(-1);
	}

	MAT_t * new = alloc_matrix(mat->rows, mat->cols-1);

	for (int row = 0; row < new->rows; row++){
		for (int col = 0, skip = 0; col < new->cols; col++){
			
			if (!skip && col == index) skip = 1;

			new->entries[row][col] = mat->entries[row][col + skip];
		}
	}

	return new;
}

MAT_t * mat_trn(MAT_t * mat){
	MAT_t * mat_tp = alloc_matrix(mat->cols, mat->rows);

	for (int row = 0; row < mat->rows; row++){
		for (int col = 0; col < mat->cols; col++){
			mat_tp->entries[col][row] = mat->entries[row][col];
		}
	}

	return mat_tp;
}

//adapted from template
//optimize for sparse matrices
MAT_t * mat_inv(MAT_t * mat) {

	if (mat->rows != mat->cols) {
		puts("Error. Nonsquare input for call to invert_matrix().");
		exit(-1);
	}

	if (mat_det(mat) == 0) {
		puts("Error. Non-invertable matrix A found with det(A) = 0 on call to invert_matrix().");
		exit(-1);
	}

	MAT_t * inv = identity(mat->rows);

	double temp;

	int i, j, k;

	/*---------------LoGiC starts here------------------*/	//procedure to make the matrix mat to unit matrix
	for(k=0;k<mat->rows;k++)								//by some row operations,and the same row operations of
	{														//Unit mat. I gives the inverse of matrix mat
		temp=mat->entries[k][k];							//'temp' stores the mat[k][k] value so that [mat[k][k] will not change
		for(j=0;j<mat->rows;j++)							//during the operation mat[i][j]/=mat[k][k] when i=j=k
		{
			mat->entries[k][j]/=temp;						//it performs the following row operations to make mat to unit matrix
			inv->entries[k][j]/=temp;						//R0=R0/mat[0][0],similarly for I also R0=R0/mat[0][0]
		}													//R1=R1-R0*mat[1][0] similarly for I
		for(i=0;i<mat->rows;i++)							//R2=R2-R0*mat[2][0]		,,
		{
			temp=mat->entries[i][k];						//R1=R1/mat[1][1]
			for(j=0;j<mat->rows;j++)						//R0=R0-R1*mat[0][1]
			{												//R2=R2-R1*mat[2][1]
				if(i==k)
					break;									//R2=R2/mat[2][2]
				mat->entries[i][j]-=mat->entries[k][j]*temp;//R0=R0-R2*mat[0][2]
				inv->entries[i][j]-=inv->entries[k][j]*temp;//R1=R1-R2*mat[1][2]
			}
		}
	}
	return inv;
}

MAT_t * lin_reg(MAT_t * mat){
	//Given:
	//  mat = [a|b] where a and b are point vectors
	//  and A = horizontal, right-to-left polynomial expansion of vector a
	//	x ~= ((AT*A)^-1) * (AT*b)

	if (mat->cols != 2){
		puts("Expected two column data matrix in call to mat_least_squares");
		exit(-1);
	}

	MAT_t * b = get_col(mat, mat->cols-1);

	MAT_t * A = mat_drop_col(mat, mat->cols-1);

	MAT_t * AT = mat_trn(A);

	MAT_t * ATA = mat_mul(AT, A);

	MAT_t * pseudo_inv = mat_inv(ATA);

	MAT_t * ATb = mat_mul(AT, b);

	MAT_t * X = mat_mul(pseudo_inv, ATb);

	delete_matrix(b);
	delete_matrix(A);
	delete_matrix(AT);
	delete_matrix(ATA);
	delete_matrix(pseudo_inv);
	delete_matrix(ATb);

	return X;
}