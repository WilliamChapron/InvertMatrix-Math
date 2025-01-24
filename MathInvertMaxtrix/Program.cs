using System;


class MatrixOperations
{
    // Take a base matrix, and remove a line and a column from specified, used in determinant calculation
    public double[,] getMinorMatrix(double[,] matrix, int removedI, int removedJ)
    {
        int size = matrix.GetLength(0);

        double[,] minor = new double[size - 1, size - 1];

        // minor Row/Column index for avoid take i/j that are index of 3x3 matrix , where as we minimize this matrix to 2x2 for example 
        int r = 0, c;
        for (int i = 0; i < size; i++)
        {
            if (i == removedI) continue;  // don't copy value on removed row

            c = 0;
            for (int j = 0; j < size; j++)
            {
                if (j == removedJ) continue;  // don't copy value on removed column

                minor[r, c] = matrix[i, j];
                c++;
            }
            r++;
        }
        return minor;
    }

    public double calculateDet(double[,] matrix)
    {
        int size = matrix.GetLength(0);

        // matrix 1x1
        if (size == 1)
            return matrix[0, 0];

        // matrix 2x2
        if (size == 2)
            return (matrix[0, 0] * matrix[1, 1]) - (matrix[0, 1] * matrix[1, 0]);

        double determinant = 0;

        // 
        int j = 0;
        for (int i = 0; i < size; i++)
        {
            double[,] minor = getMinorMatrix(matrix, i, 0);
            double cofactor = Math.Pow(-1, i + j) * calculateDet(minor); // cofactor is result of a calcul using sub determinant, using only one box/case of the matrix
            determinant += matrix[i, j] * cofactor;
        }

        return determinant;
    }





    public double[,] transposeMatrix(double[,] matrix)
    {
        int rows = matrix.GetLength(0); 
        int cols = matrix.GetLength(1); 

        double[,] transposed = new double[cols, rows]; 

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                transposed[j, i] = matrix[i, j];
            }
        }

        return transposed;
    }

    public double[,] calculateComatrix(double[,] matrix)
    {
        int size = matrix.GetLength(0);
        double[,] comatrix = new double[size, size];

        // Calculate cofactor matrix
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                double[,] minor = getMinorMatrix(matrix, i, j);
                comatrix[i, j] = Math.Pow(-1, i + j) * calculateDet(minor); // IDK WHY BUT WE DO NOT NEED TO MULTIPLY BY ELEMENT HIMSELF
            }
        }

        // transpose to obtain comatrix
        return transposeMatrix(comatrix);
    }

    public double[,] multiplyMatrixByDouble(double[,] matrix, double scalar)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        double[,] result = new double[rows, cols];

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                result[i, j] = matrix[i, j] * scalar;
            }
        }

        return result;
    }

    public void dislayMatrix(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                Console.Write(matrix[i, j] + " ");
            }
            Console.WriteLine();
        }
    }

}

    class Program
{
    static void Main()
    {
        MatrixOperations matrixOperations = new MatrixOperations();
        double[,] matrix = new double[3, 3]
        {
            { 1, 2, 1 },
            { 3, 5, 0 },
            { 4, 2, 6 }
        };


        double det = matrixOperations.calculateDet(matrix);
        double[,] comatrix = matrixOperations.calculateComatrix(matrix);


        double[,] invertMatrix = matrixOperations.multiplyMatrixByDouble(comatrix, 1 / det);

        matrixOperations.dislayMatrix(invertMatrix); 
    }
}

