using System;
using System.Drawing;


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

    public double deter(double[,] matrix)
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
            double cofactor = Math.Pow(-1, i + j) * deter(minor); // cofactor is result of a calcul using sub determinant, using only one box/case of the matrix
            determinant += matrix[i, j] * cofactor;
        }

        return determinant;
    }





    public double[,] tran(double[,] matrix)
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

    public double[,] com(double[,] matrix)
    {
        int size = matrix.GetLength(0);
        double[,] comatrix = new double[size, size];

        // Calculate cofactor matrix
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                double[,] minor = getMinorMatrix(matrix, i, j);
                comatrix[i, j] = Math.Pow(-1, i + j) * deter(minor); // IDK WHY BUT WE DO NOT NEED TO MULTIPLY BY ELEMENT HIMSELF -> because it's comatrix calcul
            }
        }

        // transpose to obtain comatrix
        return tran(comatrix);
    }

    public double[,] multiplyMatrixByScalar(double[,] matrix, double scalar)
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

    public double[,] inverse(double[,] matrix)
    {

        
        double[,] inverseMatrix = new double[matrix.GetLength(0), matrix.GetLength(1)];

        double det = deter(matrix);
        double[,] comatrix = com(matrix);


        inverseMatrix = multiplyMatrixByScalar(comatrix, 1 / det);

        dislayMatrix(inverseMatrix);

        return inverseMatrix;
    }

    public double[,] prodMat(double[,] matrix1, double[,] matrix2)
    {
        if (matrix1.GetLength(1) == matrix2.GetLength(0))
        {
            // Determine new matrix size (m-n / p-q) New size = m*q
            int m = matrix1.GetLength(0); // Line
            int q = matrix2.GetLength(1); // Column

            int n = matrix1.GetLength(1);

            double[,] productMatrix = new double[m, q];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < q; j++)
                {
                    productMatrix[i, j] = 0;
                    Console.WriteLine($"Calcul de C[{i},{j}] :");
                    
                    // Sum final result n time
                    for (int k = 0; k < n; k++)
                    {
                        double multiplication = matrix1[i, k] * matrix2[k, j];
                        productMatrix[i, j] += multiplication;

                        Console.WriteLine($"  A[{i},{k}] ({matrix1[i, k]}) * B[{k},{j}] ({matrix2[k, j]}) = {multiplication}  -> Somme partielle: {productMatrix[i, j]}");
                    }

                    Console.WriteLine($"Résultat final C[{i},{j}] = {productMatrix[i, j]}\n");
                }
            }
            return productMatrix;   


        }
        else
        {
            Console.WriteLine("Erreur : Multiplication de matrice pas possible car n n'est pas égale a p");
            throw new Exception("Message d'erreur");
        }
    }

    public double[] prodVect(double[] vect1, double[] vect2)
    {
        double[] finalVect = new double[3];
        finalVect[0] = vect1[1]* vect2[2] - vect1[2] * vect2[1];
        finalVect[1] = vect1[2] * vect2[0] - vect1[0] * vect2[2];
        finalVect[2] = vect1[0] * vect2[1] - vect1[1] * vect2[0];

        return finalVect;
    }

    public double[] subVect(double[] vect1, double[] vect2)
    {
        double[] subVec = new double[3];
        for (int i = 0; i < 3; i++)
        {
            subVec[i] = vect1[i] - vect2[i];
        }

        return subVec;
    }

    public double[] moment(double[] force, double[] applicationVector, double[] inertieCenter)
    {
        double[] agVec = subVect(inertieCenter, applicationVector);
        //Console.WriteLine(agVec[0] + " /" + agVec[1] + " /" + agVec[2]);

        double[] moment = prodVect(agVec, force);

        return moment;
    }
}

    class Program
{
    static void Main()
    {
        MatrixOperations matrixOperations = new MatrixOperations();

        double[] force = { 2, 3, 4 };
        double[] applicationVector = { 1, 1, 1 };
        double[] inertieCenter = { 3, 3, 3 };

        double[] resultat = matrixOperations.moment(force, applicationVector, inertieCenter);

        //double[] A = { 1, 2, 3 };
        //double[] B = { 4, 5, 6 };

        //double[] resultat = matrixOperations.prodVect(A, B);

        Console.WriteLine(resultat[0] + " /" + resultat[1] + " /" + resultat[2]);

        //double[,] matrix1 = new double[2, 3]
        //{
        //    { 1, 2, 3 },
        //    { 4, 5, 6 },
        //};

        //double[,] matrix2 = new double[3, 2]
        //{
        //    { 7, 8},
        //    { 9, 10},
        //    { 11, 12},
        //};

        //matrixOperations.prodMat(matrix1, matrix2); 


        //double[,] matrix = new double[3, 3]
        //{
        //    { 1, 2, 1 },
        //    { 3, 5, 0 },
        //    { 4, 2, 6 }
        //};
        //double[,] inverseMatrix = matrixOperations.inverse(matrix);




    }
}

//-1,5 0,5 0,25
//0,9 - 0,1 - 0,15000000000000002
//0,7000000000000001 - 0,30000000000000004 0,05