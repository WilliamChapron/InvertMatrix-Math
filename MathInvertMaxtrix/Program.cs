﻿using System;
using System.Drawing;
using System.Numerics;


class MatrixOperations
{


    // TP 1

    
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
















    // TP 2
    public double[] prodVect(double[] vect1, double[] vect2)
    {
        double[] finalVect = new double[3];
        finalVect[0] = vect1[1]* vect2[2] - vect1[2] * vect2[1];
        finalVect[1] = vect1[2] * vect2[0] - vect1[0] * vect2[2];
        finalVect[2] = vect1[0] * vect2[1] - vect1[1] * vect2[0];

        return finalVect;
    }

    public double[] addVect(double[] vect1, double[] vect2)
    {

        double[] addVec = new double[3];
        for (int i = 0; i < 3; i++)
        {
            addVec[i] = vect1[i] + vect2[i];
        }

        return addVec;
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

    public double[] momentF(double[] force, double[] applicationPoint, double[] inertieCenter)
    {
        // Calcul du vecteur agVec (distance entre le centre d'inertie et le point d'application de la force)
        double[] agVec = subVect(inertieCenter, applicationPoint);

        // Affichage du vecteur agVec pour débogage
        Console.WriteLine($"Vecteur agVec (inertieCenter - applicationPoint): {string.Join(", ", agVec)}");

        // Calcul du moment en utilisant le produit vectoriel
        double[] moment = prodVect(agVec, force);

        // Affichage du moment calculé pour débogage
        Console.WriteLine($"Moment (Produit vectoriel agVec x force): {string.Join(", ", moment)}");

        return moment;
    }

    public double solve1(double f, double fp, double h)
    {
        return f + fp * h;
    }


    public void translation(double m, double h, double[] forceSum, ref double[] pos, ref double[] speed)
    {
        double[] acceleration = { forceSum[0] / m, forceSum[1] / m, forceSum[2] / m };

        double[] newSpeed = {
            solve1(speed[0], acceleration[0], h),
            solve1(speed[1], acceleration[1], h),
            solve1(speed[2], acceleration[2], h),
        };

        double[] newPosition = {
            solve1(pos[0], newSpeed[0], h),
            solve1(pos[1], newSpeed[1], h),
            solve1(pos[2], newSpeed[2], h),
        };

        // set ref
        for (int i = 0; i < 3; i++)
        {
            speed[i] = newSpeed[i];
        }
        for (int i = 0; i < 3; i++)
        {
            pos[i] = newPosition[i];
        }
    }
    public double[] multiplyMatrixVector(double[,] matrix, double[] vector)
    {

        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1); 

        double[] result = new double[rows]; 

        for (int i = 0; i < rows; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < cols; j++)
            {
                sum += matrix[i, j] * vector[j]; 
            }
            result[i] = sum; 
        }

        return result;
    }

    public void rotation(double h, double[,] forces, double[,] applicationPoints, double[] inertieCenter, double[,] inertieMatrix, ref double[] rotAngle, ref double[] rotSpeed)
    {
        Console.WriteLine("Début de la fonction de rotation");

        // Sum forces moment
        double[] totalMoment = new double[3];
        for (int i = 0; i < forces.GetLength(0); i++)
        {
            double[] force = { forces[i, 0], forces[i, 1], forces[i, 2] };
            double[] applicationPoint = { applicationPoints[i, 0], applicationPoints[i, 1], applicationPoints[i, 2] };

            Console.WriteLine($"Force {i}: {string.Join(", ", force)}");
            Console.WriteLine($"Point d'application {i}: {string.Join(", ", applicationPoint)}");

            double[] moment = momentF(force, applicationPoint, inertieCenter);
            Console.WriteLine($"Moment calculé pour la force {i}: {string.Join(", ", moment)}");

            totalMoment = addVect(totalMoment, moment);
            Console.WriteLine($"Total Moment après ajout: {string.Join(", ", totalMoment)}");
        }

        // sum of moment = I * Ω/ω -> Ω (angular acceleration = I-1 * sum of moment)
        double[,] matrixInverse = inverse(inertieMatrix);
        Console.WriteLine("Matrice d'inertie inverse: ");
        PrintMatrix(matrixInverse);

        double[] Ω = multiplyMatrixVector(matrixInverse, totalMoment);
        Console.WriteLine($"Accélération angulaire (Ω): {string.Join(", ", Ω)}");

        double[] newRotSpeed = {
        solve1(rotSpeed[0], Ω[0], h),
        solve1(rotSpeed[1], Ω[1], h),
        solve1(rotSpeed[2], Ω[2], h),
    };
        Console.WriteLine($"Vitesse angulaire calculée (nouvelle vitesse): {string.Join(", ", newRotSpeed)}");

        double[] newRotAngle = {
        solve1(rotAngle[0], newRotSpeed[0], h),
        solve1(rotAngle[1], newRotSpeed[1], h),
        solve1(rotAngle[2], newRotSpeed[2], h),
    };
        Console.WriteLine($"Angle de rotation calculé (nouvel angle): {string.Join(", ", newRotAngle)}");

        rotSpeed[0] = newRotSpeed[0];
        rotSpeed[1] = newRotSpeed[1];
        rotSpeed[2] = newRotSpeed[2];

        rotAngle[0] = newRotAngle[0];
        rotAngle[1] = newRotAngle[1];
        rotAngle[2] = newRotAngle[2];

        Console.WriteLine("Fin de la fonction de rotation");
    }

    private void PrintMatrix(double[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            Console.WriteLine($"[{matrix[i, 0]}, {matrix[i, 1]}, {matrix[i, 2]}]");
        }
    }


}

    class Program
{
    static void Main()
    {
        MatrixOperations matrixOperations = new MatrixOperations();

        double[,] forces = {
            {0, 100, 0},   
            {0, 0, 0},     
            {0, 0, 0}      
        };
        double[,] applicationPoints = {
            {1, 0, 0},     
            {0, 0, 0},     
            {0, 0, 0}     
        };
        double[] inertieCenter = { 0, 0, 0 };  
        double[,] inertieMatrix = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}     
        };
        double[] rotSpeed = { 0, 0, 0 }; 
        double[] rotAngle = { 0, 0, 0 }; 
        double h = 0.1;  

        matrixOperations.rotation(h, forces, applicationPoints, inertieCenter, inertieMatrix, ref rotAngle, ref rotSpeed);

        //double m = 2.0;
        //double h = 0.1;
        //double[] force = { 6.0, 0.0, 0.0 };
        //double[] inertieCenter = { 1.0, 1.0, 0.0 };
        //double[] position = { inertieCenter[0], inertieCenter[1], inertieCenter[2] };
        //double[] speed = { 0.0, 0.0, 0.0 };


        //matrixOperations.translation(m, h, force, ref position, ref speed);

        //Console.WriteLine("Nouvelle position G après 1 translations :");
        //Console.WriteLine($"  X: {position[0]:F2}");
        //Console.WriteLine($"  Y: {position[1]:F2}");
        //Console.WriteLine($"  Z: {position[2]:F2}");

        //Console.WriteLine("\nNouvelle vitesse v après 1 translation :");
        //Console.WriteLine($"  Vx: {speed[0]:F2}");
        //Console.WriteLine($"  Vy: {speed[1]:F2}");
        //Console.WriteLine($"  Vz: {speed[2]:F2}");



        //matrixOperations.translation(m, h, force, ref position, ref speed);


        //Console.WriteLine("Nouvelle position G après 2 translations :");
        //Console.WriteLine($"  X: {position[0]:F2}");
        //Console.WriteLine($"  Y: {position[1]:F2}");
        //Console.WriteLine($"  Z: {position[2]:F2}");

        //Console.WriteLine("\nNouvelle vitesse v après 2 translations :");
        //Console.WriteLine($"  Vx: {speed[0]:F2}");
        //Console.WriteLine($"  Vy: {speed[1]:F2}");
        //Console.WriteLine($"  Vz: {speed[2]:F2}");



        //double[] force = { 2, 3, 4 };
        //double[] applicationVector = { 1, 1, 1 };
        //double[] inertieCenter = { 3, 3, 3 };

        //double[] resultat = matrixOperations.moment(force, applicationVector, inertieCenter);

        ////double[] A = { 1, 2, 3 };
        ////double[] B = { 4, 5, 6 };

        ////double[] resultat = matrixOperations.prodVect(A, B);

        //Console.WriteLine(resultat[0] + " /" + resultat[1] + " /" + resultat[2]);

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