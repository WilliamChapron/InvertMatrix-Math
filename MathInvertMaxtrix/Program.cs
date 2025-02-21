using System;
using System.Drawing;
using System.Numerics;
using System.Reflection;


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
        double[] agVec = subVect(inertieCenter, applicationPoint);

        Console.WriteLine($"Vecteur agVec (inertieCenter - applicationPoint): {string.Join(", ", agVec)}");

        double[] moment = prodVect(agVec, force);

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
        // Sum forces moment
        double[] totalMoment = new double[3];
        for (int i = 0; i < forces.GetLength(0); i++)
        {
            double[] force = { forces[i, 0], forces[i, 1], forces[i, 2] };
            double[] applicationPoint = { applicationPoints[i, 0], applicationPoints[i, 1], applicationPoints[i, 2] };

            double[] moment = momentF(force, applicationPoint, inertieCenter);

            totalMoment = addVect(totalMoment, moment);
        }

        // sum of moment = I * Ω/ω -> Ω (angular acceleration = I-1 * sum of moment)
        double[,] matrixInverse = inverse(inertieMatrix);

        double[] Ω = multiplyMatrixVector(matrixInverse, totalMoment);

        double[] newRotSpeed = {
            solve1(rotSpeed[0], Ω[0], h),
            solve1(rotSpeed[1], Ω[1], h),
            solve1(rotSpeed[2], Ω[2], h),
        };

        double[] newRotAngle = {
            solve1(rotAngle[0], newRotSpeed[0], h),
            solve1(rotAngle[1], newRotSpeed[1], h),
            solve1(rotAngle[2], newRotSpeed[2], h),
        };

        rotSpeed[0] = newRotSpeed[0];
        rotSpeed[1] = newRotSpeed[1];
        rotSpeed[2] = newRotSpeed[2];

        rotAngle[0] = newRotAngle[0];
        rotAngle[1] = newRotAngle[1];
        rotAngle[2] = newRotAngle[2];
    }

    private void PrintMatrix(double[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            Console.WriteLine($"[{matrix[i, 0]}, {matrix[i, 1]}, {matrix[i, 2]}]");
        }
    }

    public double[,] centre_inert(double[,] points)
    {
        double totalMass = 0;
        double sumX = 0, sumY = 0, sumZ = 0;

        for (int i = 0; i < points.GetLength(0); i++)
        {
            double mass = points[i, 0];
            double x = points[i, 1];
            double y = points[i, 2];
            double z = points[i, 3];

            totalMass += mass;
            sumX += x * mass;
            sumY += y * mass;
            sumZ += z * mass;
        }

        double centerX = sumX / totalMass;
        double centerY = sumY / totalMass;
        double centerZ = sumZ / totalMass;

        return new double[,] { { centerX, centerY, centerZ } };
    }

    public double[,] matrice_inert(double[,] points)
    {
        double[,] I = new double[3, 3];


        double A = 0, B = 0, C = 0;

        double D = 0, E = 0, F = 0;

        for (int i = 0; i < points.GetLength(0); i++)
        {
            double x = points[i, 0];  // X
            double y = points[i, 1];  // Y
            double z = points[i, 2];  // Z
            double m = points[i, 3];  // m


            A += m * (Math.Pow(y, 2) + Math.Pow(z, 2));
            B += m * (Math.Pow(x, 2) + Math.Pow(z, 2));  
            C += m * (Math.Pow(x, 2) + Math.Pow(y, 2)); 

            D += m * y * z;
            E += m * x * z;
            F += m * x * y;
        }

        // Diagonal
        I[0, 0] = A;
        I[1, 1] = B;
        I[2, 2] = C;



        // Non diagonal
        I[0, 1] = -F;
        I[1, 0] = -F;

        I[0, 2] = -E;
        I[2, 0] = -E;

        I[1, 2] = -D;
        I[2, 1] = -D;



        return I;
    }

    public double[,] displace_matrix(double[,] I_O, double totalMass, double[] oPoint, double[] aPoint)
    {
        double a = aPoint[0] - oPoint[0];
        double b = aPoint[1] - oPoint[1];
        double c = aPoint[2] - oPoint[2];

        // final matrix (moved to A)
        double[,] I_A = new double[3, 3];

        // calculate displacement matrix adding matrix O, the both are made in same time 

        // line 1
        I_A[0, 0] = I_O[0, 0] + totalMass * (Math.Pow(b, 2) + Math.Pow(c, 2));
        I_A[0, 1] = I_O[0, 1] - totalMass * a * b;
        I_A[0, 2] = I_O[0, 2] - totalMass * a * c;

        // line 2
        I_A[1, 0] = I_O[1, 0] - totalMass * a * b;
        I_A[1, 1] = I_O[1, 1] + totalMass * (Math.Pow(a, 2) + Math.Pow(c, 2));
        I_A[1, 2] = I_O[1, 2] - totalMass * b * c;

        // line 3
        I_A[2, 0] = I_O[2, 0] - totalMass * a * c;
        I_A[2, 1] = I_O[2, 1] - totalMass * b * c;
        I_A[2, 2] = I_O[2, 2] + totalMass * (Math.Pow(a, 2) + Math.Pow(b, 2));

        return I_A;
    }

    // TP 4

    public double[,] pave_plein(int n, double a, double b, double c, double x0, double y0, double z0)
    {
        // approximation of number of points by direction (x,y,z)
        int iterations = (int)Math.Pow(n, 1.0 / 3.0); 

        // x,y,z lists
        double[,] result = new double[n, 3];



        // calculate offset
        double dx = a / (iterations - 1);
        double dy = b / (iterations - 1);
        double dz = c / (iterations - 1);

        int index = 0;

        // coherent loop :)
        for (int i = 0; i < iterations; i++)
        {
            for (int j = 0; j < iterations; j++)
            {
                for (int k = 0; k < iterations; k++)
                {
                    // base pos + (iteration value between 0 and "iterations" multiply by offset dx/dy or dz)
                    result[index, 0] = x0 + i * dx;  
                    result[index, 1] = y0 + j * dy;  
                    result[index, 2] = z0 + k * dz;
                    // add all variations possible of points for each axes coordinate

                    index++;
                }
            }
        }

        return result;
    }

    public static long factoriel(int n)
    {
        long result = 1;
        for (int i = 1; i <= n; i++)
        {
            result *= i;
        }
        return result;
    }

    public double ld_cosinus(double x)
    {
        // begin by 1
        double result = 1;
        // also set to 1 :)
        double term = 1;

        int n = 2;
        int k = -1;

        // Approximation until x^6 because begin x^2
        while (n <= 6)
        {
            term *= k * (Math.Pow(x, n) / factoriel(n)); 
            result += term;

            // increment factoriel
            n += 2;
            // change sign
            k *= -1;
        }

        return result;
    }

    public double ld_sinus(double x)
    {
        // begin by x
        double result = x;
        // also set to x :)
        double term = x;

        int n = 3;
        int k = -1;

        // Approximation until x^5 because begin x^3
        while (n <= 5)
        {
            term *= k * (Math.Pow(x, n) / factoriel(n));
            result += term;

            // increment factoriel
            n += 2;
            // change sign
            k *= -1;
        }

        return result;
    }
}

    class Program
{
    static void Main()
    {
        MatrixOperations matrixOperations = new MatrixOperations();

        //int n = 27;
        //double a = 3.0; 
        //double b = 3.0; 
        //double c = 3.0; 
        //double x0 = 0.0; 
        //double y0 = 0.0; 
        //double z0 = 0.0; 
        //double[,] points = matrixOperations.pave_plein(n, a, b, c, x0, y0, z0);
        //for (int i = 0; i < points.GetLength(0); i++)
        //{
        //    Console.WriteLine($"Point {i + 1}: X = {points[i, 0]}, Y = {points[i, 1]}, Z = {points[i, 2]}");
        //}

        //double[,] I_O = new double[,] 
        //{
        //    { 6.0, -2.0, -1.0 },
        //    { -2.0, 6.0, -1.0 },
        //    { -1.0, -1.0, 9.0 }
        //};
        //double m = 10.0; 
        //double[] O = { 0.0, 0.0, 0.0 }; 
        //double[] A = { 1.0, 1.0, 1.0 };  

        //double[,] I_A = matrixOperations.displace_matrix(I_O, m, O, A);

        //Console.WriteLine("Matrice d'inertie déplacée :");
        //for (int i = 0; i < 3; i++)
        //{
        //    for (int j = 0; j < 3; j++)
        //    {
        //        Console.Write(I_A[i, j] + " ");
        //    }
        //    Console.WriteLine();
        //}




        //double[,] points = new double[,]
        //{
        //            { 0, 0, 0, 2 },  
        //            { 1, 0, 0, 2 },  
        //            { 1, 1, 0, 2 },  
        //            { 0, 1, 0, 2 },  
        //            { 0.5, 0.5, 1, 2 }  
        //};

        //double[,] matriceInertie = matrixOperations.matrice_inert(points);


        //Console.WriteLine("Matrice d'inertie :");
        //Console.WriteLine($"[{matriceInertie[0, 0]} , {matriceInertie[0, 1]} , {matriceInertie[0, 2]}]");
        //Console.WriteLine($"[{matriceInertie[1, 0]} , {matriceInertie[1, 1]} , {matriceInertie[1, 2]}]");
        //Console.WriteLine($"[{matriceInertie[2, 0]} , {matriceInertie[2, 1]} , {matriceInertie[2, 2]}]");



        //double[,] points = new double[,]
        //{
        //{ 1, 1, 2, 3 }, 
        //{ 1, 2, 5, 6 },  
        //{ 2, 3, 8, 9 }  
        //};

        //double[,] centre = matrixOperations.centre_inert(points);
        //Console.WriteLine("Centre de masse (inertie) : X = " + centre[0, 0] + ", Y = " + centre[0, 1] + ", Z = " + centre[0, 2]);


        //double[,] forces = {
        //    {0, 100, 0},   
        //    {0, 0, 0},     
        //    {0, 0, 0}      
        //};
        //double[,] applicationPoints = {
        //    {1, 0, 0},     
        //    {0, 0, 0},     
        //    {0, 0, 0}     
        //};
        //double[] inertieCenter = { 0, 0, 0 };  
        //double[,] inertieMatrix = {
        //    {1, 0, 0},
        //    {0, 1, 0},
        //    {0, 0, 1}     
        //};
        //double[] rotSpeed = { 0, 0, 0 }; 
        //double[] rotAngle = { 0, 0, 0 }; 
        //double h = 0.1;  

        //matrixOperations.rotation(h, forces, applicationPoints, inertieCenter, inertieMatrix, ref rotAngle, ref rotSpeed);

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