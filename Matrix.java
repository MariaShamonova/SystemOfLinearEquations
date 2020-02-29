package com.metanit;
import java.text.NumberFormat;
import java.util.Random;

public class Matrix {

    public double[][] data;
    public double[] matrixB;

    public int N;

    public Matrix(int size) {
        N = size;

        data = new double[N][];
        matrixB = new double[N];
        for (int i = 0; i < size; i++) {
            data[i] = new double[N];
            for (int j = 0; j < size; j++) {
                data[i][j] = 0.0;
            }
        }
    }

    public Matrix Mult(Matrix matrix) {
        Matrix result = new Matrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    result.data[i][j] += data[i][k] * matrix.data[k][j];
                }
            }
        }
        return result;
    }

    public void Generate(double q) {
        double sum;
        for (int i = 0; i < N; i++) {
            sum = 0;;
            for (int j = 0; j < N; j++) {
                if (i != j) data[i][j] = (int) (Math.random() * (20 + 1)) - 10;
                sum+= Math.abs(data[i][j]);
            }
            data[i][i] = sum * q;
            matrixB[i] = (int) (Math.random() * (20 + 1)) - 10;
        }
    }

    public void Print() {
        NumberFormat nf = NumberFormat.getNumberInstance();
        nf.setMaximumFractionDigits(2);
        nf.setGroupingUsed(false);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.print(nf.format(data[i][j]) + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    public void PrintArray(double[] array) {
        NumberFormat nf = NumberFormat.getNumberInstance();
        nf.setMaximumFractionDigits(2);
        nf.setGroupingUsed(false);

        for (int i = 0; i < N; i++)
        {
            System.out.println(array[i]);
        }
        System.out.println();
    }

    public void SetE() {
        for (int i = 0; i < N; i++)
            data[i][i] = 1;
    }

    public void Copy(Matrix matrix) {
        for (int i = 0; i < N; i++)
            data[i] = matrix.data[i].clone();
    }

    public void SwapRow(int indexRow1, int indexRow2) {
        double tmp;
        for (int j = 0; j < N; j++) {
            tmp = data[indexRow1][j];
            data[indexRow1][j] = data[indexRow2][j];
            data[indexRow2][j] = tmp;
        }
    }

    private int SearchMaxElement(int indexDiag) {
        int resultIndex = indexDiag;
        double maxElement = data[indexDiag][indexDiag];
        for (int i = indexDiag; i < N; i++) {
            if (Math.abs(data[i][indexDiag]) > Math.abs(maxElement)) {
                maxElement = data[i][indexDiag];
                resultIndex = i;
            }
        }
        return resultIndex;
    }

    public void LUP(Matrix L, Matrix U, Matrix P) {
        U.Copy(this);
        P.SetE();

        for (int i = 0; i < N; i++) {
            int iSwap = U.SearchMaxElement(i);
            U.SwapRow(i, iSwap);
            L.SwapRow(i, iSwap);
            P.SwapRow(i, iSwap);

            for (int j = i; j < N; j++) {
                double alpha = U.data[j][i] / U.data[i][i];
                L.data[j][i] = alpha;

                if (i != j) {
                    for (int k = i; k < N; k++) {
                        U.data[j][k] -= U.data[i][k] * alpha;
                    }
                }
            }
        }
    }

    public double DetU(Matrix U) {
        double determinant = 1.0;

        for (int i = 0; i < N; i++) {
            determinant *= U.data[i][i];
        }
        return determinant;
    }

    public int DetP(Matrix P) {
        int count = 0;
        for (int i = 0; i < N; i++) {
            if (P.data[i][i] == 0) count++;
        }
        return count;
    }

    public double Determinant(Matrix U, Matrix P) {
        double result = DetU(U) * (double) Math.pow((-1), DetP(P));
        return result;
    }

    public Matrix ReverseL() {
        Matrix L_reverse = new Matrix(this.N);
        for (int i = 0; i < this.N; i++)
        {
            for (int j = this.N - 1; j >= 0; j--)
            {
                if (i < j) L_reverse.data[i][j] = 0;
                if (i == j) L_reverse.data[i][j] = 1 / this.data[i][j];
                if (i > j)
                {
                    double tempSum = 0;
                    for (int k = j; k <= i - 1; k++)
                    {
                        tempSum += this.data[i][k] * L_reverse.data[k][j];
                    }
                    L_reverse.data[i][j] = (-1) * (tempSum / this.data[i][i]);
                }
            }
        }
        return L_reverse;
    }

    public Matrix ReverseU() {
        Matrix U_reverse = new Matrix(this.N);
        for (int i = this.N - 1; i >= 0; i--)
        {
            for (int j = 0; j < this.N; j++)
            {
                if (i > j) U_reverse.data[i][j] = 0;
                if (i == j) U_reverse.data[i][j] = 1 / this.data[i][j];
                if (i < j)
                {
                    double tempSum = 0;
                    for (int k = j; k >= i + 1; k--)
                    {
                        tempSum += this.data[i][k] * U_reverse.data[k][j];
                    }
                    U_reverse.data[i][j] = (-1) * (tempSum / this.data[i][i]);
                }
            }
        }
        return U_reverse;
    }

    static public Matrix Reverse(Matrix L, Matrix U, Matrix P) {
        Matrix L_reverse = new Matrix(L.N);
        for (int i = 0; i < L.N; i++)
        {
            for (int j = L.N - 1; j >= 0; j--)
            {
                if (i < j) L_reverse.data[i][j] = 0;
                if (i == j) L_reverse.data[i][j] = 1 / L.data[i][j];
                if (i > j)
                {
                    double tempSum = 0;
                    for (int k = j; k <= i - 1; k++)
                    {
                        tempSum += L.data[i][k] * L_reverse.data[k][j];
                    }
                    L_reverse.data[i][j] = (-1) * (tempSum / L.data[i][i]);
                }
            }
        }

        Matrix U_reverse = new Matrix(U.N);
        for (int i = U.N - 1; i >= 0; i--)
        {
            for (int j = 0; j < U.N; j++)
            {
                if (i > j) U_reverse.data[i][j] = 0;
                if (i == j) U_reverse.data[i][j] = 1 / U.data[i][j];
                if (i < j)
                {
                    double tempSum = 0;
                    for (int k = j; k >= i + 1; k--)
                    {
                        tempSum += U.data[i][k] * U_reverse.data[k][j];
                    }
                    U_reverse.data[i][j] = (-1) * (tempSum / U.data[i][i]);
                }
            }
        }
        return U_reverse.Mult(L_reverse).Mult(P);
    }

    public double[] SLAY(Matrix ReverseU, Matrix P, double[] B) {

        double[] matrixY = new double[N];
        double[] matrixX = new double[N];

        Matrix tmp = new Matrix(N);
        tmp = this.Mult(P);

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++){
                matrixY[i] +=  tmp.data[i][j] * B[j];
            }
        }

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++){
                matrixX[i] += ReverseU.data[i][j] * matrixY[j];
            }
        }
        return matrixX;
    }

    public double MethodСonditional () {
        double sum, maximumSum = 0;

        for(int i = 0; i < N; i++)
        {
            sum = 0;
            for(int j = 0; j < N; j++)
            {
                sum += Math.abs(this.data[i][j]);
            }
            if (sum > maximumSum) maximumSum = sum;
        }
        return maximumSum;
    }

    public double[] MethodZeidel(double eps){


        double[] previousVariableValues = new double[N];

        while(true){
            double[] currentVariableValues = new double[N];
            for (int i = 0; i < N; i++) {
                // Инициализируем i-ую неизвестную значением
                // свободного члена i-ой строки матрицы
                currentVariableValues[i] = matrixB[i];

                // Вычитаем сумму по всем отличным от i-ой неизвестным
                for (int j = 0; j < N; j++) {
                    // При j < i можем использовать уже посчитанные
                    // на этой итерации значения неизвестных
                    if (j < i) {
                        currentVariableValues[i] -= this.data[i][j] * currentVariableValues[j];
                    }

                    // При j > i используем значения с прошлой итерации
                    if (j > i) {
                        currentVariableValues[i] -= this.data[i][j] * previousVariableValues[j];
                    }
                }

                // Делим на коэффициент при i-ой неизвестной
                currentVariableValues[i] /= this.data[i][i];
            }

            double error = 0.0;

            for (int i = 0; i < N; i++) {
                error += Math.abs(currentVariableValues[i] - previousVariableValues[i]);
            }
            if (error < eps) {
                break;
            }

            previousVariableValues = currentVariableValues;
        }

        return previousVariableValues;
    }

    public double[] MethodYacoby(double eps) {

        double[] F = new double[N];
        double[] B = new double[N];

        for (int i = 0; i < N; i++){
            F[i] = matrixB[i];
            B[i] = matrixB[i];
        }

        double[] temp = new double[N];
        double norm;

        do {

            for (int i = 0; i < N; i++) {

                temp[i] =- F[i];

                for (int g = 0; g < N; g++) {

                    if (i != g)

                        temp[i] += this.data[i][g] * B[g];
                }

                temp[i] /= -this.data[i][i];
            }

            norm = Math.abs(B[0] - temp[0]);

            for (int h = 0; h < N; h++) {

                if (Math.abs(B[h] - temp[h]) > norm)

                    norm = Math.abs(B[h] - temp[h]);

                B[h] = temp[h];

            }

        } while (norm > eps);

        return B;
    }

}
