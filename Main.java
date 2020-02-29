package com.metanit;

import java.io.FileReader;
import java.text.NumberFormat;
import java.util.Scanner;
import java.util.Random;
//import java.util.Matr--;


public class Main {

    public static void main(String[] args) throws Exception {


        FileReader reader = new FileReader("C:\\Users\\honey\\IdeaProjects\\matrix\\src\\com\\metanit\\text.txt");
        Scanner scan = new Scanner(reader);

        String str = new String(scan.nextLine());
        System.out.println(str);

        NumberFormat nf = NumberFormat.getNumberInstance();
        nf.setMaximumFractionDigits(2);
        nf.setGroupingUsed(false);

        int i = 0, count = 0;
        for (String retval : str.split(" ")) {
            count++;
        }
        double[] array = new double[count];

        for (String retval : str.split(" ")) {
            array[i] = Double.parseDouble(retval);
            i++;
        }

        //Объявление матриц
        int N = (int)array[0];
        double q = array[1];
        Matrix matrixA = new Matrix(N);
        matrixA.Generate(q);

        matrixA.Print();
        System.out.println();

        //Раздложение матрицы
        Matrix matrixU = new Matrix(N);
        Matrix matrixL = new Matrix(N);
        Matrix matrixP = new Matrix(N);

        matrixA.LUP(matrixL, matrixU, matrixP);

        System.out.println("U: ");
        matrixU.Print();

        System.out.println();
        System.out.println("L: ");
        matrixL.Print();

        System.out.println();
        System.out.println("P: ");
        matrixP.Print();
        System.out.println();


        //Результирующая матрица
        System.out.println("P * A: ");
        matrixP.Mult(matrixA).Print();
        System.out.println();
        System.out.println("L * U: ");
        matrixL.Mult(matrixU).Print();

        double matrixDeterminant = matrixA.Determinant(matrixU, matrixP);

        double detU = matrixP.DetU(matrixU);
        System.out.println();
        System.out.println("Determinant U: ");
        System.out.println(detU);

        System.out.println();
        System.out.println("Determinant A: ");
        System.out.println(matrixDeterminant);


        Matrix reverseA = matrixA.Reverse(matrixL, matrixU, matrixP);
        System.out.println();
        System.out.println("Reverse A: ");
        reverseA.Print();

        matrixA.Mult(reverseA).Print();
        reverseA.Mult(matrixA).Print();

        Matrix reverseL = matrixL.ReverseL();
        Matrix reverseU = matrixU.ReverseU();

        double[] X = new double[N];
        X = reverseL.SLAY(reverseU, matrixP, matrixA.matrixB);
        matrixA.PrintArray(X);

        double normA = matrixA.MethodСonditional();
        double normReverseA = reverseA.MethodСonditional();

        double digitCond = normA * normReverseA;
        System.out.println();
        System.out.println("Norma A: ");
        System.out.println(digitCond);

        double EPS = Math.pow(10, (-6));
        double[] solution = new double[N];
        solution = matrixA.MethodZeidel(EPS);
        System.out.println();
        System.out.println("Methos Zeydely: ");
        matrixA.PrintArray(solution);

        double[] solutionYacoby = new double[N];
        solutionYacoby = matrixA.MethodYacoby(EPS);
        System.out.println();
        System.out.println("Methos Yacoby: ");
        matrixA.PrintArray(solutionYacoby);



    }
}
