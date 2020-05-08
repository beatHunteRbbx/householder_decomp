import java.util.*;

public class Main {

    private static Double[][] matrix = new Double[][] {
            {2.0, 1.2, -1.0, 1.0},
            {1.2, 0.5, 2.0, -1.0},
            {-1.0, 2.0, -1.5, 0.2},
            {1.0, -1.0, 0.2, 1.5}
    };
    private static Double[][] E = new Double[][] {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0, 1.0}
    };


    public static void main(String[] args) {
        System.out.println("A0:");
        printMatrix(matrix);
        System.out.println();
        Double[][] tempMatrix = matrix.clone();
        Double[][] answerMatrix = new Double[matrix.length][matrix[0].length];
        for (int i = 0; i < 1000; i++) {
            Double[][] A = doHouseholderDecomp(tempMatrix);
            //вычисляем разность между диагональными элементами прошлой матрицы и текущей матрицы
            Double[][] deltaMatrix = subtractMatrices(tempMatrix, A);

            //запоминаем текущую матрицу для следующего шага
            tempMatrix = A.clone();

            //вывод матриц в консоль
            int numb = i+1;
            int prevNumb = numb-1;
            System.out.println("A" + numb + ":");
            printMatrix(A);
            System.out.println("A"+prevNumb+" - "+"A"+numb+":");
            printMatrix(deltaMatrix);
            System.out.println();

            //если разность диагональных элементов прошлой матрицы и
            //текущей <= 0.01, то берём текущую матрицу как финальную
            //у которой собственные значения находятся на её диагонали
            if (isMatrixValid(deltaMatrix)) {
                answerMatrix = A.clone();
                break;
            }
        }

        ArrayList<Double> matrixLambdaList = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++) {
            matrixLambdaList.add(answerMatrix[i][i]);
        }

        System.out.println("Собственные значения матрицы A0:");
        int number = 1;
        for (Double lambda : matrixLambdaList) {
            System.out.println(number++ + ": " + lambda);
        }
    }

    private static Double[][] doHouseholderDecomp(Double[][] matrix) {
        //----------------Шаг 1: обнуляются элементы 1ого стобца под диагональным элементом а11-------------------------
        Double[][] vector = new Double[matrix.length][1];
        double sum = 0;
        //формируем вектор V
        for (int i = 0; i < matrix.length; i++) sum += Math.pow(matrix[i][0], 2);
        vector[0][0] = matrix[0][0] + Math.pow(sum, 0.5);
        vector[1][0] = matrix[1][0];
        vector[2][0] = matrix[2][0];
        vector[3][0] = matrix[3][0];

        //Вычисление H1
        Double[][] H1 = calculateH(vector);

        //Вычисление А1
        Double[][] A1 = multiplyMatrices(H1, matrix);
        //--------------------------------------------Конец Шага 1------------------------------------------------------

        //----------------Шаг 2: обнуляются элементы 2ого стобца под диагональным элементом а22-------------------------
        //формируем вектор V
        sum = 0;
        for (int i = 1; i < matrix.length; i++) sum += Math.pow(matrix[i][1], 2);
        vector[0][0] = 0.0;
        vector[1][0] = matrix[1][1] + Math.pow(sum, 0.5);
        vector[2][0] = matrix[2][0];
        vector[3][0] = matrix[3][0];

        //Вычисление H2
        Double[][] H2 = calculateH(vector);

        //Вычисление А2
        Double[][] A2 = multiplyMatrices(H2, A1);
        //--------------------------------------------Конец Шага 2------------------------------------------------------

        //----------------Шаг 3: обнуляются элементы 3ого стобца под диагональным элементом а33-------------------------
        //формируем вектор V
        sum = 0;
        for (int i = 2; i < matrix.length; i++) sum += Math.pow(matrix[i][1], 2);
        vector[0][0] = 0.0;
        vector[1][0] = 0.0;
        vector[2][0] = matrix[2][2] - Math.pow(sum, 0.5);
        vector[3][0] = matrix[3][0];

        //Вычисление H3
        Double[][] H3 = calculateH(vector);

        //Вычисление А3
        Double[][] A3 = multiplyMatrices(H3, A2);
        //--------------------------------------------Конец Шага 3------------------------------------------------------

        Double[][] Q = multiplyMatrices(multiplyMatrices(H1, H2), H3);
        Double[][] R = A3;

        return multiplyMatrices(R, Q);
    }

    //Вычисление H1
    private static Double[][] calculateH(Double[][] vector) {
        Double[][] tVector = transpose(vector);
        Double[][] numerator = multiplyMatrices(vector, tVector);       //числитель дроби
        Double[][] denominator = multiplyMatrices(tVector, vector);     //знаменатель дроби

        Double[][] fraction = multiplyMatrixOnNumb(numerator, 1/denominator[0][0]); //дробь
        Double[][] H = subtractMatrices(E, multiplyMatrixOnNumb(fraction, 2));
        return H;
    }
    private static Double[][] transpose(Double[][] array) {
        Double[][] transposed = new Double[array[0].length][array.length];
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                transposed[j][i] = array[i][j];
            }
        }
        return transposed;
    }


    private static void printMatrix(Double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            System.out.print("(\t");
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.printf("%.2f\t", matrix[i][j]);
            }
            System.out.println(")");
        }
    }

    private static Double[][] multiplyMatrixOnNumb(Double[][] matrix, double numb) {
        Double[][] resultMatrix = new Double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                resultMatrix[i][j] = matrix[i][j] * numb;
            }
        }
        return resultMatrix;
    }

    private static Double[][] multiplyMatrices(Double[][] matrix1, Double[][] matrix2) {
        if (matrix1[0].length != matrix2.length) {
            throw new IllegalArgumentException(
                    "Такие матрицы нельзя перемножить, так как количество столбцов матрицы matrix1 не равно " +
                    "количеству строк матрицы matrix2."
            );
        }
        else {
            Double[][] resultMatrix = new Double[matrix1.length][matrix2[0].length];
            for (int i = 0; i < matrix1.length; i++) {
                for (int j = 0; j < matrix2[0].length; j++) {
                    resultMatrix[i][j] = 0.0;
                }
            }

            for (int i = 0; i < matrix1.length; i++) {
                for (int j = 0; j < matrix2[0].length; j++) {
                    for (int k = 0; k < matrix2.length; k++) {
                        resultMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
                    }
                }
            }
            return resultMatrix;
        }
    }

    private static Double[][] subtractMatrices(Double[][] matrix1, Double[][] matrix2) {
        if (matrix1.length != matrix2.length || matrix1[0].length != matrix2[0].length) {
            throw new IllegalArgumentException(
                    "Такие матрицы нельзя вычесть, так как matrix1 и matrix2 имеют не одинакове " +
                            "количество строк и столбцов"
            );
        }
        else {
            Double[][] resultMatrix = new Double[matrix1.length][matrix2[0].length];
            for (int i = 0; i < matrix1.length; i++) {
                for (int j = 0; j < matrix2[0].length; j++) {
                    resultMatrix[i][j] = matrix1[i][j] - matrix2[i][j];
                }
            }
            return resultMatrix;
        }
    }

    private static boolean isMatrixValid(Double[][] matrix) {
        boolean isValid = false;
        if (    Math.abs(matrix[0][0]) < 0.01 && Math.abs(matrix[1][1]) < 0.01 &&
                Math.abs(matrix[2][2]) < 0.01 && Math.abs(matrix[3][3]) < 0.01) isValid = true;
        return isValid;
    }
}
