extern crate libc;

use std::io;
use std::os::raw::{c_double, c_int};

#[link(name = "determinant", kind = "dylib")]
extern "C" {
    fn determinant(matrix: *mut *mut c_double, n: c_int) -> c_double;
}

fn main() {
    // Вводим размерность матрицы с клавиатуры
    println!("Введите размерность системы:");
    let mut dimension = String::new();
    io::stdin()
        .read_line(&mut dimension)
        .expect("Не удалось прочитать строку");
    let dimension: usize = dimension.trim().parse().expect("Некорректное число");

    // Создаем векторы для хранения коэффициентов матрицы и свободных членов
    let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; dimension]; dimension];
    let mut constants: Vec<f64> = vec![0.0; dimension];

    // Вводим коэффициенты матрицы и свободные члены с клавиатуры
    println!("Введите коэффициенты матрицы:");
    for i in 0..dimension {
        for j in 0..dimension {
            let mut coefficient = String::new();
            io::stdin()
                .read_line(&mut coefficient)
                .expect("Не удалось прочитать строку");
            let coefficient: f64 = coefficient.trim().parse().expect("Некорректное число");
            matrix[i][j] = coefficient;
        }
    }

    println!("Введите свободные члены:");
    for i in 0..dimension {
        let mut constant = String::new();
        io::stdin()
            .read_line(&mut constant)
            .expect("Не удалось прочитать строку");
        let constant: f64 = constant.trim().parse().expect("Некорректное число");
        constants[i] = constant;
    }

    // Создаем C-совместимую матрицу, чтобы передать ее C-функции
    let mut c_matrix: Vec<Vec<c_double>> = matrix
        .iter()
        .map(|row| row.iter().copied().map(|x| x as c_double).collect())
        .collect();
    let c_matrix_ptr: *mut *mut c_double = c_matrix
        .iter_mut()
        .map(|row| row.as_mut_ptr())
        .collect::<Vec<_>>()
        .as_mut_ptr();

        //println!("{:?}", matrix);

    // Вычисляем определитель с помощью C-функции
    let matrix_determinant = unsafe { determinant(c_matrix_ptr, dimension as c_int) };

    // Проверяем, что определитель не равен нулю (иначе система не имеет единственного решения)
    if matrix_determinant.abs() < 1e-10 {
        println!("Система уравнений не имеет единственного решения");
        return;
    }

    // Создаем вектор для хранения решений
    let mut solutions: Vec<f64> = vec![0.0; dimension];

    // Решаем систему уравнений методом Крамера
    for i in 0..dimension {
        let mut modified_matrix = matrix.clone();
        
        // Заменяем i-ый столбец на вектор свободных членов
        for j in 0..dimension {
            modified_matrix[j][i] = constants[j];
        }
        
        let c_matrix_i_ptr: *mut *mut c_double = modified_matrix
            .iter_mut()
            .map(|row| row.as_mut_ptr())
            .collect::<Vec<_>>()
            .as_mut_ptr();
        
        // Вычисляем определитель измененной матрицы
        let determinant_i = unsafe { determinant(c_matrix_i_ptr, dimension as c_int) };
        
        // Вычисляем решение и добавляем его в вектор решений
        solutions[i] = determinant_i / matrix_determinant;
    }

    // Выводим решение
    println!("Решение системы уравнений:");
    for (i, solution) in solutions.iter().enumerate() {
        println!("x{} = {}", i + 1, solution);
    }
}
