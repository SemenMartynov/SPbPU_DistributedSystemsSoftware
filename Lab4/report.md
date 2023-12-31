# Сравнение оптимизаций

## Часть 1. Выбор подходящего уровня оптимизации

Запуск скрипта тестирования `$ bash test.sh livermorec/livermorec.cpp`

| Ключи оптимизации                     | Размер исп. файла  | Время исполнения  |
|:------------------------------------- |:------------------:| -----------------:|
| -O0                                   | 156K               | 0m20.744s         |
| -Os                                   | 40K                | 0m3.329s          |
| -O1                                   | 44K                | 0m3.646s          |
| -O2                                   | 48K                | 0m3.346s          |
| -O3                                   | 56K                | 0m2.976s          |
| -O2 -march=native                     | 48K                | 0m3.191s          |
| -O3 -march=native                     | 56K                | 0m2.955s          |
| -O2 -march=native -funroll-loops      | 64K                | 0m3.127s          |
| -O3 -march=native -funroll-loops      | 68K                | 0m2.853s          |

## Часть 2. Выбор системного метода оптимизации

Оптимальной опцией, с точки зрения наименьшего времени исполнения, является `-O3 -march=native -funroll-loops`

Выбор опций межпроцедурной оптимизации:  
-`-fipa-sra` включена на уровнях `-O2`, `-Os` и `-O3`  
-`-fipa-ra` включена на уровнях `-O2`, `-Os` и `-O3`  
-`-fipa-pure-const` включена с уровня `-O1` и выше  
-`-fipa-reference` включена с уровня `-O1` и выше  
-`-fipa-reference-addressable` включена с уровня `-O1` и выше  
-`-fipa-stack-alignment` включена по умолчанию  
-`-fipa-pta` по умолчанию не включена, подходит для анализа  
-`-fipa-profile` включена с уровня `-O1` и выше  
-`-fipa-modref` включена с уровня `-O1` и выше  
-`-fipa-cp` включена на уровнях `-O2`, `-Os` и `-O3`  
-`-fipa-cp-clone` включена на уровне `-O3`  
-`-fipa-bit-cp` включена на уровне `-O2`  
-`-fipa-vrp` включена на уровне `-O2`  
-`-fipa-icf` включена на уровне `-O2`  
-`-fipa-strict-aliasing` включена на уровне `-O2` и `-O3`  

Выбор опций оптимизацией времени компоновки:  
-`-flto` не включена по умолчанию
> Следует использовать `-flto=auto` (https://www.mail-archive.com/gcc-patches@gcc.gnu.org/msg258441.html)

Запуск скрипта тестирования `$ bash test2.sh livermorec/livermorec.cpp`

| Ключи оптимизации                                     | Размер исп. файла  | Время исполнения  |
|:----------------------------------------------------- |:------------------:| -----------------:|
| -O3 -march=native -funroll-loops -fipa-pta            | 64K                | 0m2.743s          |
| -O3 -march=native -funroll-loops -flto=auto           | 64K                | 0m2.704s          |
| -O3 -march=native -funroll-loops -fipa-pta -flto=auto | 60K                | 0m2.628s          |

Оптимизация с обратной связью
```sh
$ g++ -O3 -march=native -funroll-loops -fprofile-generate livermorec/livermorec.cpp -o test
$ ./test
$ g++ -O3 -march=native -funroll-loops -fprofile-use livermorec/livermorec.cpp -o test
```

| Ключи оптимизации                                | Размер исп. файла  | Время исполнения  |
|:------------------------------------------------ |:------------------:| -----------------:|
| -O3 -march=native -funroll-loops -fprofile-use   | 56K                | 0m2.655s          |

Оптимизация с оптимальной опцией, межпроцедурной оптимизацией, оптимизацией времени компоновки и с оптимизацией с обратной связью
```sh
$ g++ -O3 -march=native -funroll-loops -fipa-pta -flto=auto -fprofile-generate livermorec/livermorec.cpp -o test
$ ./test
$ g++ -O3 -march=native -funroll-loops -fipa-pta -flto=auto -fprofile-use livermorec/livermorec.cpp -o test
```

| Ключи оптимизации                                                     | Размер исп. файла  | Время исполнения  |
|:--------------------------------------------------------------------- |:------------------:| -----------------:|
| -O3 -march=native -funroll-loops -fipa-pta -flto=auto -fprofile-use   | 52K                | 0m2.677s          |

## Выводы

1. Наиболее заметным был прирост производительности и сокращение размера исполняемого файла при переходе от `-O0` к `-Os` (Optimize for size)  
2. Размер файла и время его исполнения не коррелируют друг с другом на 95 процентили: можно отдельно добиться минимального размера или максимального быстродействия  
3. Большинство опций межпроцедурной оптимизации уже встроены в существующий набор стандартных оптимизаций `-O{n}`, нет смысла их указывать дополнительно  
4. Наилучший результат по производительности был получен набором `-O3 -march=native -funroll-loops -fipa-pta -flto=auto`, но это результат для конкретного кода  
5. Оптимизация с обратной связью дала хороший результат, но не лучший (причина может быть в погрешности измерений.)