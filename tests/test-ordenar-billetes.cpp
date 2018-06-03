#include <gtest/gtest.h>
#include "../src/tp3.h"
#include <chrono>
#include <random>

fajo ordenar_billetes_naif(const fajo & falsos, const fajo & a_ordenar) {
    // conteo de cantidad de falsos por año
    std::map<int, int> contador;
    for (auto & b: falsos)
        contador[int(b)]++;

    fajo resultado;
    resultado.reserve(a_ordenar.size());
    // clasifico billetes a ordenar
    for (auto & b : a_ordenar)
        if (std::find_if(falsos.begin(), falsos.end(), [&](const billete& b2) {return b == b2;} ) != falsos.end()) {
            resultado.push_back(b);
        } else {
            resultado.push_back(billete(b.numero_de_serie, contador[int(b)]));
        }

    std::sort(resultado.rbegin(), resultado.rend());

    return resultado;
}

TEST(ordenar_billetes, test_basico) {
    // 2014: 3 falsos
    // 2015: 1 falso
    // 2016: 2 falsos
    // 2017: 2 falsos
    fajo falsos = { billete(1012014),
                    billete(1012015),
                    billete(1012016),
                    billete(1012017),
                    billete(1052014),
                    billete(1072014),
                    billete(1022016),
                    billete(1152017) };

    // Billetes a ordenar, con alguno falso en medio
    fajo ordenar= {
                    billete(1022015),
                    billete(1032017),
                    billete(1062014),
                    billete(1012014),  // falso
                    billete(1072014),  // falso
                    billete(1032016),
                    billete(1012016),  //falso
                    billete(1112017) };

    auto v = ordenar_por_probabilidad(falsos, ordenar);

    // orden esperado:  falsos, 2014, 2016/2017 (ordenados por nro de serie), 2015
    EXPECT_EQ(v, fajo   ({ billete(1072014),    //falso
                           billete(1012016),    //falso
                           billete(1012014),    //falso
                           billete(1062014),    // p = 3
                           billete(1112017),    // p = 2
                           billete(1032017),    // p = 2
                           billete(1032016),    // p = 2
                           billete(1022015)     // p = 1
                           }));

}

TEST(ordenar_billetes, test_billetera_grande) {
    fajo todos, falsos;
    todos.reserve(26*100000);
    int count_falsos = 0;

    /**
      * Genero 100k Billetes cada año entre 1991 y 2016 (260k billetes total)
      * - Un billete será falso si su nro de serie es múltiplo de (año - 1979)
      * Esto hace que los billetes falsos mas viejos sean mas frecuentes:
      * 1991: multiplos de 12 falsos
      * 1992: todos los multiplos de 13 son falsos
      * ...
      * 2016: Todos los múltiplos de 37 son falsos
      */

    for (int anio = 1991; anio < 2017; ++anio) {
        for (int nserie = 10000; nserie <= 19999; ++nserie) {
            todos.push_back(billete(anio + 10000 * nserie));
            // Reparto los billetes falsos de manera decreciente en los años
            if (nserie % (anio-1989) == 0) {
                ++count_falsos;
                falsos.push_back(billete(anio + 10000 * nserie));
            }
        }
    }
    auto t0 = std::chrono::system_clock::now();
    fajo res = ordenar_por_probabilidad(falsos, todos);
    auto t1 = std::chrono::system_clock::now();

    EXPECT_EQ(res.front(), billete(199991996));
    EXPECT_EQ(res[count_falsos-1].probabilidad_falso, probabilidad_max);
    EXPECT_LE(res[count_falsos].probabilidad_falso, probabilidad_max);
    EXPECT_EQ(res[res.size()/2], billete(168912003));
    EXPECT_EQ(res.back(), billete(100002016));

    std::cout << "T0: " << (t1-t0).count() << std::endl;

    auto t2 = std::chrono::system_clock::now();
    auto res2 = ordenar_billetes_naif(falsos, todos);
    auto t3 = std::chrono::system_clock::now();

    std::cout << "T1: " << (t3-t2).count() << std::endl;
    auto d = double((t1-t0).count()) / (t3-t2).count();
    std::cout << "T0/T1: " << d << std::endl;

    EXPECT_EQ(res, res2);
    EXPECT_LE(d, 0.2);
}

TEST(ordenar_billetes, test_pocos_anios_muchos_falsos) {
    fajo todos, falsos;
    todos.reserve(10000);
    int count_falsos = 0;

    /**
      * Genero 10000 Billetes en 1991
      * - Un billete será falso si su nro de serie es múltiplo de 2
      */

    for (int anio = 1991; anio < 1992; ++anio) {
        for (int nserie = 10000; nserie <= 19999; ++nserie) {
            todos.push_back(billete(anio + 10000 * nserie));
            if (nserie % 2 == 0) {
                ++count_falsos;
                falsos.push_back(billete(anio + 10000 * nserie));
            }
        }
    }
    // Cuadruplico el set de falsos con billetes que no están en "todos"
    int j = falsos.size();
    for (int i = 0; i < j; ++i) {
        falsos.push_back(billete(falsos[i].numero_de_serie+100000000000));
        falsos.push_back(billete(falsos[i].numero_de_serie+200000000000));
        falsos.push_back(billete(falsos[i].numero_de_serie+300000000000));
    }
    auto t0 = std::chrono::system_clock::now();
    fajo res = ordenar_por_probabilidad(falsos, todos);
    auto t1 = std::chrono::system_clock::now();

    EXPECT_EQ(res.front(), billete(199981991));
    EXPECT_EQ(res[count_falsos-1].probabilidad_falso, probabilidad_max);
    EXPECT_LE(res[count_falsos].probabilidad_falso, probabilidad_max);
    EXPECT_EQ(res[res.size()/2], billete(199991991));
    EXPECT_EQ(res.back(), billete(100011991));

    std::cout << "T0: " << (t1-t0).count() << std::endl;

    auto t2 = std::chrono::system_clock::now();
    auto res2 = ordenar_billetes_naif(falsos, todos);
    auto t3 = std::chrono::system_clock::now();

    std::cout << "T1: " << (t3-t2).count() << std::endl;
    auto d = double((t1-t0).count()) / (t3-t2).count();
    std::cout << "T0/T1: " << d << std::endl;

    EXPECT_EQ(res, res2);
    EXPECT_LE(d, 0.2);
}

TEST(ordenar_billetes, test_muchos_anios_pocos_falsos) {
    fajo todos, falsos;
    todos.reserve(2015000);
    int count_falsos = 0;

    /**
      * Genero 1000 Billetes cada año entre 1 y 2016 (2M billetes total)
      * - Un billete será falso si su nro de serie es 1101
      */

    for (int anio = 1; anio < 2017; ++anio) {
        for (int nserie = 1000; nserie <= 1999; ++nserie) {
            todos.push_back(billete(anio + 10000 * nserie));
            // Reparto los billetes falsos de manera decreciente en los años
            if (nserie == 1101) {
                ++count_falsos;
                falsos.push_back(billete(anio + 10000 * nserie));
            }
        }
    }
    auto t0 = std::chrono::system_clock::now();
    fajo res = ordenar_por_probabilidad(falsos, todos);
    auto t1 = std::chrono::system_clock::now();

    EXPECT_EQ(res.front(), billete(11012016));
    EXPECT_EQ(res[count_falsos-1].probabilidad_falso, probabilidad_max);
    EXPECT_LE(res[count_falsos].probabilidad_falso, probabilidad_max);
    EXPECT_EQ(res[res.size()/2], billete(15002016));
    EXPECT_EQ(res.back(), billete(10000001));

    std::cout << "T0: " << (t1-t0).count() << std::endl;

    auto t2 = std::chrono::system_clock::now();
    auto res2 = ordenar_billetes_naif(falsos, todos);
    auto t3 = std::chrono::system_clock::now();

    std::cout << "T1: " << (t3-t2).count() << std::endl;
    auto d = double((t1-t0).count()) / (t3-t2).count();
    std::cout << "T0/T1: " << d << std::endl;

    EXPECT_EQ(res, res2);
    EXPECT_LE(d, 0.25);
}
