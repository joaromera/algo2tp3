#ifndef TP3_IMPL_H
#define TP3_IMPL_H

#include "tp3.h"

#include <limits>
#include <algorithm>
#include <map>

using std::map;

///////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 1
///////////////////////////////////////////////////////////////////////////////

template <typename iterator, typename bucket>
vector<bucket> generar_buckets(iterator input_begin, iterator input_end) {
    long minimo = std::numeric_limits<signed long>::max();
    long maximo = std::numeric_limits<signed long>::min();

    for (auto it = input_begin; it != input_end; ++it) {
        if (int(*it) > maximo) {
            maximo = int(*it);
        }
        if (int(*it) < minimo) {
            minimo = int(*it);
        }
    }

    long diferencia = maximo - minimo + 1;
    vector<bucket> buckets(diferencia);

    for (auto it = input_begin; it != input_end; ++it) {
        long pos = int(*it) - minimo;
        buckets[pos].insert(buckets[pos].end(), *it);
    }

    return buckets;
}

template <typename bucket>
vector<typename bucket::value_type> aplanar_buckets(const std::vector<bucket> & B) {
    vector<typename bucket::value_type> res;
    for (auto b : B) {
        for (auto t : b) {
            res.push_back(t);
        }
    }
    return res;
}

///////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 2
///////////////////////////////////////////////////////////////////////////////

fajo ordenar_por_probabilidad(const fajo& falsos_conocidos, const fajo & a_ordenar) {
    auto bucketYear = generar_buckets<fajo::const_iterator, std::set<billete>>(falsos_conocidos.begin(), falsos_conocidos.end());
    long minimo = *bucketYear[0].begin();

    fajo result;
    result.reserve(a_ordenar.size());

    for (auto b : a_ordenar) {
        billete tmp = b;
        if (bucketYear[int(b)-minimo].count(b.numero_de_serie) == 1) {
            tmp.probabilidad_falso = probabilidad_max;
        } else {
            auto prob = bucketYear[int(b)-minimo].size();
            if (prob == 0) prob = 1;
            tmp.probabilidad_falso = prob;
        }
        result.push_back(tmp);
    }

    std::sort(result.begin(), result.end());
    std::reverse(result.begin(), result.end());

    return result;
}

///////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 3
///////////////////////////////////////////////////////////////////////////////

inline Matriz sumarMatriz(const Matriz& A, int i_a, int j_a, const Matriz& B, int i_b, int j_b, int size) {
    Matriz result = crear(size, 0);

    for (int i = i_a; i < size; i++) {
        for (int j = j_a; j < size; j++) {
            result[i][j] = A[i][j] + B[i + i_b][j + j_b];
        }
    }
    return result;
}

inline Matriz restaMatriz(const Matriz& A, int i_a, int j_a, const Matriz& B, int i_b, int j_b, int size) {
    Matriz result = crear(size, 0);

    for (int i = i_a; i < size; i++) {
        for (int j = j_a; j < size; j++) {
            result[i][j] = A[i][j] - B[i + i_b][j + j_b];
        }
    }
    return result;
}

inline Matriz reconstruir(const Matriz& A11, const Matriz& A12, const Matriz& A21, const Matriz& A22, int N) {
    int size = N * 2;
    Matriz result = crear(size, 0);

    for (int i = 0; i < size / 2; i++) {
        for (int j = 0; j < size / 2; j++) {
            result[i][j] = A11[i][j];
        }
    }

    for (int i = 0; i < size / 2; i++) {
        for (int j = size / 2; j < size; j++) {
            result[i][j] = A12[i][j - (size / 2)];
        }
    }

    for (int i = size / 2; i < size; i++) {
        for (int j = 0; j < size / 2; j++) {
            result[i][j] = A21[i - (size / 2)][j];
        }
    }

    for (int i = size / 2; i < size; i++) {
        for (int j = size / 2; j < size; j++) {
            result[i][j] = A22[i - (size / 2)][j - (size / 2)];
        }
    }

    return result;
}

inline Matriz submatriz(const Matriz& A, int f, int c, int N) {
    Matriz result = crear(N, 0);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i][j] = A[i + f][j + c];
        }
    }
    return result;
}

inline Matriz multiplicar_strassen(const Matriz& A, const Matriz& B, int K) {
    if (A.size() <= K) return multiplicar(A,B);

    int size = A.size() / 2;

    Matriz a11 = submatriz(A,0,0,size);
    Matriz a12 = submatriz(A,0,size,size);
    Matriz a21 = submatriz(A,size,0,size);
    Matriz a22 = submatriz(A,size,size,size);

    Matriz b11 = submatriz(B,0,0,size);
    Matriz b12 = submatriz(B,0,size,size);
    Matriz b21 = submatriz(B,size,0,size);
    Matriz b22 = submatriz(B,size,size,size);

    Matriz M1 = multiplicar_strassen(sumarMatriz(a11,0,0,a22,0,0,size), sumarMatriz(b11,0,0,b22,0,0,size), K);
    Matriz M2 = multiplicar_strassen(sumarMatriz(a21,0,0,a22,0,0,size), b11, K);
    Matriz M3 = multiplicar_strassen(a11, restaMatriz(b12,0,0,b22,0,0,size), K);
    Matriz M4 = multiplicar_strassen(a22, restaMatriz(b21,0,0,b11,0,0,size), K);
    Matriz M5 = multiplicar_strassen(sumarMatriz(a11,0,0,a12,0,0,size), b22, K);
    Matriz M6 = multiplicar_strassen(restaMatriz(a21,0,0,a11,0,0,size), sumarMatriz(b11,0,0,b12,0,0,size), K);
    Matriz M7 = multiplicar_strassen(restaMatriz(a12,0,0,a22,0,0,size), sumarMatriz(b21,0,0,b22,0,0,size), K);

    Matriz C11 = sumarMatriz(M1,0,0,M4,0,0,size);
    C11 = restaMatriz(C11,0,0,M5,0,0,size);
    C11 = sumarMatriz(C11,0,0,M7,0,0,size);

    Matriz C21 = sumarMatriz(M2,0,0,M4,0,0,size);

    Matriz C12 = sumarMatriz(M3,0,0,M5,0,0,size);

    Matriz C22 = restaMatriz(M1,0,0,M2,0,0,size);
    C22 = sumarMatriz(C22,0,0,M3,0,0,size);
    C22 = sumarMatriz(C22,0,0,M6,0,0,size);

    Matriz result = reconstruir(C11, C12, C21, C22, size);

    return result;
}

#endif // TP3_IMPL_H
