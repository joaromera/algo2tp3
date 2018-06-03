//
// Created by gdperez on 01/11/17.
//

#ifndef TP3_AED2_2C2017_H
#define TP3_AED2_2C2017_H

#include <vector>
#include <list>
#include <set>
#include <ostream>

using std::vector;
using std::list;
using std::ostream;

/////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 1: BUCKET SORT GENERICO

/** @brief Generador de buckets en base a un conjunto de valores de input
 *
 *  El iterator::value_type (o sea, el tipo de los valores iterados) debe definir lo siguiente:
 *   - operator int ()                    : devuelve un valor int asociado al value_type
 *   - bool operator<(const value_type &) : comparación entre dos valores
 *
 *  El tipo bucket deber tener operaciones similares a std::set :
 *    - insert(iterator pos, const value_type &)
 *    - iterator begin(), end()
 *    - count(const value_type&)
 *  Dependiendo del tipo de bucket elegido se puede conseguir buckets ordenados o no
 *
 *  @tparam iterator   tipo de los iteradores de entrada
 *  @tparam bucket     tipo del contenedor usado para almacenar un bucket
 *  @param input_begin iterador apuntando al comenizo de los datos de entrada
 *  @param input_end   iterador apuntando más allá del final de los datos de entrada
 */
template <typename iterator, typename bucket>
vector<bucket> generar_buckets(iterator input_begin, iterator input_end);

/**
 * @brief aplanar_buckets genera un vector volcando el contenido de los buckets recibidos
 *
 * Los buckets se deben volcar en el orden definido por el contenedor bucket (begin() ... end() )
 *
 * @tparam bucket tipo del contenedor usado para almacenar un bucket
 * @param B vector de buckets a aplanar
 * @return un vector de bucket::value_type con los valores aplanados
 */
template <typename bucket>
vector<typename bucket::value_type> aplanar_buckets(const std::vector<bucket> & B);


// Bucket sort simple sobre listas
template <typename iterador>
vector<typename iterador::value_type> bucket_sort_list(iterador desde, iterador hasta) {
    return aplanar_buckets(generar_buckets<iterador, list<typename iterador::value_type> >(desde, hasta));
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 2: ORDENAMIENTO DE BILLETES POR PROBABILIDAD DE SER FALSOS

#include <limits>

// Si un billete es falso, esta constante representa su probabilidad (el número más grande posible)
static const unsigned long long probabilidad_max = std::numeric_limits<unsigned long long>::max();

struct billete {
    // Un billete se construye con su numero de serie y se asume falso
    billete(unsigned long long serie) : numero_de_serie(serie), probabilidad_falso(probabilidad_max) {}
    // Constructor alternativo para especifiar un valor de probabilidad
    billete(unsigned long long serie, long prob) : numero_de_serie(serie), probabilidad_falso(prob) {}

    /// \brief operator int devuelve el int correspondiente al año de emisión de un billete
    operator int() const { return numero_de_serie % 10000; }

    /// \brief operator < compara dos billetes por su probabilidad de ser falsos
    ///  Se desempata por numero de serie
    bool operator<(const billete& otro) const { return probabilidad_falso < otro.probabilidad_falso ||
                (probabilidad_falso == otro.probabilidad_falso && numero_de_serie < otro.numero_de_serie) ; }

    /// \brief operator == compara numeros de serie
    bool operator==(const billete& otro) const { return numero_de_serie == otro.numero_de_serie; }

    /// Formato del numero de serie: NNNNNNNNNNAAAA
    /// donde AAAA es el año de emisión y NNNNNNNNNN es el número de secuencia dentro del mismo año
    unsigned long long numero_de_serie;

    /// La probabilidad de ser falso de un billete se define como la cantidad de billetes falsos
    /// conocidos para ese año.
    /// Hay al menos un billete falso emitido por año
    unsigned long long probabilidad_falso;

    friend ostream & operator<<(ostream& os, const billete& b) {
        os << b.numero_de_serie << "(";
        if (b.probabilidad_falso == probabilidad_max) os << "(FALSO)";
        else os << "(p=" << b.probabilidad_falso << ")";
        return os;
    }
};
typedef vector<billete> fajo;

/**
 * @brief ordenar_por_probabilidad ordena billetes según su probabilidad de ser falsos
 * Se sabe que se imprimen M billetes en total por año
 *
 * Todos los billetes de ambos fajos tienen inicialmente probabilidad = probabilidad_max
 *
 * NOTA: "inline" es para poder definir la función en un .h
 *
 * @param falsos_conocidos lista de F billetes falsos ya identificados por el banco
 * @param a_ordenar        lista de N billetes a clasificar, asignar probabilidad, y ordenar
 * @return un vector con los a_ordenar en orden decreciente de probabilidad de falsedad
 *
 * \complexity: O(F log(M) + N log(M*N))
 */
inline fajo ordenar_por_probabilidad(const fajo& falsos_conocidos, const fajo & a_ordenar);


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/// EJERCICIO 3: MULTIPLICACIÓN DE MATRICES POR DIVIDE & CONQUER USANDO STRASSEN

/**
 * Una matriz es un vector de vectores donde
 * - el vector externo almacena las filas
 * - cada vector interno almacena las columnas de cada fila
 */
typedef vector<vector<double> > Matriz;

inline Matriz crear(size_t N, double valor) {
    return Matriz(N, vector<double>(N, valor));
}

/**
 * @brief multiplicación tradicional de matrices A y B iterativa
 *
 * @param A
 * @param B
 * @return el producto de la matriz
 *
 * @complexity O(N^3)
 */
inline Matriz multiplicar(const Matriz& A, const Matriz& B) {
    Matriz C(A.size(), vector<double>(A.size(),0.0) );
    for (size_t i=0; i < A.size(); ++i)
        for (size_t j=0; j < A.size(); ++j)
            for (size_t k=0; k < B.size(); ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

/**
 * Multiplicación de matrices por Divide & Conquer.
 * El Divide & Conquer naif tiene complejidad O(N³).
 * Implementar una versión con complejidad estrictamente mejor.
 * El algoritmo de Strassen plantea un producto inteligente de submatrices que mejora la complejidad a
 * O(N^2.83). Parece poco pero es diferencia!
 *
 * SE PUEDE ASUMIR QUE N ES POTENCIA DE 2
 *
 * @param M1 matriz de N filas por N columnas
 * @param M2 matriz de N filas por N columnas
 * @param K mínimo tamaño para resolver por multiplicación tradicional
 * @return una matriz de N x N con el producto de M1 x M2
 */
inline Matriz multiplicar_strassen(const Matriz& M1, const Matriz& M2, int K);

// TODO: Escribir implementaciones en tp3_impl.h
#include "tp3_impl.h"

#endif //TP3_AED2_2C2017_H
