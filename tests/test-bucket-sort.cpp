//
// Created by gdperez on 01/11/17.
//

#include <gtest/gtest.h>
#include "../src/tp3.h"
#include <string>
#include <vector>
#include <list>
#include <set>

using namespace std;

TEST(bucket_sort, ordenar_numeros) {
    typedef vector<int>::const_iterator iter_t;
    typedef list<int> bucket_t;

    // Prueba básica
    vector<int> v1 = { 1, 2, 3, 4, 5 }, v2 = {5, 4, 3, 2, 1}, v3 = { 3, 2, 1, 4, 5 };

    // Ordenado    
    EXPECT_EQ(bucket_sort_list(v1.begin(), v1.end()), vector<int>({ 1, 2, 3, 4, 5 }) );

    // Ordenado al revés
    auto v = bucket_sort_list(v2.begin(), v2.end());
    EXPECT_EQ(v, v1);

    // Desordenado
    v = bucket_sort_list(v3.begin(), v3.end());
    EXPECT_EQ(v, v1);

    // Negativos, prueba mas chica
    vector<int> v4 = {1,-1,0};
    v = bucket_sort_list(v4.begin(), v4.end());
    EXPECT_EQ(v, vector<int>({-1,0,1}) );

    // Numeros no consecutivos y negativos
    vector<int> v5 = {-2, 111, -10, 25, 33};
    v = bucket_sort_list(v5.begin(), v5.end());
    EXPECT_EQ(v, vector<int>({-10,-2,25,33,111}) );
}

TEST(bucket_sort, generar_buckets_lista) {
    typedef vector<int>::const_iterator iter_t;
    typedef list<int> bucket_t;

    vector<int> v = {1,2,2,3,4,5,5,5,6,8};

    auto B = generar_buckets<iter_t,bucket_t>(v.begin(), v.end());
    EXPECT_EQ(B.size(), 8);
    EXPECT_EQ(B[0].size(),1);
    EXPECT_EQ(B[0].front(),1);
    EXPECT_EQ(B[4].size(),3);
    EXPECT_EQ(B[4].front(),5);

    // Numeros no consecutivos y negativos
    vector<int> v4 = {-2, 111, -10, 25, 33};
    auto B2 = generar_buckets<iter_t,bucket_t>(v4.begin(), v4.end());
    EXPECT_EQ(B2.size(), 122);
    EXPECT_EQ(B2[0].size(),1);
    EXPECT_EQ(B2[0].front(),-10);
    EXPECT_EQ(B2[8].size(),1);
    EXPECT_EQ(B2[8].front(),-2);
    EXPECT_EQ(B2[121].size(),1);
    EXPECT_EQ(B2[121].front(),111);
}

struct alumno {
    string nombre;
    int anio_lu, nro_lu;

    // Agrupo por año
    operator int () const { return anio_lu; }

    // Ordeno empates por nombre
    bool operator <(const alumno& otro) const {
        return anio_lu < otro.anio_lu || (anio_lu == otro.anio_lu && nombre < otro.nombre);
    }
    // comparación de ambas coordenadas
    bool operator ==(const alumno& otro) const {
        return anio_lu == otro.anio_lu && nombre == otro.nombre && nro_lu == otro.nro_lu;
    }
};

TEST(bucket_sort, ordenar_tuplas_lista) {

    // Vector vacio
    vector<alumno> vacio;
    bucket_sort_list(vacio.begin(), vacio.end());
    EXPECT_EQ(vacio, vector<alumno>());

    // un sólo bucket con un repetido (debe mantener orden relativo)
    vector<alumno> a1 = { { "pedro", 2017, 223},
                         { "pedro", 2017, 123},
                        { "juan",  2016, 456} };
    auto a = bucket_sort_list(a1.begin(), a1.end());
    EXPECT_EQ(a, vector<alumno>({ { "juan",  2016, 456}, { "pedro", 2017, 223}, { "pedro", 2017, 123}}) );

    // varios buckets de un elemento c/u
    vector<alumno> a2 = { { "juan",  2015, 456},
                         { "pedro", 2018, 223},
                         { "pablo", 2017, 123}};
    a = bucket_sort_list(a2.begin(), a2.end());
    EXPECT_EQ(a, vector<alumno>({ { "juan",  2015, 456}, { "pablo", 2017, 123}, { "pedro", 2018, 223}}) );

    // varios buckets de varios elementos c/u, con cantidades distintas cada bucket
    vector<alumno> a3 = { { "juan",  2018, 456},
                          { "pedro", 2015, 223},
                          { "pablo", 2017, 456},
                          { "julio", 2015, 223},
                          { "mario", 2017, 123},
                          { "pablo", 2017, 123}};

    a = bucket_sort_list(a3.begin(), a3.end());
    EXPECT_EQ(a, vector<alumno>({
                                   { "pedro", 2015, 223},
                                   { "julio", 2015, 223},
                                   { "pablo", 2017, 456},
                                   { "mario", 2017, 123},
                                   { "pablo", 2017, 123},
                                   { "juan",  2018, 456}
                                 }));
}

TEST(bucket_sort, generar_buckets_set) {

    // varios buckets de varios elementos c/u, con cantidades distintas cada bucket
    vector<alumno> a3 = { { "juan",  2018, 456},
                          { "pedro", 2015, 223},
                          { "pablo", 2017, 456},
                          { "julio", 2015, 223},
                          { "mario", 2017, 123},
                          { "pablo", 2017, 123}};

    // Uso std::set como bucket: cada bucket va a resultar ordenado
    auto B = generar_buckets<vector<alumno>::const_iterator, set<alumno> >(a3.begin(), a3.end());
    EXPECT_EQ(B.size(), 4);
    EXPECT_EQ(*(B[0].begin()), alumno({"julio", 2015, 223}));
    EXPECT_EQ(B[1].size(), 0);
    EXPECT_EQ(*(B[2].begin()), alumno({"mario", 2017, 123}));

}
