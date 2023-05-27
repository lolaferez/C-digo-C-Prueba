// MODELO DE RED NEURONAL USANDO EL ALGORITMO DE METROPOLIS
// MODELO DE HOPFIELD

#include <iostream>
#include <cmath>
#include <fstream>
#include "gsl_rng.h" //Libreria para generación de números aleatorios

gsl_rng *tau;

using namespace std;

#define N 70     // Dimension de la matriz cuadrada
#define T 0.0001 // Temperatura del sistema
#define monte 90 // Numero de pasos montecarlo a realizar
#define pat 1    // Cantidad de patrones para aprender
#define c 1 //Para comenzar la matriz ordenada c=0; aleatoria c=1

int main()
{

    extern gsl_rng *tau;  // Puntero al estado del número aleatorio
    int semilla = 135294; // Semilla del generador de números aleatorios

    ifstream patron1; // Ficheros para lectura de los patrones iniciales
    ofstream fichero; // Fichero para la escritura de datos iterativos

    int i, j, i2, j2, k, q, n, m, a[pat], cuenta;           // Multiples contadores y el valor a caracteristico de cada patron
    int M[N][N], P[N][N][pat];                                 // Matriz de la red neuronal y del patron
    double w[N][N], suma, suma2, tetha[N][N], H, aux, valor, p, eps; // Para la interaccion entre neuronas

    tau = gsl_rng_alloc(gsl_rng_taus); // Inicializamos el puntero
    gsl_rng_set(tau, semilla);         // Inicializamos la semilla

    fichero.open("red.txt");
    patron1.open("patron1.txt");


    if (c == 0) // Llenamos la matriz ordenada
    {
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
                M[i][j] = 1;
        }
    }
    else // Llenamos la matriz con valor de actividad neuronal 1 y 0 aleatoriamente
    {
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                M[i][j] = gsl_rng_uniform_int(tau, 2); // Numero random entero o 0 o 1
            }
        }
    }

    // Guardamos la matriz inicial en el archivo
    for (i = 0; i < N; i++)
    { // en las primeras posiciones del archivo van los datos y dejamos un espacio en blanco entre instante e instante
        for (j = 0; j < N; j++)
        {
            if (j < (N - 1))
            {
                fichero << M[i][j] << ", ";
            }
            else
            {
                fichero << M[i][j] << endl;
            }
        }
    }

    // Cambio de linea para nueva iteracion
    fichero << endl;

    // Ahora, llenamos la matriz con los datos del primer patron deseado
    while (!patron1.eof())
    { // Para comprobar que el fichero contiene datos.
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                patron1 >> P[i][j][0];
            }
        }
    }

    /* RESTO DE PATRONES
    while (!patron2.eof())
    { // Para comprobar que el fichero contiene datos.
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                patron2 >> P[i][j][1];
            }
        }

    }*/

    for (k = 0; k < pat; k++)
    {
        a[k] = 0;
    }

    // Calculo el valor de a del patron k-esimo
    for (k = 0; k < pat; k++)
    {
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                a[k] += P[i][j][k];
            }
        }
    }

    cuenta = 0;
    for (q = 1; q <= (monte * N * N); q++)
    {
        // Elegimos un punto a azar de la red
        n = gsl_rng_uniform_int(tau, N);
        m = gsl_rng_uniform_int(tau, N);

        // Calculamos la energía en ese determinado M[n][m] paso a paso

        // Comenzamos con el calculo del primer termino de la expresion para el incremento de energia

        suma = 0.0;

        for (k = 0; k < pat; k++)
        {
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    if ((i == n) && (j == m))
                    {
                        w[i][j] = 0.0;
                    }
                    else
                    {
                        w[i][j] = 1 / (N * N) * (P[n][m][k] - a[k]) * (P[i][j][k] - a[k]);
                        suma = suma + w[i][j] * M[i][j];
                    }
                }
            }
        }

        suma2 = 0.0;

        // Calculo del segundo termino de la expresion
        for (k = 0; k < pat; k++) // Para cada patron
        {
            for (i = 0; i < N; i++) // Para cada neurona del sistema
            {
                for (j = 0; j < N; j++)
                {
                    tetha[i][j]=0.0;
                    aux=0.0;
                    
                    for (i2 = 0; i2 < N; i2++) //Para calcular el sumatorio de interacciones w de una
                    {                          // determinada neurona [i][j] con las demas
                        for(j2=0;j2<N;j2++)
                        {
                            if ((i2 == i) && (j2 == j))
                            {
                                aux = 0.0;
                            }
                            else
                            {
                                aux = (1.0) / (N * N) * (P[i2][j2][k] - a[k]) * (P[i][j][k] - a[k]);
                                tetha[i][j]=tetha[i][j]+aux;
                            }
                        }
                    }

                    suma2 = suma2 + M[i][j]*tetha[i][j]/(2.0);
                }
            }
        }

        // Ahora si calculamos el  incremento del hamiltoniano del sistema para una determinada configuracion
        if (M[n][m] == 0)
        {
            H = suma2 - suma;
        }
        else
        {
            H = suma - suma2;
        }

        // Una vez tenemos la energia de la configuracion, seguimos con el algoritmo de metropolis
        // Recordemos que no se considera el estado para temperatura 0, ya que esta no es accesible

        valor = exp(-H / T);
        if ((1.0) <= valor)
        {
            p = 1.0;
        }
        else
        {
            p = valor;
        }

        eps = gsl_rng_uniform(tau); // Numero aleatorio real entre 0 y 1

        if ((eps <= p) && (M[n][m] == 0))
        {
            M[n][m] = 1;
        }

        if ((eps <= p) && (M[n][m] == 1))
        {
            M[n][m] = 0;
        }

        // Ahora guardo la matriz cada paso montecarlo
        if (cuenta == (N * N))
        {
            for (i = 0; i <= (N - 1); i++)
            { // en las primeras posiciones del archivo van los datos y dejamos un
              // espacio en blanco entre instante e instante

                if (i < (N - 1))
                {
                    for (j = 0; j <= (N - 1); j++)
                    {
                        if (j < (N - 1))
                        {
                            fichero << M[i][j] << ", ";
                        }
                        else
                            fichero << M[i][j] << endl;
                    }
                }
                else
                {
                    fichero << endl;
                }
            }

            cuenta = 0;
        }

        cuenta = cuenta + 1;
    }

    fichero.close();
    patron1.close();
}