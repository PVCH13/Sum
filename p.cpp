#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm> 
#include <sys/time.h>

//-------------------------------------------------------------
// Función para obtener el tiempo actual en segundos
// (sirve para medir el tiempo de ejecución del programa)
//-------------------------------------------------------------
double seconds()
{
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  
  return sec;
}
//-------------------------------------------------------------
// Función que llena el centro de la matriz con temperatura 100.0
// Parámetros:
//   Nx, Ny         -> dimensiones de la matriz
//   centro_medios  -> tamaño del cuadrado central
//   current        -> puntero a la matriz actual
//-------------------------------------------------------------
int llenar_centro( int Nx, int Ny, double centro_medios, double *current){
  for(int i = ((Nx/2)-centro_medios);i<(Nx/2)+centro_medios;++i) {
    for(int j =(Ny/2)-centro_medios;j<((Ny/2)+centro_medios);++j) {
    current[(i * Ny) + j] = 100.0;// asigna valor 100 en el centro
    }
  }
  return 0;
}


//-------------------------------------------------------------
// Función que recalcula los valores de temperatura (difusión térmica)
// usando un esquema explícito en diferencias finitas.
// Parámetros:
//   Nx, Ny -> dimensiones de la matriz
//   current -> matriz actual
//   prox    -> matriz siguiente
//   alpha   -> coeficiente de difusión
//   dx, dy  -> tamaño de paso espacial
//   dt      -> paso de tiempo
//-------------------------------------------------------------
int recalcular_matriz( int Nx, int Ny, double *current, double *prox, double alpha, double dx, double dy, double dt){
  // Recorre la matriz (excepto bordes)
  for(int i = 1; i < Nx-1; ++i){
    for(int j = 1; j < Ny-1; ++j){
      prox[(i*Ny) + j] = current[(i * Ny) + j]+
      alpha*dt*((current[((i+1)*Ny)+j]-2*current[(i * Ny) + j]+current[((i-1)*Ny)+j])/(dx*dx)+
      (current[(i * Ny) + j+1]-2*current[(i * Ny) + j]+current[(i * Ny) + j-1])/(dy*dy));
    }
  }
  
  // Condiciones de frontera (bordes reflejan el valor interior)
  for (int j = 0; j < Ny; ++j) {
    prox[0 * Ny + j] = prox[1 * Ny + j];         // borde izquierdo
    prox[(Nx - 1) * Ny + j] = prox[(Nx - 2) * Ny + j];   // borde derecho
  }
  for (int i = 0; i < Nx; ++i) {
    prox[i * Ny + 0] = prox[i * Ny + 1];         // borde inferior
    prox[i * Ny + (Ny - 1)] = prox[i * Ny + (Ny - 2)];   // borde superior
  }
  
  return 0;
}

//-------------------------------------------------------------
// Función que guarda en un archivo los valores de una matriz 2D con sus coordenadas.
// Escribe cada elemento en el archivo "data.dat" en el formato: j  i  valor.
// Se agrega una línea en blanco entre filas para facilitar la visualización.
// Parámetros:
//   Nx, Ny  -> dimensiones de la matriz
//   current -> matriz actual almacenada en un arreglo lineal (orden fila-major)
// Retorna:
//   0 si la operación fue exitosa
//   1 si no se pudo abrir el archivo para escritura
//-------------------------------------------------------------
int guardar_datos( int Nx, int Ny, double *current){
  std::ofstream fout("data.dat");
  if (!fout.is_open()) {
    std::cerr << "Error: no se pudo abrir data.dat para escritura.\n";
    return 1;
  }

  // Escribe cada valor de la matriz con sus coordenadas
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      fout << j << " " << i << " " << current[(i * Ny) + j] << "\n";  // separa por espacios
    }
    fout << "\n";// línea en blanco entre filas (útil para visualización)
  }
  fout.close();
  return 0;
}

//-------------------------------------------------------------
// Función principal que ejecuta la simulación de difusión térmica en 2D.
// Lee las dimensiones de la matriz y el número de iteraciones desde la línea de comandos.
// Inicializa el dominio, calcula parámetros físicos (dx, dy, dt) y realiza la evolución temporal
// aplicando difusión y manteniendo una región caliente en el centro.
// También mide el tiempo total de la simulación y guarda los resultados finales en un archivo.
// Parámetros:
//   argc, argv -> argumentos de entrada: Nx (filas), Ny (columnas), Nt (iteraciones)
// Retorna:
//   0 si la ejecución fue correcta
//   1 si hay error en los argumentos
//-------------------------------------------------------------
int main(int argc, char* argv[]) {
   // Verifica argumentos de entrada
   if(argc != 4){
    std::cerr <<"Debe escribir tres argumentos numerode filas numero de columnas numero de iteraciones" << std::endl;
    exit(1);
  } 
   
  // Dimensiones físicas del dominio
  double Lx=1.0;
  double Ly=1.0;
  int Nx; //Numero de filas
  int Ny; //Numero de columnas
  int Nt; //Numero de iteraciones (tiempo)
  
  // Conversión segura de argumentos
  try{
    Nx=std::stoi(argv[1]);
    Ny=std::stoi(argv[2]);
    Nt = std::stoi(argv[3]);
  }catch(const std::out_of_range& e) {
    std::cerr << "Error: el número está fuera del rango de int." << std::endl;
  } 
  
  // Cálculo de pasos espaciales y temporales
  double dx= Lx/Nx;
  double dy= Ly/Ny;
  double alpha=0.01; // coeficiente de difusión
  double dt=(0.25*std::min(dx, dy)*std::min(dx, dy))/alpha; // paso de tiempo estable
  
  // Reserva dinámica de memoria
  double *current= new double[Nx*Ny];
  double *prox=new double[Nx*Ny];
  double *temp;
  
  // Inicializa la matriz con ceros
  for(int i = 0; i < Nx; ++i) {
      for(int j = 0; j < Ny; ++j) {
          current[(i * Ny) + j] = 0;
      }
  }

  // Define el tamaño del cuadrado caliente en el centro
  double centro_medios= (Nx/10)/2;
  
  // Mide el tiempo de simulación
  double time_1 = seconds();
  
  // Bucle principal de tiempo (simulación)
  for(int n = 0; n < Nt; ++n){
    llenar_centro(Nx, Ny, centro_medios, current); // mantiene el centro caliente
    recalcular_matriz(Nx, Ny, current, prox, alpha, dx, dy, dt); // aplica difusión
    
    // Intercambia punteros (sin copiar matrices)
    temp=current;
    current=prox;
    prox=temp;	
  }
  double time_2 = seconds();
  std::cout << "# Time: " << time_2 - time_1 << std::endl;
  
  // Guarda los resultados en archivo
  guardar_datos( Nx, Ny, current);
  
  // Libera la memoria dinámica
  delete[] current;
  delete[] prox;
  return 0;
}
